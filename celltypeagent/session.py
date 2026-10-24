"""Annotation session with an append-only version history.

A version is fully described by (sample, region state, constraints, corrections). Corrections are cells the
user labeled by hand: they are locked to that label and also retrain the expression model, so every other
cell is re-predicted (the model iterates). Expression-model outputs are cached per set of corrections.
Session directory (work/sessions/<session_id>/):
  session.json        sample id, creation time
  versions/vNNNN.json one file per version: parent, action, message, regions, constraints, corrections
  base/<key>.npy      expression-model log-probabilities for one set of corrections ("init" = none)
  selections/         areas the user selected on the map (cell ids)
  verdicts.jsonl      reviewer verdicts on conflict cells; suggestions.jsonl rule suggestions
  figures/, exports/
"""
import fcntl, hashlib, json, threading, time, uuid
from datetime import datetime
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy.spatial import cKDTree
from . import config, vocab, markers, regions as regmod, constraints as cmod, conflicts as confmod, render


def _tag(s):
    return hashlib.md5(str(s).encode()).hexdigest()[:8]


def _now():
    return datetime.now().isoformat(timespec="seconds")


def _clean(v):
    return None if v in ("nan", "", None) else v


def panel_genes(sample_id=None):
    """Genes of a sample's panel (column order of its counts.npz). Without a sample: the first imported sample's
    panel (all samples of a dataset normally share one panel)."""
    if sample_id is None and config.SAMPLES.exists():
        sample_id = next((d.name for d in sorted(config.SAMPLES.iterdir()) if (d / "genes.txt").exists()), None)
    if sample_id is not None and (config.SAMPLES / str(sample_id) / "genes.txt").exists():
        return (config.SAMPLES / str(sample_id) / "genes.txt").read_text().split()
    return pd.read_csv(config.HOME / "genes.csv", header=None)[0].tolist()  # lab-server bundles


def list_samples():
    rows = []
    for d in sorted(config.SAMPLES.iterdir()) if config.SAMPLES.exists() else []:
        if (d / "meta.json").exists():
            m = json.load(open(d / "meta.json"))
            rows.append({k: _clean(m.get(k)) for k in ["sample_id", "description", "disease", "location", "inflammation",
                                                         "ulceration", "n_cells", "has_pathologist_regions"]}
                        | {"has_he": m.get("he_path") is not None})
    return pd.DataFrame(rows)


def list_sessions():
    rows = []
    if config.SESSIONS.exists():
        for d in sorted(config.SESSIONS.iterdir()):
            if (d / "session.json").exists():
                s = json.load(open(d / "session.json"))
                rows.append({k: s[k] for k in ("session_id", "sample_id", "created_at")}
                            | {"method": s.get("method", "reference"), "n_versions": len(list((d / "versions").glob("v*.json")))})
    return pd.DataFrame(rows)


class Session:
    def __init__(self, session_id):
        self.id = session_id
        self.dir = config.SESSIONS / session_id
        info = json.load(open(self.dir / "session.json"))
        self.sample_id = info["sample_id"]
        sd = config.SAMPLES / self.sample_id
        self.meta = json.load(open(sd / "meta.json"))
        ds_path = config.SAMPLES / "dataset.json"  # reference taxonomy (lab server); marker sessions do not need it
        ds = json.load(open(ds_path)) if ds_path.exists() else {}
        self.scheme = info.get("scheme", {"kind": "reference"})
        if self.scheme["kind"] == "markers":
            # cell types, groups and lineages all come from the user's marker tree
            # cell types are the leaves of the marker tree (+ Unassigned); lineages and parent types are groups.
            # Sessions made before this rule ("leaf_types" missing) keep every node as a type.
            self.node_names = self.scheme.get("nodes", self.scheme["types"])
            self.types = self.scheme["types"]
            self.lineage_of_type = self.scheme["top"]
            self.type_groups = {k: [t for t in v if t in self.types]
                                for k, v in markers.descendants(self.scheme["tree"]).items()}
            self.lineages = [n["name"] for n in self.scheme["tree"]] + [markers.UNASSIGNED]
        else:
            self.types = ds["types"]
            self.lineage_of_type = ds["lineage_of_type"]
            self.type_groups = vocab.type_groups(self.types, ds["lineage_of_type"])
            self.lineages = list(vocab.LINEAGE_GROUPS.values())
        self.sample_dir = sd
        self.panel = panel_genes(self.sample_id)
        self.he_um_per_px = float(self.meta.get("he_um_per_px", config.HE_UM_PER_PX))
        if self.scheme["kind"] != "markers":
            self.node_names = self.types
        self.cells = pd.read_parquet(sd / "cells.parquet")
        self._counts = None
        self._cache, self._bases, self._nbr_frac_cache = {}, {}, {}
        self._lock = threading.RLock()
        for sub in ("base", "selections"):
            (self.dir / sub).mkdir(exist_ok=True)

    @property
    def counts(self):
        if self._counts is None:
            self._counts = sp.load_npz(self.sample_dir / "counts.npz").tocsr()
        return self._counts

    # ---------- creation / persistence ----------
    @classmethod
    def create(cls, sample_id, region_source="none", method="markers", markers_text=None, resolution=None,
               carry_from=None):
        """method 'markers' (default): hierarchical marker annotation (Leiden + marker scores) of this sample; the data
        is treated as having no cell type labels. markers_text defaults to config.DEFAULT_MARKERS.
        method 'reference': expression model trained on the dataset's labeled slides (benchmarks only)."""
        sd = config.SAMPLES / str(sample_id)
        if not (sd / "meta.json").exists():
            raise ValueError(f"unknown sample {sample_id!r}; see list_samples()")
        meta = json.load(open(sd / "meta.json"))
        if region_source == "auto":  # older clients
            region_source = "none"
        if region_source == "pathologist" and not meta["has_pathologist_regions"]:
            raise ValueError("this sample has no pathologist annotation")
        scheme = {"kind": "reference"}
        if method == "markers":
            parsed = markers.parse(markers_text or config.DEFAULT_MARKERS.read_text(), panel_genes(sample_id))
            if not parsed["ok"]:
                raise vocab.VocabError("marker list has errors: " + "; ".join(parsed["errors"][:8]))
            scheme = {"kind": "markers", "text": parsed["canonical"], "tree": parsed["tree"], "types": parsed["leaves"],
                      "nodes": parsed["types"], "leaf_types": True,
                      "parent": parsed["parent"], "top": parsed["top"],
                      "resolution": markers.AUTO if resolution in (None, "", markers.AUTO) else float(resolution)}
        elif method != "reference":
            raise ValueError("method must be 'reference' or 'markers'")
        sid = f"{sample_id}_{time.strftime('%Y%m%d-%H%M%S')}_{uuid.uuid4().hex[:4]}"
        d = config.SESSIONS / sid
        (d / "versions").mkdir(parents=True)
        (d / "figures").mkdir(); (d / "exports").mkdir()
        json.dump({"session_id": sid, "sample_id": str(sample_id), "created_at": _now(), "method": method,
                   "scheme": scheme}, open(d / "session.json", "w"), indent=1)
        s = cls(sid)
        state = {"version": 1, "parent": None, "created_at": _now(), "action": "create",
                 "message": f"session created ({method} annotation); "
                            + ("no regions yet" if region_source == "none" else "regions from the pathologist annotation"),
                 "regions": regmod.default_state(region_source), "constraints": [],
                 "corrections": [], "marker_step": {"ops": []}}
        s._write(state)
        if carry_from:
            s._carry_over(carry_from)
        return s

    def _carry_over(self, sid):
        """Copy regions, priors and corrections of another session of the same sample into a new version.
        Items that name cell types or regions this session does not have are dropped and listed in the message."""
        old = Session(sid)
        if old.sample_id != self.sample_id:
            raise ValueError("can only carry over from a session of the same sample")
        ost = old.state()
        regions = ost["regions"]
        names = regmod.names(regions, regmod.assign(self.cells, regions))
        kept, dropped = [], []
        for c in ost["constraints"]:
            try:
                kept.append(cmod.normalize(c, self.types, self.type_groups, names))
            except vocab.VocabError as e:
                dropped.append(f"{c['id']} ({e})")
        corr = [c for c in ost.get("corrections", []) if c["label"] in self.types]
        dropped += [f"{c['id']} (type {c['label']} not in the new list)" for c in ost.get("corrections", []) if c["label"] not in self.types]
        msg = (f"carried over from {sid}: {len(regions.get('ops', []))} region polygons, {len(kept)} priors, "
               f"{len(corr)} corrections" + (f"; dropped: {'; '.join(dropped)}" if dropped else ""))
        return self._commit("carry_over", msg, regions=regions, constraints=kept, corrections=corr)

    def _vpath(self, v):
        return self.dir / "versions" / f"v{v:04d}.json"

    def _write(self, state):
        json.dump(state, open(self._vpath(state["version"]), "w"), indent=1)

    def versions(self):
        return sorted(int(p.stem[1:]) for p in (self.dir / "versions").glob("v*.json"))

    def state(self, version=None):
        v = version or self.versions()[-1]
        return json.load(open(self._vpath(v)))

    def _commit(self, action, message, regions=None, constraints=None, expected_parent=None, corrections=None,
                marker_step=None):
        """Append a version. The lock keeps version numbers unique when the agent and the review page
        (or two reviewers) commit at the same time; expected_parent rejects edits made on a stale version."""
        with open(self.dir / "versions" / ".lock", "w") as lk:
            fcntl.flock(lk, fcntl.LOCK_EX)
            cur = self.state()
            if expected_parent is not None and cur["version"] != expected_parent:
                raise ValueError(f"session moved on to v{cur['version']} (you were on v{expected_parent}); reload and retry")
            new = {"version": cur["version"] + 1, "parent": cur["version"], "created_at": _now(), "action": action,
                   "message": message, "regions": regions if regions is not None else cur["regions"],
                   "constraints": constraints if constraints is not None else cur["constraints"],
                   "corrections": corrections if corrections is not None else cur.get("corrections", []),
                   "marker_step": marker_step if marker_step is not None else cur.get("marker_step", {"ops": []})}
            self._write(new)
        return new["version"]

    # ---------- computation ----------
    @staticmethod
    def _mkey(ops):
        return "init" if not ops else _tag(json.dumps(ops, sort_keys=True))

    def _base_key(self, st):
        corr, ops = st.get("corrections", []), st.get("marker_step", {}).get("ops", [])
        if not corr and not ops:
            return "init"
        return _tag(json.dumps([[[c["label"], sorted(c["cells"])] for c in corr], ops], sort_keys=True))

    def _locks(self, corrections):
        """Rows and type indices of hand-labeled cells; a later correction of the same cell wins."""
        last = {}
        for c in corrections:
            for cell in c["cells"]:
                last[int(cell)] = self.types.index(c["label"])
        if not last:
            return np.array([], int), np.array([], int)
        rows = np.array([self.row_of(c) for c in last])
        return rows, np.array(list(last.values()))

    def base(self, version=None, progress=None):
        """Expression-model output for a version (its marker step and corrections): computed once, cached in memory
        and on disk."""
        st = version if isinstance(version, dict) else self.state(version)
        key = self._base_key(st)
        with self._lock:
            if key not in self._bases:
                f = self.dir / "base" / f"{key}.npy"
                if f.exists():
                    logp = np.load(f)
                else:
                    logp = self._train_base(st, progress)
                    np.save(f, logp)
                pred, conf = cmod.posterior(logp, [], None, self.types)
                self._bases[key] = {"key": key, "logp": logp, "pred": pred, "conf": conf}
            elif progress:
                progress(1.0, "expression model ready")
            return self._bases[key]

    def _init_model(self, st, say):
        """Expression model before hand corrections: trained on the marker-step labels of this version (marker mode)
        or on the reference slides (reference mode)."""
        from . import model
        mkey = self._mkey(st.get("marker_step", {}).get("ops", []))
        f = self.dir / "base" / ("init_model.npz" if mkey == "init" else f"model_{mkey}.npz")
        if f.exists():
            return model.load(f)
        if self.scheme["kind"] == "markers":
            lab = self.marker_labels(st, lambda fr, msg: say(0.02 + 0.6 * fr, msg))
            say(0.65, "training expression model on the marker annotation")
            m = model.train_pseudo(self.counts, lab, len(self.types), progress=lambda fr, msg: say(0.65 + 0.25 * fr, msg))
        else:
            say(0.02, "loading reference cells")
            model.reference()
            m = model.train(self.meta["slide_id"], progress=lambda fr, msg: say(0.05 + 0.85 * fr, msg))
        model.save(m, f)
        return m

    # ----- marker step (versioned: each local re-clustering is an op in state["marker_step"]["ops"]) -----
    def _marker_tree(self, ops):
        """The marker tree with the marker edits of `ops` applied (edits accumulate in op order)."""
        import copy
        tree = copy.deepcopy(self.scheme["tree"])
        for op in ops:
            for t, mk in (op.get("markers") or {}).items():
                node = markers._find(tree, t)
                node["markers"] = markers.parse(f"{t}: {mk}", self.panel)["tree"][0]["markers"]
        return tree

    def _save_mres(self, mkey, res):
        names = sorted(set(res["cluster_of"].tolist()))
        code = {n: i for i, n in enumerate(names)}
        ids = list(res["members"])
        rows = [res["members"][i] for i in ids]
        np.savez(self.dir / "base" / f"markers_{mkey}.npz",
                 labels=np.array([self.node_names.index(x) for x in res["labels"]], np.int16),
                 cluster=np.array([code[c] for c in res["cluster_of"]], np.int32),
                 member_rows=np.concatenate(rows).astype(np.int32) if rows else np.zeros(0, np.int32),
                 member_offsets=np.cumsum([0] + [len(r) for r in rows]))
        json.dump({"cluster_names": names, "member_ids": ids, "records": res["records"], "nodes": res["nodes"]},
                  open(self.dir / "base" / f"markers_{mkey}.json", "w"))

    def _load_mres(self, mkey):
        f = self.dir / "base" / f"markers_{mkey}.npz"
        if not f.exists():
            return None
        z, j = np.load(f), json.load(open(self.dir / "base" / f"markers_{mkey}.json"))
        off = z["member_offsets"]
        return {"labels": np.array(self.node_names, dtype=object)[z["labels"]],
                "cluster_of": np.array(j["cluster_names"], dtype=object)[z["cluster"]],
                "members": {i: z["member_rows"][off[k]:off[k + 1]] for k, i in enumerate(j["member_ids"])},
                "records": j["records"], "nodes": j["nodes"]}

    def marker_result(self, st=None, progress=None):
        """Marker-step result of a version: the full annotation, then each re-clustering op applied in order.
        Every prefix of the op list is cached, so a new op only re-runs its own part."""
        st = st if isinstance(st, dict) else self.state(st)
        ops = st.get("marker_step", {}).get("ops", [])
        with self._lock:
            if not hasattr(self, "_mres"):
                self._mres = {}
            res = None
            for i in range(len(ops), -1, -1):  # longest cached prefix
                key = self._mkey(ops[:i])
                res = self._mres.get(key) or self._load_mres(key)
                if res is not None:
                    self._mres[key] = res
                    break
            if res is None:
                res = markers.annotate(self.counts, self.panel, self.scheme["tree"], progress,
                                       resolution=self.scheme.get("resolution", markers.RESOLUTION))
                self._save_mres("init", res)
                self._mres["init"] = res
                i = 0
            for j in range(i, len(ops)):
                res = self._apply_recluster(res, ops[j], ops[:j + 1], progress)
                key = self._mkey(ops[:j + 1])
                self._save_mres(key, res)
                self._mres[key] = res
            return res

    def _apply_recluster(self, res, op, ops_so_far, progress=None):
        return markers.rerun(res, self.counts, self.panel, self._marker_tree(ops_so_far), op["target"],
                             op["resolution"], prefix=f"{op['id']}:", progress=progress)

    def marker_labels(self, st=None, progress=None):
        """Marker-step label index per cell for a version."""
        idx = {t: i for i, t in enumerate(self.types)}  # -1: left at a lineage / parent level, not a cell type
        return np.array([idx.get(x, -1) for x in self.marker_result(st, progress)["labels"]], np.int64)

    def clusters(self, version=None):
        """Marker-step clusters with their assignment, marker scores and current label composition."""
        if self.scheme["kind"] != "markers":
            return []
        st = self.state(version)
        mres = self.marker_result(st)
        lab = np.array(self.types)[self.result(version)["pred"]]
        out = []
        for r in mres["records"]:
            rows = mres["members"].get(r["id"], np.zeros(0, int))
            out.append(dict(r, n_cells_now=int(len(rows)),
                            labels_now=pd.Series(lab[rows]).value_counts().head(4).to_dict() if len(rows) else {}))
        return out

    def marker_levels(self, version=None):
        """Leiden resolution chosen at every level of the marker step (with the candidates tried)."""
        return self.marker_result(self.state(version))["nodes"] if self.scheme["kind"] == "markers" else []

    def _within_parent(self, logp, st):
        """Cells the marker step left at a lineage or parent level get a cell type below that level only: the other
        types are down-weighted by 1e-4, so the expression model picks among the right subtypes."""
        lab = self.marker_result(st)["labels"]
        groups = markers.descendants(self.scheme["tree"])
        lp = logp.astype(np.float32)
        for node in set(lab.tolist()) - set(self.types):
            rows = np.where(lab == node)[0]
            allowed = np.isin(self.types, groups.get(node, []))
            lp[np.ix_(rows, np.where(~allowed)[0])] += np.float32(np.log(1e-4))
        lp -= lp.max(1, keepdims=True)
        return (lp - np.log(np.exp(lp).sum(1, keepdims=True))).astype(np.float16)

    def _train_base(self, st, progress=None):
        """Model before corrections (see _init_model); hand corrections then fine-tune it."""
        say = progress or (lambda f, m: None)
        from . import model
        rows, labs = self._locks(st.get("corrections", []))
        m = self._init_model(st, say if not len(rows) else (lambda f, msg: say(0.3 * f, msg)))
        if len(rows):
            if self.scheme["kind"] == "markers":
                rX, ry, _ = model.pseudo_pool(self.counts, self.marker_labels(st), len(self.types))
            else:
                rX, ry = model.reference_replay(self.meta["slide_id"])
            m = model.finetune(m, self.counts[rows], labs, rX, ry, progress=lambda fr, msg: say(0.3 + 0.6 * fr, msg))
        say(0.92, "predicting all cells")
        logp = model.predict(m, self.counts)
        if self.scheme.get("leaf_types"):
            logp = self._within_parent(logp, st)
        say(1.0, "expression model ready")
        return logp

    def knn(self):
        """Row indices of every cell's K nearest neighbors (itself excluded)."""
        if not hasattr(self, "_nn"):
            xy = self.cells[["x", "y"]].to_numpy()
            self._nn = cKDTree(xy).query(xy, k=cmod.NEIGHBORS_K + 1)[1][:, 1:]
        return self._nn

    def _neighbor_frac_fn(self, base):
        def frac(type_names):
            """Fraction of each cell's K nearest neighbors whose expression-model label is in type_names."""
            key = (base["key"], tuple(sorted(type_names)))
            if key not in self._nbr_frac_cache:
                hit = np.isin(base["pred"], [self.types.index(t) for t in type_names])
                self._nbr_frac_cache[key] = hit[self.knn()].mean(1)
            return self._nbr_frac_cache[key]
        return frac

    def _ctx(self, region, base):
        return {"region": region, "depth_um": self.cells.depth_um.to_numpy(), "neighbor_frac": self._neighbor_frac_fn(base)}

    def _compute(self, st, constraints=None):
        base = self.base(st)
        constraints = st["constraints"] if constraints is None else constraints
        region = regmod.assign(self.cells, st["regions"])
        ctx = self._ctx(region, base)
        pred, conf, top_i, top_p = cmod.posterior(base["logp"], constraints, ctx, self.types, top=3)
        rows, labs = self._locks(st.get("corrections", []))
        if len(rows):
            pred[rows], conf[rows] = labs, 1.0
            top_i[rows, 0], top_p[rows, 0] = labs, 1.0
        return {"region": region, "pred": pred, "conf": conf, "top_idx": top_i, "top_p": top_p, "ctx": ctx,
                "base": base, "locked_rows": rows}

    def result(self, version=None):
        with self._lock:
            st = self.state(version)
            v = st["version"]
            if v not in self._cache:
                c = self._compute(st)
                b = c["base"]
                cf = confmod.find(self.cells, c["region"], b["pred"], b["conf"], c["pred"], st["constraints"],
                                  self.types, c["ctx"])
                cf["cell"] = self.cells.cell.to_numpy()[cf.row.to_numpy(dtype=int)] if len(cf) else []
                self._cache[v] = {"state": st, "region": c["region"], "pred": c["pred"], "conf": c["conf"],
                                  "top_idx": c["top_idx"], "top_p": c["top_p"], "base_pred": b["pred"],
                                  "base_conf": b["conf"], "locked_rows": c["locked_rows"], "conflicts": cf}
            return self._cache[v]

    def lineage_labels(self, version=None):
        """Lineage (top-level type) of every cell's current label: the dataset's own lineage split is not used."""
        lut = np.array([self.lineage_of_type.get(t, t) for t in self.types], dtype=object)
        return lut[self.result(version)["pred"]]

    def lineage_colors(self):
        return render.lineage_palette(self.lineages)

    def labels(self, version=None):
        r = self.result(version)
        return np.array(self.types)[r["pred"]]

    # ---------- queries ----------
    def vocabulary(self):
        return {"method": self.scheme["kind"], "cell_types": self.types,
                "type_groups": {k: v for k, v in self.type_groups.items()},
                "lineages": self.lineages, "lineage_of_type": self.lineage_of_type,
                "marker_list": self.scheme.get("text"),
                "regions": self.region_names(), "free_regions": vocab.FREE_REGIONS, "kinds": cmod.KINDS,
                "modes": cmod.MODES, "strengths": cmod.STRENGTH, "verdicts": vocab.VERDICTS,
                "depth_note": "depth = distance to the mucosa core (epithelial-marker density); 0 inside crypts",
                "region_note": "regions are drawn and named by the user; cells outside every region are 'unknown'"}

    def region_names(self, version=None):
        reg = self.state(version)["regions"]
        return regmod.names(reg, regmod.assign(self.cells, reg))

    def summary(self, version=None):
        r = self.result(version)
        st = r["state"]
        lab = pd.Series(np.array(self.types)[r["pred"]])
        return {"session_id": self.id, "sample": {k: _clean(v) for k, v in self.meta.items() if k != "he_path"},
                "version": st["version"], "latest_version": self.versions()[-1], "message": st["message"],
                "region_source": st["regions"]["source"],
                "region_composition": pd.Series(r["region"]).value_counts().to_dict(),
                "constraints": [{"id": c["id"], "kind": c.get("kind", "region"), "text": c["text"], "author": c.get("author"),
                                 "rule": cmod.describe(c)} for c in st["constraints"]],
                "corrections": [{"id": c["id"], "label": c["label"], "n_cells": len(c["cells"]), "text": c.get("text", ""),
                                 "author": c.get("author")} for c in st.get("corrections", [])],
                "n_changed_vs_expression_only": int((r["pred"] != r["base_pred"]).sum()),
                "n_changed_since_first": int((r["pred"] != self.result(self.versions()[0])["pred"]).sum()),
                "has_regions": bool(self.region_names(version)),
                "n_conflicts": int(len(r["conflicts"])),
                "top_labels": lab.value_counts().head(12).to_dict()}

    def diff(self, a, b):
        ra, rb = self.result(a), self.result(b)
        ch = ra["pred"] != rb["pred"]
        t = np.array(self.types)
        tr = pd.DataFrame({"from": t[ra["pred"][ch]], "to": t[rb["pred"][ch]], "region": rb["region"][ch]})
        reg_changed = int((ra["region"] != rb["region"]).sum())
        return {"from_version": a, "to_version": b, "n_label_changes": int(ch.sum()), "n_region_changes": reg_changed,
                "top_transitions": tr.value_counts().head(15).reset_index().to_dict("records")}

    def history(self):
        rows = []
        for v in self.versions():
            st = self.state(v)
            rows.append({"version": v, "parent": st["parent"], "action": st["action"], "message": st["message"],
                         "n_constraints": len(st["constraints"]), "n_corrections": len(st.get("corrections", [])),
                         "created_at": st["created_at"]})
        return rows

    def conflict_summary(self, by=("constraint", "region"), version=None):
        return confmod.summarize(self.result(version)["conflicts"], by)

    # ---------- constraint actions ----------
    def _next_id(self, constraints):
        used = {c["id"] for c in constraints}
        i = 1
        while f"C{i}" in used:
            i += 1
        return f"C{i}"

    def _check_depth(self, c):
        if c.get("kind") == "depth" and not np.isfinite(self.cells.depth_um.to_numpy()).any():
            raise vocab.VocabError("depth priors need an epithelial-marker depth estimate, which this sample does not have "
                                   "(its panel lacks EPCAM / KRT8); use a region prior instead")
        return c

    def preview_constraint(self, spec):
        """Validate and show the effect of a new constraint without committing it."""
        cur = self.state()
        c = self._check_depth(cmod.normalize(dict(spec, id=spec.get("id") or self._next_id(cur["constraints"])), self.types, self.type_groups, self.region_names()))
        before = self.result()
        after = self._compute(cur, cur["constraints"] + [c])
        region, pred, ctx = after["region"], after["pred"], after["ctx"]
        ch = pred != before["pred"]
        t = np.array(self.types)
        tr = pd.DataFrame({"from": t[before["pred"][ch]], "to": t[pred[ch]], "region": region[ch]})
        affected = np.isin(before["pred"], [self.types.index(x) for x in c["types_resolved"]])
        return {"interpretation": cmod.describe(c), "based_on_version": cur["version"],
                "constraint": {k: c[k] for k in cmod.SPEC_FIELDS if k in c},
                "types_resolved": c["types_resolved"],
                "condition_true_cells": int(cmod.condition_mask(c, ctx).sum()),
                "n_cells_currently_of_these_types": int(affected.sum()),
                "n_cells_of_these_types_penalized": int((affected & cmod.penalized_cells(c, ctx)).sum()),
                "n_label_changes": int(ch.sum()), "top_transitions": tr.value_counts().head(10).reset_index().to_dict("records")}

    def add_constraint(self, spec, message=None, expected_parent=None):
        cur = self.state()
        c = self._check_depth(cmod.normalize(dict(spec, id=spec.get("id") or self._next_id(cur["constraints"])), self.types, self.type_groups, self.region_names()))
        if c["id"] in {x["id"] for x in cur["constraints"]}:
            raise ValueError(f"constraint id {c['id']} already exists; use update_constraint")
        by = f" by {c['author']}" if c.get("author") else ""
        return self._commit("add_constraint", message or f"add {c['id']}{by}: {c['text']}",
                            constraints=cur["constraints"] + [c], expected_parent=expected_parent or cur["version"])

    def update_constraint(self, cid, changes, message=None, expected_parent=None):
        cur = self.state()
        out, found = [], False
        for c in cur["constraints"]:
            if c["id"] == cid:
                found = True
                c = cmod.normalize(c | changes, self.types, self.type_groups, self.region_names())
            out.append(c)
        if not found:
            raise ValueError(f"no constraint {cid}")
        return self._commit("update_constraint", message or f"update {cid}: {changes}", constraints=out,
                            expected_parent=expected_parent or cur["version"])

    def remove_constraint(self, cid, message=None, expected_parent=None):
        cur = self.state()
        out = [c for c in cur["constraints"] if c["id"] != cid]
        if len(out) == len(cur["constraints"]):
            raise ValueError(f"no constraint {cid}")
        return self._commit("remove_constraint", message or f"remove {cid}", constraints=out,
                            expected_parent=expected_parent or cur["version"])

    # ---------- region actions (regions are drawn by the user) ----------
    def set_region_source(self, source, message=None):
        """Start from no regions ('none') or the dataset's pathologist annotation; drawn regions are kept."""
        if source == "pathologist" and not self.meta["has_pathologist_regions"]:
            raise ValueError("this sample has no pathologist annotation")
        cur = self.state()
        reg = dict(regmod.default_state(source), ops=cur["regions"].get("ops", []))
        return self._commit("set_region_source", message or f"start regions from {source}", regions=reg,
                            expected_parent=cur["version"])

    def add_region(self, name, polygon_px, author=None, expected_parent=None):
        """Save a polygon drawn on the map (H&E level-0 pixels) as (part of) region `name`."""
        name = regmod.validate_name(name)
        rows = regmod.cells_in_polygon(self.cells.px.to_numpy(), self.cells.py.to_numpy(), polygon_px)
        if not len(rows):
            raise vocab.VocabError("the polygon contains no cells")
        cur = self.state()
        ops = cur["regions"].get("ops", [])
        i = 1
        while f"R{i}" in {o.get("id") for o in ops}:
            i += 1
        op = {"op": "polygon", "id": f"R{i}", "name": name, "polygon": [[float(x), float(y)] for x, y in polygon_px],
              "cells": self.cells.cell.to_numpy()[rows].tolist(), "author": author, "created_at": _now()}
        reg = dict(cur["regions"], ops=ops + [op])
        return self._commit("add_region", f"R{i}: {len(rows)} cells drawn as region {name}", regions=reg,
                            expected_parent=expected_parent or cur["version"])

    def rename_region(self, old, new, message=None):
        """Rename a region everywhere (all its polygons and every constraint that names it)."""
        new = regmod.validate_name(new)
        cur = self.state()
        ops = [dict(o, name=new) if o.get("name") == old else o for o in cur["regions"].get("ops", [])]
        if ops == cur["regions"].get("ops", []):
            raise ValueError(f"no drawn region named {old!r}")
        cons = []
        for c in cur["constraints"]:
            if old in c.get("regions_resolved", []) or old in c.get("regions", []):
                c = dict(c, regions=[new if x == old else x for x in c.get("regions", [])],
                         regions_resolved=[new if x == old else x for x in c.get("regions_resolved", [])])
            cons.append(c)
        return self._commit("rename_region", message or f"rename region {old} -> {new}",
                            regions=dict(cur["regions"], ops=ops), constraints=cons, expected_parent=cur["version"])

    def region_users(self, rid):
        """Priors that would lose their region if polygon rid were deleted (it is the region's last polygon)."""
        cur = self.state()
        ops = cur["regions"].get("ops", [])
        op = next((o for o in ops if o.get("id") == rid), None)
        if op is None:
            raise ValueError(f"no region polygon {rid}")
        still = any(o.get("name") == op["name"] for o in ops if o.get("id") != rid) or op["name"] in regmod.assign(
            self.cells, dict(cur["regions"], ops=[])).tolist()
        users = [] if still else [c for c in cur["constraints"] if op["name"] in c.get("regions_resolved", [])]
        return op, users

    def delete_region(self, rid, message=None, remove_priors=False):
        """Delete one drawn polygon. If it is the last polygon of a region that priors use, it is refused unless
        remove_priors, which removes those priors in the same version (Undo restores both)."""
        cur = self.state()
        op, users = self.region_users(rid)
        if users and not remove_priors:
            raise vocab.VocabError(f"region {op['name']} is used by {', '.join(c['id'] for c in users)}; "
                                   "change or remove those priors first")
        drop = {c["id"] for c in users}
        msg = f"delete {rid} (region {op['name']})" + (f" and priors {', '.join(sorted(drop))}" if drop else "")
        return self._commit("delete_region", message or msg,
                            regions=dict(cur["regions"], ops=[o for o in cur["regions"]["ops"] if o.get("id") != rid]),
                            constraints=[c for c in cur["constraints"] if c["id"] not in drop],
                            expected_parent=cur["version"])

    def regions_info(self, version=None):
        """Drawn polygons and per-region cell counts of a version."""
        r = self.result(version)
        counts = pd.Series(r["region"]).value_counts().to_dict()
        polys = [{k: o[k] for k in ("id", "name", "polygon", "author", "created_at") if k in o} | {"n_cells_drawn": len(o["cells"])}
                 for o in r["state"]["regions"].get("ops", []) if o["op"] == "polygon"]
        return {"source": r["state"]["regions"].get("source", "none"), "polygons": polys,
                "regions": [{"name": n, "n_cells": int(counts.get(n, 0))} for n in self.region_names(version)],
                "n_unknown": int(counts.get("unknown", 0))}

    def revert(self, version, message=None):
        cur = self.state()
        old = self.state(version)
        return self._commit("revert", message or f"revert to v{version}", regions=old["regions"],
                            constraints=old["constraints"], corrections=old.get("corrections", []),
                            marker_step=old.get("marker_step", {"ops": []}), expected_parent=cur["version"])

    # ---------- hand corrections ----------
    def _cells_for(self, region, cells, only_types, version=None):
        r = self.result(version)
        if region:
            if region not in self.region_names(version):
                raise vocab.VocabError(f"no region named {region!r}; regions: {self.region_names(version)}")
            cells = self.cells.cell.to_numpy()[r["region"] == region].tolist()
        cells = [int(c) for c in cells or []]
        if only_types:
            keep = vocab.resolve(only_types, self.types, self.type_groups, "cell type")
            lab = np.array(self.types)[r["pred"][[self.row_of(c) for c in cells]]]
            cells = [c for c, l in zip(cells, lab) if l in keep]
        if not cells:
            raise ValueError("no cells match")
        return cells

    def label_cells(self, cell_type, region=None, cells=None, only_types=None, text="", author=None):
        """Hand-label cells. They are locked to cell_type and retrain the expression model for all other cells."""
        if self.scheme["kind"] == "markers" and cell_type in self.type_groups and cell_type not in self.types:
            raise vocab.VocabError(f"{cell_type} is a lineage / parent type; label the cells with one of its cell types: "
                                   + ", ".join(self.type_groups[cell_type]))
        label = vocab.resolve([cell_type], self.types, {}, "cell type")[0]
        cur = self.state()
        cells = self._cells_for(region, cells, only_types)
        corr = cur.get("corrections", [])
        i = 1
        while f"L{i}" in {c["id"] for c in corr}:
            i += 1
        new = {"id": f"L{i}", "cells": cells, "label": label, "text": text, "author": author, "region": region,
               "only_types": only_types}
        return self._commit("label_cells", f"L{i}: {len(cells)} cells labeled {label}" + (f" ({text})" if text else ""),
                            corrections=corr + [new], expected_parent=cur["version"])

    def remove_correction(self, cid, message=None):
        cur = self.state()
        corr = [c for c in cur.get("corrections", []) if c["id"] != cid]
        if len(corr) == len(cur.get("corrections", [])):
            raise ValueError(f"no correction {cid}")
        return self._commit("remove_correction", message or f"remove {cid}", corrections=corr, expected_parent=cur["version"])

    # ---------- local re-clustering of the marker step ----------
    def _recluster_op(self, target, resolution, marker_edits, text, author, cur):
        if self.scheme["kind"] != "markers":
            raise ValueError("re-clustering needs a marker-based annotation")
        if resolution not in (markers.AUTO, None):
            resolution = float(resolution)
            if not 0.1 <= resolution <= 10:
                raise vocab.VocabError("resolution must be 'auto' or a number between 0.1 and 10")
        resolution = resolution if resolution is not None else markers.AUTO
        ops = cur.get("marker_step", {}).get("ops", [])
        mres = self.marker_result(cur)
        valid = {"root"} | {r["id"] for r in mres["records"]} | set(self.node_names)
        if target not in valid:
            raise vocab.VocabError(f"unknown target {target!r}: use 'root', a cell type, or a cluster id from list_clusters")
        edits = {}
        for t, mk in (marker_edits or {}).items():
            t = vocab.resolve([t], self.node_names, {}, "cell type")[0]
            p = markers.parse(f"{t}: {mk}", self.panel)
            if not p["ok"]:
                raise vocab.VocabError("; ".join(p["errors"]))
            edits[t] = markers.to_text(p["tree"]).split(":", 1)[1].strip()
        i = 1
        while f"K{i}" in {o["id"] for o in ops}:
            i += 1
        return {"op": "recluster", "id": f"K{i}", "target": target, "resolution": resolution, "markers": edits,
                "text": text, "author": author}

    @staticmethod
    def _op_sig(op):
        return json.dumps({k: op[k] for k in ("id", "target", "resolution", "markers")}, sort_keys=True)

    def preview_recluster(self, target, resolution=markers.AUTO, marker_edits=None):
        """Re-cluster one level ('root' or a type) or one cluster without committing; report what would change in the
        marker step (the expression model is retrained only on commit)."""
        cur = self.state()
        op = self._recluster_op(target, resolution, marker_edits, "", None, cur)
        ops = cur.get("marker_step", {}).get("ops", [])
        before = self.marker_result(cur)
        after = self._apply_recluster(before, op, ops + [op])
        # the commit reuses exactly this result: Leiden on a GPU-built graph is not bit-reproducible, and the
        # automatic resolution choice can flip between near-equal candidates
        self._previews = {(cur["version"], self._op_sig(op)): after}
        new_recs = [r for r in after["records"] if r["id"].startswith(op["id"] + ":")]
        rows = np.unique(np.concatenate([after["members"][r["id"]] for r in new_recs])) if new_recs else np.zeros(0, int)
        tr = pd.DataFrame({"from": before["labels"][rows], "to": after["labels"][rows]})
        tr = tr[tr["from"] != tr["to"]]
        levels = [n for n in after["nodes"] if n.get("prefix") == op["id"] + ":"]
        return {"based_on_version": cur["version"], "target": target, "n_cells": int(len(rows)),
                "n_marker_labels_changed": int(len(tr)),
                "transitions": tr.value_counts().head(10).reset_index().to_dict("records"),
                "levels": [{k: n[k] for k in ("node", "n_cells", "resolution", "tried", "assigned_frac")} for n in levels],
                "new_clusters": [{k: r[k] for k in ("id", "parent", "n_cells", "assigned", "margin")}
                                 for r in sorted(new_recs, key=lambda r: -r["n_cells"])[:12]],
                "marker_edits": op["markers"]}

    def recluster(self, target, resolution=markers.AUTO, marker_edits=None, text="", author=None):
        """Commit a local re-clustering. Priors, regions and corrections are kept; the expression model is retrained
        on the new marker-step labels (and corrections re-applied) the next time the version is computed."""
        cur = self.state()
        op = self._recluster_op(target, resolution, marker_edits, text, author, cur)
        ms = {"ops": cur.get("marker_step", {}).get("ops", []) + [op]}
        previewed = getattr(self, "_previews", {}).get((cur["version"], self._op_sig(op)))
        if previewed is not None:
            with self._lock:
                key = self._mkey(ms["ops"])
                self._save_mres(key, previewed)
                if not hasattr(self, "_mres"):
                    self._mres = {}
                self._mres[key] = previewed
        res = "auto" if op["resolution"] == markers.AUTO else f"{op['resolution']:g}"
        return self._commit("recluster", f"{op['id']}: re-clustered {target} (resolution {res})"
                            + (f", markers edited for {', '.join(op['markers'])}" if op["markers"] else "")
                            + (f" ({text})" if text else ""), marker_step=ms, expected_parent=cur["version"])

    def cell_info(self, row, version=None):
        r = self.result(version)
        t = np.array(self.types)
        c = self.cells.iloc[row]
        hand = bool(np.isin(row, r["locked_rows"]))
        if hand:  # show what the expression model itself thinks of a hand-labeled cell
            lp = self.base(r["state"])["logp"][row].astype(float)
            idx = np.argsort(-lp)[:3]
            top = [{"type": t[i], "p": float(np.exp(lp[i]))} for i in idx]
        else:
            top = [{"type": t[i], "p": float(p)} for i, p in zip(r["top_idx"][row], r["top_p"][row])]
        return {"cell": int(c.cell), "label": str(t[r["pred"][row]]), "confidence": float(r["conf"][row]), "top": top,
                "expression_label": str(t[r["base_pred"][row]]), "expression_conf": float(r["base_conf"][row]),
                "hand_labeled": hand, "region": str(r["region"][row]),
                "depth_um": round(float(c.depth_um), 1),
                "total_counts": int(c.total_counts), "x_um": round(float(c.x), 1), "y_um": round(float(c.y), 1)} \
            | self._cluster_info(row, r["state"])

    def _cluster_info(self, row, st):
        if self.scheme["kind"] != "markers":
            return {}
        mres = self.marker_result(st)
        cid = str(mres["cluster_of"][row])
        r = next((x for x in mres["records"] if x["id"] == cid), {})
        return {"cluster": cid, "marker_label": str(mres["labels"][row]), "cluster_scores": r.get("scores", {}),
                "cluster_assigned": r.get("assigned"), "cluster_resolution": r.get("resolution")}

    # ---------- rule suggestions (free text from reviewers, turned into constraints by the agent) ----------
    def add_suggestion(self, text, author, cell=None, group=None, version=None):
        if not text.strip():
            raise ValueError("suggestion text is empty")
        sid = f"S{uuid.uuid4().hex[:6]}"
        rec = {"id": sid, "event": "created", "text": text.strip(), "author": author, "cell": cell, "group": group,
               "version": version or self.versions()[-1], "time": _now()}
        self._append("suggestions.jsonl", rec)
        return sid

    def resolve_suggestion(self, sid, status, note="", resolved_by="agent"):
        if status not in ("accepted", "rejected", "open"):
            raise ValueError("status must be accepted, rejected or open")
        if sid not in {s["id"] for s in self.suggestions()}:
            raise ValueError(f"no suggestion {sid}")
        self._append("suggestions.jsonl", {"id": sid, "event": status, "note": note, "author": resolved_by, "time": _now()})

    def suggestions(self, status=None):
        p = self.dir / "suggestions.jsonl"
        if not p.exists():
            return []
        out = {}
        for line in open(p):
            e = json.loads(line)
            if e["event"] == "created":
                out[e["id"]] = e | {"status": "open", "resolution_note": ""}
            elif e["id"] in out:
                out[e["id"]].update(status=e["event"], resolution_note=e.get("note", ""), resolved_by=e["author"])
        rows = list(out.values())
        return [r for r in rows if status is None or r["status"] == status]

    def _append(self, name, rec):
        with open(self.dir / name, "a") as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            f.write(json.dumps(rec) + "\n")
            fcntl.flock(f, fcntl.LOCK_UN)

    # ---------- review ----------
    def review_list(self, gid, version=None):
        """All conflict cells of a group in a fixed random order, identical for every reviewer."""
        sel = confmod.select(self.result(version)["conflicts"], gid)
        return sel.sample(frac=1, random_state=0) if len(sel) else sel

    def review_queue(self, gid, n=8, version=None):
        sel = self.review_list(gid, version).head(n)
        done = self.verdicts()
        reviewed = set(done.cell) if len(done) else set()
        return [r | {"already_reviewed": r["cell"] in reviewed} for r in sel.to_dict("records")]

    def row_of(self, cell):
        if not hasattr(self, "_row_of"):
            self._row_of = pd.Series(np.arange(len(self.cells)), index=self.cells.cell.to_numpy())
        return int(self._row_of[int(cell)])

    def record_verdicts(self, cells, verdict, note="", reviewer="agent", version=None, group=None):
        if verdict not in vocab.VERDICTS:
            raise vocab.VocabError(f"verdict must be one of {list(vocab.VERDICTS)}")
        v = version or self.versions()[-1]
        with open(self.dir / "verdicts.jsonl", "a") as f:
            fcntl.flock(f, fcntl.LOCK_EX)  # several reviewers may write at once
            for c in cells:
                f.write(json.dumps({"cell": int(c), "verdict": verdict, "note": note, "reviewer": reviewer,
                                    "version": v, "group": group, "time": _now()}) + "\n")
            fcntl.flock(f, fcntl.LOCK_UN)
        return len(cells)

    def verdicts(self, latest=False):
        """All verdict records; latest=True keeps only the last verdict per (cell, reviewer)."""
        p = self.dir / "verdicts.jsonl"
        if not p.exists() or p.stat().st_size == 0:
            return pd.DataFrame(columns=["cell", "verdict", "note", "reviewer", "version", "group", "time"])
        v = pd.read_json(p, lines=True)
        return v.drop_duplicates(["cell", "reviewer"], keep="last") if latest else v

    # ---------- figures ----------
    def fig_region_map(self, version=None, conflict_group=None):
        r = self.result(version)
        hl = None
        if conflict_group is not None:
            hl = confmod.select(r["conflicts"], conflict_group).row.to_numpy(dtype=int) if conflict_group else r["conflicts"].row.to_numpy(dtype=int)
        out = self.dir / "figures" / f"map_v{r['state']['version']}_{_tag(conflict_group)}.png"
        return render.region_map(self.cells, r["region"], out, title=f"sample {self.sample_id} v{r['state']['version']} "
                                 f"regions ({r['state']['regions']['source']})", highlight=hl,
                                 highlight_label=conflict_group or "conflicts")

    def fig_conflicts(self, gid, n=6, version=None):
        q = self.review_queue(gid, n, version)
        if not q:
            raise ValueError(f"no conflicts in group {gid!r}")
        labels = self.labels(version)
        out = self.dir / "figures" / f"conflicts_{_tag(gid)}.jpg"
        return render.conflict_sheet(self.cells, q, self.meta["he_path"], labels, out, title=gid,
                                     lineage=self.lineage_labels(version), palette=self.lineage_colors(),
                                     um_per_px=self.he_um_per_px), q

    # ---------- export ----------
    def cell_ids(self):
        """Original Xenium cell ids (obs names of the source h5ad files), one per row of self.cells."""
        if "cell_id" in self.cells:
            return self.cells.cell_id.to_numpy()
        if not hasattr(self, "_cell_ids"):  # lab-server bundles: ids from the merged obs table
            import pyarrow.parquet as pq
            names = pq.read_table(config.WORK / "obs_all.parquet", columns=["NAME"]).column(0).to_numpy()
            self._cell_ids = names[self.cells.cell.to_numpy()]
        return self._cell_ids

    def type_path(self, t):
        """Types from the top level down to t (marker tree), or [lineage, t] for the reference taxonomy."""
        if self.scheme["kind"] == "markers":
            path = [t]
            while self.scheme["parent"].get(path[0]):
                path.insert(0, self.scheme["parent"][path[0]])
            return path
        return [self.lineage_of_type.get(t, t), t]

    def export_table(self, version=None):
        """One row per cell: original cell id first, then region, cell type, lineage, the full type path, confidence,
        expression-only label, marker cluster, whether it was hand-labeled, and coordinates."""
        r = self.result(version)
        types = np.array(self.types, dtype=object)
        lab = types[r["pred"]]
        paths = {t: self.type_path(t) for t in self.types}
        depth = max(len(p) for p in paths.values())
        df = pd.DataFrame({"cell_id": self.cell_ids(),
                           "region": np.where(r["region"] == "unknown", "", r["region"]).astype(object),
                           "cell_type": lab,
                           "lineage": np.array([paths[t][0] for t in self.types], dtype=object)[r["pred"]]})
        for k in range(depth):
            lut = np.array([paths[t][k] if k < len(paths[t]) else "" for t in self.types], dtype=object)
            df[f"level_{k + 1}"] = lut[r["pred"]]
        df["confidence"] = r["conf"].round(4)
        df["expression_label"] = types[r["base_pred"]]
        if self.scheme["kind"] == "markers":
            df["marker_cluster"] = self.marker_result(r["state"])["cluster_of"]
        df["hand_labeled"] = np.isin(np.arange(len(self.cells)), r["locked_rows"])
        df["x_um"], df["y_um"] = self.cells.x.round(2).to_numpy(), self.cells.y.round(2).to_numpy()
        return df

    def export(self, version=None):
        r = self.result(version)
        v = r["state"]["version"]
        d = self.dir / "exports" / f"v{v:04d}"
        d.mkdir(parents=True, exist_ok=True)
        self.export_table(version).to_csv(d / "cells.csv", index=False)
        json.dump(r["state"], open(d / "state.json", "w"), indent=1)
        r["conflicts"].to_csv(d / "conflicts.csv", index=False)
        lines = [f"# Session {self.id}, sample {self.sample_id}, version {v}", ""]
        for h in self.history():
            if h["version"] <= v:
                lines.append(f"- v{h['version']} ({h['created_at']}, {h['action']}): {h['message']}")
        lines += ["", "## Constraints"] + [f"- {cmod.describe(c)}  (\"{c['text']}\")" for c in r["state"]["constraints"]]
        (d / "CHANGELOG.md").write_text("\n".join(lines) + "\n")
        return str(d)
