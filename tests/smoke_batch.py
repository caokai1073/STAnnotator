"""stagent batch: a manifest checked before anything runs, two samples annotated in worker processes, the CSV files
and the summary written, a second run keeps the steps done, a sample that fails does not stop the others, shards
split the rows. Runs in a temporary folder (no H&E: the registration is tested on real data, not here)."""
import csv, json, os, shutil, subprocess, sys, tempfile
from pathlib import Path
EX = Path(__file__).resolve().parents[1] / "examples" / "cadinu_merfish" / "cadinu_merfish_D0_m3a_1_slice_2.h5ad"
tmp = Path(tempfile.mkdtemp(prefix="stagent_batch_"))
os.environ["STAGENT_SETTINGS_DIR"] = str(tmp / "cfg")  # (never the user's settings)


def write(rows, name="m.csv"):
    p = tmp / name
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return p


def batch(*args, env=None):
    r = subprocess.run([sys.executable, "-m", "stagent.cli", "batch", *map(str, args)], capture_output=True, text=True,
                       env=dict(os.environ, STAGENT_HOME=str(tmp / "unused"), **(env or {})))
    return r.returncode, r.stdout + r.stderr


try:
    # the manifest's problems, all at once, before anything runs
    bad = write([{"sample": "a b", "h5ad": str(EX), "organ": ""}, {"sample": "x", "h5ad": "/no/such.h5ad", "organ": "colon"},
                 {"sample": "x", "h5ad": str(EX), "organ": "colon"}], "bad.csv")
    rc, out = batch(bad, "--out", tmp / "o0")
    assert rc != 0 and "a b" in out and "not found" in out and "twice" in out and "no marker list" in out, out
    assert not (tmp / "o0").exists(), "ran with a bad manifest"
    # two samples (one by its own marker list file, one by --organ), one that fails (a wrong list)
    lst = tmp / "list.txt"
    lst.write_text((EX.parent / "mouse_colon_merfish_markers.txt").read_text())
    (tmp / "broken.txt").write_text("this is not a marker list\n")
    m = write([{"sample": "s1", "h5ad": str(EX), "markers": str(lst)}, {"sample": "s2", "h5ad": str(EX), "markers": ""},
               {"sample": "s3", "h5ad": str(EX), "markers": str(tmp / "broken.txt")}])
    rc, out = batch(m, "--out", tmp / "o", "--jobs", "2", "--organ", "colon")
    o = tmp / "o"
    assert rc == 1, out  # (one failed)
    summ = {r["sample"]: r for r in csv.DictReader(open(o / "summary.csv"))}
    assert summ["s1"]["status"] == "ok" and summ["s2"]["status"] == "ok" and summ["s3"]["status"] == "failed", summ
    assert summ["s3"]["error"], summ["s3"]
    for s in ("s1", "s2"):
        a = list(csv.DictReader(open(o / "results" / s / "annotation.csv")))
        assert len(a) == int(summ[s]["n_cells"]) > 30000 and {"cell_id", "cell_type", "lineage", "x_um", "y_um"} <= set(a[0]), a[0]
        p = list(csv.DictReader(open(o / "results" / s / "cell_positions.csv")))
        assert len(p) == len(a) and set(p[0]) == {"cell_id", "x_um", "y_um", "px", "py"}
        assert not (o / "results" / s / "registration.csv").exists()  # (no H&E)
        q = json.load(open(o / "results" / s / "qc.json"))
        assert q["n_cells"] == len(a) and summ[s]["registration"] == "none (no H&E)"
    assert summ["s1"]["session"] != summ["s2"]["session"]
    # again: the steps done are kept (nothing annotated again); the failed one is tried again
    rc, out = batch(m, "--out", o, "--organ", "colon")
    log = (o / "logs" / "s1.log").read_text()
    assert log.count("annotate: done before, kept") == 1 and log.count("annotate ...") == 1, log[-800:]
    assert {r["sample"]: r for r in csv.DictReader(open(o / "summary.csv"))}["s1"]["session"] == summ["s1"]["session"]
    # the marker list changed: only annotate and export run again (not import, not register)
    lst.write_text(lst.read_text() + "\n")
    lst.write_text(lst.read_text().replace("\n\n", "\n# edited\n", 1))
    rc, out = batch(m, "--out", o, "--organ", "colon", "--only", "s1")
    log = (o / "logs" / "s1.log").read_text().split("==== ")[-2]
    assert "import: done before, kept" in log and "register: done before, kept" in log and "annotate ..." in log, log
    # a process that dies (as out of memory) fails its own sample only, and keeps its finished steps
    rc, out = batch(m, "--out", tmp / "o4", "--jobs", "2", "--organ", "colon", "--only", "s1", "s2",
                    env={"STAGENT_BATCH_TEST_CRASH": "s1"})
    summ = {r["sample"]: r for r in csv.DictReader(open(tmp / "o4" / "summary.csv"))}
    assert summ["s1"]["status"] == "failed" and "died" in summ["s1"]["error"] and summ["s2"]["status"] == "ok", summ
    assert "import" in json.load(open(tmp / "o4" / "status" / "s1.json"))["steps"], "finished steps lost"
    # --drop-data removes only the sample's own data (s2 and s2_b: one name the start of the other)
    m2 = write([{"sample": "s2", "h5ad": str(EX), "markers": ""}, {"sample": "s2_b", "h5ad": str(EX), "markers": ""}], "m2.csv")
    rc, out = batch(m2, "--out", tmp / "o5", "--jobs", "2", "--organ", "colon", "--drop-data")
    summ = {r["sample"]: r for r in csv.DictReader(open(tmp / "o5" / "summary.csv"))}
    assert summ["s2"]["status"] == summ["s2_b"]["status"] == "ok", summ
    assert not list((tmp / "o5" / "data" / "sessions").glob("*")) and not list((tmp / "o5" / "data" / "samples").glob("*"))
    rc, out = batch(m2, "--out", tmp / "o5", "--organ", "colon")  # (results there: nothing to do again)
    assert "annotate ..." not in (tmp / "o5" / "logs" / "s2.log").read_text().split("==== ")[-2]
    # names not in the manifest, a malformed shard: refused with a message
    rc, out = batch(m, "--out", tmp / "o6", "--organ", "colon", "--only", "nope")
    assert rc != 0 and "not in the manifest" in out, out
    rc, out = batch(m, "--out", tmp / "o6", "--organ", "colon", "--shard", "x")
    assert rc != 0 and "I/N" in out and "Traceback" not in out, out
    # shards: rows 0, 2 / row 1
    rc, out = batch(m, "--out", tmp / "o2", "--shard", "1/2", "--organ", "colon", "--only", "s2")
    assert [r["sample"] for r in csv.DictReader(open(tmp / "o2" / "summary.csv"))] == ["s2"], out
    rc, out = batch(m, "--out", tmp / "o3", "--shard", "0/2", "--organ", "colon", "--only", "s2")
    assert "0 samples" in out, out
    print("ok")
finally:
    shutil.rmtree(tmp, ignore_errors=True)
