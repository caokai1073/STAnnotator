"""A project of two slices (the example MERFISH slice cut in halves): import, joint annotation, joint embedding with
Harmony, joint niches, composition. Runs in a temporary data folder."""
import os, shutil, tempfile, time
import h5py
import numpy as np
EX = os.path.join(os.path.dirname(__file__), "..", "examples", "cadinu_merfish")
home = tempfile.mkdtemp(prefix="celltypeagent_joint_")
os.environ["CELLTYPEAGENT_HOME"] = home
from celltypeagent import importer, joint, projects
try:
    h5 = os.path.join(home, "two.h5ad")
    shutil.copy(os.path.join(EX, "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), h5)
    with h5py.File(h5, "r+") as h:  # a categorical obs column "half": left / right of the middle
        x = h["obs"]["x"][:]
        g = h["obs"].create_group("half")
        g.attrs.update({"encoding-type": "categorical", "encoding-version": "0.2.0", "ordered": False})
        g["codes"] = (x >= (x.min() + x.max()) / 2).astype(np.int8)
        g.create_dataset("categories", data=np.array(["left", "right"], dtype=object), dtype=h5py.string_dtype())
    t = time.time()
    importer.import_project("demo", "h5ad", [h5], sections=["left", "right"], section_column="half", progress=lambda *a: None)
    p = projects.list_projects()[0]
    print("project:", p["project_id"], p["samples"], [p["slices"][s]["metadata"] for s in p["samples"]], f"{time.time() - t:.0f}s")
    t = time.time()
    joint.annotate("demo", open(os.path.join(EX, "mouse_colon_merfish_markers.txt")).read())
    print("joint annotation:", {s: sid for s, sid in joint.latest_sessions("demo").items()}, f"{time.time() - t:.0f}s")
    t = time.time()
    e = joint.embed("demo", harmony=True)
    print("joint embedding:", {k: e[k] for k in e if k not in ("samples",)}, f"{time.time() - t:.0f}s")
    t = time.time()
    joint.niches("demo", n_niches=5)
    print("joint niches:", f"{time.time() - t:.0f}s")
    for mode in ("lineage", "niche", "leiden"):
        c = joint.composition("demo", mode)
        print(mode, {r["section"]: (r["n_cells"], {k: r["shares"].get(k, 0) for k in c["keys"][:3]}) for r in c["slices"]})
finally:
    shutil.rmtree(home)
