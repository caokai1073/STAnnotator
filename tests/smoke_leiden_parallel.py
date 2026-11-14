"""Leiden at several resolutions in worker processes (markers._leiden_many) gives the clusters of running them one by
one, to the bit: the marker annotation of the bundled mouse colon example with the workers forced on (threshold
lowered, the spawn check overridden) against the same with them off; and the switch: automatic by default (workers by the free
memory left after the annotation's own needs), off with the setting leiden_parallel (the app's "Parallel
clustering") or STAGENT_LEIDEN_PARALLEL=0.
Run: python tests/smoke_leiden_parallel.py  (a temporary data folder removed at the end)."""
import os, shutil, tempfile
from pathlib import Path


def main():
    home = Path(tempfile.mkdtemp(prefix="stagent_leiden_"))
    os.environ["STAGENT_HOME"] = str(home)
    os.environ["XDG_CONFIG_HOME"] = str(home / "cfg")  # (the settings of this run: not the user's)
    os.environ.pop("STAGENT_LEIDEN_PARALLEL", None)
    import numpy as np
    from stagent import importer, markers, config
    from stagent.session import Session
    # automatic by default where workers are safe: as many as the free memory allows after the annotation's own
    # needs; the setting (the app's switch) turns it off; the variable overrides it
    safe = markers._spawn_safe
    markers._spawn_safe = lambda: True  # (as in the app)
    on = markers.leiden_workers(700_000, 9e6, 4, 0.0, 64.0)
    assert on >= 2 or (os.cpu_count() or 1) < 3, "not on by default"
    assert markers.leiden_workers(700_000, 9e6, 4, 0.0, 3.0) == 0, "workers beyond the free memory"
    assert markers.leiden_workers(700_000, 9e6, 4, 60.0, 64.0) == 0, "workers in the memory the annotation needs"
    assert markers.leiden_workers(700_000, 9e6, 4, 0.0, 6.0) == 2, "not as many as fit"  # (0.8 x 6 GB / 1.8 GB)
    config.write_settings(leiden_parallel=False)
    assert markers.leiden_workers(700_000, 9e6, 4, 0.0, 64.0) == 0, "the switch off ignored"
    config.write_settings(leiden_parallel=True)
    os.environ["STAGENT_LEIDEN_PARALLEL"] = "0"
    assert markers.leiden_workers(700_000, 9e6, 4, 0.0, 64.0) == 0, "STAGENT_LEIDEN_PARALLEL=0 ignored"
    markers._spawn_safe = safe
    os.environ["STAGENT_LEIDEN_PARALLEL"] = "1"  # (this script is guarded: workers may start)
    ex = Path(__file__).resolve().parents[1] / "examples" / "cadinu_merfish"
    try:
        importer.import_sample("M", str(ex / "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), progress=lambda *a: None)
        text = (ex / "mouse_colon_merfish_markers.txt").read_text()
        markers.LEIDEN_PARALLEL_MIN = 1000
        assert markers.leiden_workers(32002, 400000, 4) >= 2, "no workers on this computer (memory or cores)"
        out, used = {}, []
        orig = markers.leiden_workers
        markers.leiden_workers = lambda *a: used.append(orig(*a)) or used[-1]  # (the workers each level got)
        for par in (False, True):
            markers.LEIDEN_PARALLEL = par
            s = Session.create("M", method="markers", markers_text=text, knn="exact")
            res = s.marker_result()
            out[par] = (np.asarray(res["labels"]).astype(str), np.asarray(res["cluster_of"]).astype(str),
                        res["nodes"], np.asarray(s.labels()).astype(str))
        assert max(used) >= 2, "the workers were never used"
        a, b = out[False], out[True]
        assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1]), "marker step differs with workers"
        assert a[2] == b[2], "resolutions chosen differ with workers"
        assert np.array_equal(a[3], b[3]), "final labels differ with workers"
        print(f"{len(a[3]):,} cells: the same clusters and labels with Leiden in worker processes")
        print(f"the switch: automatic by default ({on} workers for 700,000 cells with 64 GB free), off, overridden")
        print("ok")
    finally:
        shutil.rmtree(home, ignore_errors=True)


if __name__ == "__main__":
    main()
