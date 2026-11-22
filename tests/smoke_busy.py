"""Runs that must not overlap: an annotation waits for a run that moves the sample's cells (an alignment, a
registration check, a filter, an import: once, an annotation beside an automatic alignment froze the app, OpenBLAS
threads spinning in np.matmul while the neighbor search changed BLAS's thread count); the neighbor search's change of
BLAS's thread count and the other heavy products share one gate (markers.blas_gate); /api/annotate refuses fields that
are not text; an annotation an app left half made goes to the trash at the next start. Runs in a temporary data
folder."""
import json, os, shutil, tempfile, threading, time
EX = os.path.join(os.path.dirname(__file__), "..", "examples", "cadinu_merfish")
home = tempfile.mkdtemp(prefix="stagent_busy_")
os.environ["STAGENT_HOME"] = home; os.environ["STAGENT_SETTINGS_DIR"] = home + "/cfg"
from stagent import importer, config, markers, webapp
from starlette.testclient import TestClient
try:
    h5 = os.path.join(home, "a.h5ad")
    shutil.copy(os.path.join(EX, "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), h5)
    importer.import_sample("s1", h5, progress=lambda *a: None)
    c = TestClient(webapp.build_app(8799), base_url="http://127.0.0.1:8799")
    marker_text = open(os.path.join(EX, "mouse_colon_merfish_markers.txt")).read()
    # not text: 400, not 500
    r = c.post("/api/annotate", json={"sample_id": "s1", "markers_text": {"a": 1}})
    assert r.status_code == 400 and "text" in r.json()["error"], (r.status_code, r.text)
    # an alignment running on the sample: Annotate waits (409), the placeholder the page makes on opening does not
    webapp.JOBS["autoalign:s1"] = {"stage": "searching", "progress": 0.1, "done": False, "error": None, "kind": "autoalign"}
    r = c.post("/api/annotate", json={"sample_id": "s1", "markers_text": marker_text})
    assert r.status_code == 409 and "busy" in r.json()["error"], (r.status_code, r.text)
    r = c.post("/api/annotate", json={"sample_id": "s1", "method": "none", "mode": "new"})
    assert r.status_code == 200, r.text
    # an embedding of the sample does not hold an annotation back
    webapp.JOBS["autoalign:s1"]["done"] = True
    webapp.JOBS["embed:s1"] = {"stage": "UMAP", "progress": 0.5, "done": False, "error": None, "kind": "embed"}
    assert not webapp._placing("s1")
    webapp.JOBS["embed:s1"]["done"] = True
    # the gate: while another thread holds it, the neighbor search's change of BLAS's thread count waits
    held, entered = threading.Event(), []
    def hold():
        with markers.blas_gate():
            held.set(); time.sleep(1.0)
    t = threading.Thread(target=hold); t.start(); held.wait(5)
    t0 = time.time()
    with markers._one_blas_thread():
        entered.append(time.time() - t0)
    t.join()
    assert entered[0] >= 0.8, f"the search changed BLAS's threads while another thread held the gate ({entered[0]:.2f} s)"
    with markers.blas_gate():  # (the same thread: re-entered, no deadlock)
        with markers._one_blas_thread():
            pass
    # an annotation an app left half made (its app ended meanwhile): to the trash at the next start; one an app that
    # still runs is computing stays
    for name, pid in (("s1_left", 2 ** 22 + 12345), ("s1_running", os.getppid())):
        d = config.SESSIONS / name
        d.mkdir(parents=True)
        json.dump({"session_id": name, "sample_id": "s1", "method": "markers"}, open(d / "session.json", "w"))
        (d / webapp.UNFINISHED).write_text(str(pid))
    gone = webapp._drop_unfinished()
    assert gone == ["s1_left"] and not (config.SESSIONS / "s1_left").exists() and (config.SESSIONS / "s1_running").exists(), gone
    shutil.rmtree(config.SESSIONS / "s1_running")  # (a bare folder of this test, not a whole annotation)
    # an annotation that completes leaves no mark
    r = c.post("/api/annotate", json={"sample_id": "s1", "markers_text": marker_text}).json()
    assert "session_id" in r, r
    sid = r["session_id"]
    while not (j := c.get(f"/api/job/{sid}").json())["done"]:
        time.sleep(1)
    assert not j.get("error") and not (config.SESSIONS / sid / webapp.UNFINISHED).exists(), j
    print("smoke_busy ok")
finally:
    shutil.rmtree(home, ignore_errors=True)
