"""The analyses run in the app (report.log_analysis): one line per run in every annotation, with its settings (no cell
lists, no marker list, no words to the agent), seconds and how it ended (done, error, stopped: webapp._logged around
a background run); the log keeps its newer half past its size; a bug report and a usage report carry the lines of
their time (analysis_runs.jsonl). Runs in a temporary data folder."""
import json, os, shutil, tempfile, threading, time, zipfile

home = tempfile.mkdtemp(prefix="stagent_anlog_")
os.environ["STAGENT_HOME"] = home; os.environ["STAGENT_SETTINGS_DIR"] = home + "/cfg"
from stagent import report, jobs, webapp

try:
    lines = lambda: [json.loads(x) for x in report.RUNS.read_text(encoding="utf-8").splitlines()]
    report.log_analysis("compare", "de", {"session": "s1", "version": 3, "by": "cell type", "group": "B cell",
                                          "selection": {"s1": list(range(250)), "s2": list(range(50))},
                                          "markers_text": "T cell,,CD3E\nB cell,,MS4A1\n\n", "message": "my words",
                                          "genes": ["G%d" % i for i in range(400)]},
                        seconds=1.234, session="s1")
    r = lines()[-1]
    assert r["section"] == "compare" and r["kind"] == "de" and r["outcome"] == "done" and r["seconds"] == 1.2, r
    st = r["settings"]
    assert st["selection_cells"] == 300 and st["marker_lines"] == 2 and st["message"] is True, st
    assert "session" not in st and "version" not in st and "selection" not in st and "markers_text" not in st
    assert st["by"] == "cell type" and isinstance(st["genes"], str) and len(st["genes"]) < report.SETTING_MAX + 40
    assert "my words" not in report.RUNS.read_text(encoding="utf-8")

    # a background run noted when it ends: done, error, stopped
    def run(job, work):
        webapp.JOBS[job] = {"stage": "starting", "progress": 0.0, "done": False, "error": None}
        r = jobs.spawn(job, webapp._logged(job, work, {"section": "embed", "kind": "embed", "request": {"n_pcs": 30},
                                                       "session": "s1"}))
        return r

    def ok():
        webapp.JOBS["j1"].update(done=True)

    def bad():
        webapp.JOBS["j2"].update(error="ValueError: no cells", done=True)
    for job, work in (("j1", ok), ("j2", bad)):
        run(job, work).thread.join(10)
    go = threading.Event()

    def long():
        go.set()
        while True:
            time.sleep(0.05)
    r3 = run("j3", long)
    go.wait(5)
    jobs.stop("j3")
    r3.thread.join(10)
    got = [(x["outcome"], x.get("error")) for x in lines()[-3:]]
    assert got == [("done", None), ("error", "ValueError: no cells"), ("stopped", None)], got
    assert all(x["settings"] == {"n_pcs": 30} and x["seconds"] >= 0 for x in lines()[-3:])

    # the reports: the lines of their time, in every annotation
    out = os.path.join(home, "r.zip")
    res = report.build(out, {"version": "t"}, days=3, parts=("analyses",))
    assert res["included"]["analyses"] == 4, res
    with zipfile.ZipFile(out) as z:
        assert len(z.read("analysis_runs.jsonl").decode().splitlines()) == 4
        assert "analysis_runs.jsonl" in z.read("README.txt").decode()
    assert report.contents(3)["analyses"]["n"] == 4
    res = report.build(out, {"version": "t"}, days=3, parts=("calls",))
    with zipfile.ZipFile(out) as z:
        assert "analysis_runs.jsonl" not in z.namelist()
    # older than the report's time: left out
    old = dict(lines()[0], time="2000-01-01T00:00:00")
    with open(report.RUNS, "a", encoding="utf-8") as f:
        f.write(json.dumps(old) + "\n")
    assert len(report._recent_runs(time.time() - 86400)) == 4

    # the log keeps its newer half past its size
    keep = report.RUNS_MAX
    try:
        report.RUNS_MAX = 20000
        for i in range(200):
            report.log_analysis("compare", "composition", {"rows": "cell type", "i": i})
        assert report.RUNS.stat().st_size <= 20000 and lines()[-1]["settings"]["i"] == 199
    finally:
        report.RUNS_MAX = keep
    print("smoke_analysis_log ok")
finally:
    shutil.rmtree(home, ignore_errors=True)
