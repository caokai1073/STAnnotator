"""Background runs the user can stop (the app's Stop button under Annotate, Run PCA + UMAP, Niches, ...).

spawn() runs a job's work in a thread. stop() raises Stopped in that thread: at once in Python code, after the
current native call otherwise (one Leiden or UMAP pass), at the end of a critical() block (files being written), and
again every NUDGE_S until the run has ended (callbacks such as numba's swallow it); never while numba compiles.
What the run produced is then taken back, so the analysis is as it was before the click: a run registers how
(on_stop: e.g. delete the annotation it created), and replacing() keeps the files a run writes over until it ends
(put back when it is stopped, dropped when it completes). A new run of the same job waits for a stopped one to end.
"""
import ctypes, shutil, sys, threading, time
from contextlib import contextmanager
from pathlib import Path

WAIT_S = 30  # a new run waits this long for a stopped run of the same job to end
NUDGE_S = 1.0  # a stop is raised again this often until the run ends (code such as numba's callbacks swallows it)


class Stopped(BaseException):
    """Raised in a run's thread when the user stops it (a BaseException: `except Exception` lets it through)."""


class Run:
    def __init__(self, job):
        self.job, self.thread = job, None
        self.stop_asked = self.finished = False
        self.critical = 0
        self.undo, self.done = [], []  # on stop (in reverse order) / when it completes


RUNS = {}
_lock = threading.Lock()
_local = threading.local()


def current():
    return getattr(_local, "run", None)


def on_stop(fn):
    """Call fn if the current run is stopped (to take back what it made)."""
    r = current()
    if r is not None:
        r.undo.append(fn)


def checkpoint():
    r = current()
    if r is not None and r.stop_asked and not r.critical:
        raise Stopped()


@contextmanager
def critical():
    """A step a stop must not cut in two (writing a result's files): the stop takes effect after it."""
    r = current()
    if r is None:
        yield
        return
    with _lock:
        r.critical += 1
    try:
        yield
    finally:
        with _lock:
            r.critical -= 1
    checkpoint()


def _aside(paths, tag):
    moved = []
    for p in map(Path, paths):
        b = p.with_name(f".{p.name}.before-{tag}")
        if p.exists():
            p.replace(b)
            moved.append((p, b))
        else:
            moved.append((p, None))
    return moved


def _back(moved):
    for p, b in moved:
        if p.is_dir():
            shutil.rmtree(p, ignore_errors=True)
        else:
            p.unlink(missing_ok=True)
        if b is not None and b.exists():
            b.replace(p)


def _drop(moved):
    for _, b in moved:
        if b is not None:
            shutil.rmtree(b, ignore_errors=True) if b.is_dir() else b.unlink(missing_ok=True)


@contextmanager
def replacing(*paths):
    """Write a run's result over the files (or folders) of the last one: the old ones are kept aside until the run
    ends, so that stopping it puts them back. Not cut by a stop; an error inside puts them back at once."""
    r = current()
    with critical():
        moved = _aside(paths, id(r))
        try:
            yield
        except BaseException:
            _back(moved)
            raise
        if r is None:
            _drop(moved)
        else:
            r.undo.append(lambda: _back(moved))
            r.done.append(lambda: _drop(moved))


def _settle(r):
    """The end of a run: take back what it made if it was stopped, else drop what it kept aside."""
    while True:  # a stop that arrives just now is raised here, not lost
        try:
            with _lock:
                r.finished = True
                stopped = r.stop_asked
            break
        except Stopped:
            continue
    for fn in (reversed(r.undo) if stopped else r.done):
        try:
            fn()
        except BaseException as e:
            print(f"job {r.job}: could not {'take back' if stopped else 'clean up'}: {e}")


def spawn(job, work):
    """Run work() in a thread as the stoppable run `job` (work reports its progress and errors itself)."""
    prev = RUNS.get(job)
    if prev is not None and prev.thread is not None and prev.thread.is_alive():
        if not prev.stop_asked:
            raise RuntimeError("it is still running")
        prev.thread.join(WAIT_S)
        if prev.thread.is_alive():
            raise RuntimeError("the stopped run is still ending: try again in a moment")
    r = Run(job)
    RUNS[job] = r

    def body():
        _local.run = r
        try:
            try:
                work()
            except Stopped:
                pass
        except Stopped:  # raised again just as the run ended
            pass
        finally:
            _settle(r)
            _local.run = None
    r.thread = threading.Thread(target=body, daemon=True)
    r.thread.start()
    return r


def stop(job):
    """Stop run `job`: True when it was running (what it made is taken back), False when it had already ended."""
    r = RUNS.get(job)
    if r is None:
        return False
    with _lock:
        if r.finished or r.stop_asked:
            return not r.finished
        r.stop_asked = True
    threading.Thread(target=_nudge, args=(r,), daemon=True).start()
    return True


UNSAFE = ("/numba/", "/llvmlite/")  # code a stop must not cut into (compiling: native state could be left broken)


def _unsafe(ident):
    """True while the thread runs code of UNSAFE (e.g. numba compiling pynndescent on its first use)."""
    f = sys._current_frames().get(ident)
    while f is not None:
        if any(u in f.f_code.co_filename.replace("\\", "/") for u in UNSAFE):
            return True
        f = f.f_back
    return False


def _nudge(r):
    """Raise Stopped in the run's thread, again every NUDGE_S until it ends (not inside a critical block, whose end
    raises it, nor while it compiles)."""
    while r.thread is not None and r.thread.is_alive():
        with _lock:
            if r.finished:
                return
            if not r.critical and not _unsafe(r.thread.ident):
                ctypes.pythonapi.PyThreadState_SetAsyncExc(ctypes.c_ulong(r.thread.ident), ctypes.py_object(Stopped))
        time.sleep(NUDGE_S)
