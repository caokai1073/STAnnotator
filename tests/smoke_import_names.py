"""Import names and failed imports: slice names are never shortened (two long sections with a common prefix get
different samples), a project import that fails on its second slice leaves no sample behind and can be run again,
an import written over an existing sample replaces it as a whole, and no importing/ folder stays. Runs in a
temporary data folder."""
import os, shutil, tempfile
import h5py
import numpy as np
EX = os.path.join(os.path.dirname(__file__), "..", "examples", "cadinu_merfish")
home = tempfile.mkdtemp(prefix="stagent_names_")
os.environ["STAGENT_HOME"] = home
from stagent import config, importer, projects
quiet = lambda *a: None
try:
    h5 = os.path.join(home, "two.h5ad")
    shutil.copy(os.path.join(EX, "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), h5)
    long = "a_section_name_that_is_long_enough_to_pass_the_old_eighty_character_limit_" + "x" * 20
    secs = [long + "_left", long + "_right"]  # the same 95 first characters
    with h5py.File(h5, "r+") as h:
        x = h["obs"]["x"][:]
        g = h["obs"].create_group("half")
        g.attrs.update({"encoding-type": "categorical", "encoding-version": "0.2.0", "ordered": False})
        g["codes"] = (x >= (x.min() + x.max()) / 2).astype(np.int8)
        g.create_dataset("categories", data=np.array(secs, dtype=object), dtype=h5py.string_dtype())
    names = [projects.slice_name("demo", s) for s in secs]
    assert len(set(names)) == 2 and all(n == f"demo_{s}" for n, s in zip(names, secs)), names
    # 1. a failed second slice: the first one is removed again
    try:
        importer.import_project("demo", "h5ad", [h5], sections=[secs[0], "no_such_section"], section_column="half", progress=quiet)
        raise AssertionError("expected the import to fail")
    except ValueError as e:
        print("failed as expected:", e)
    assert not (config.SAMPLES / names[0]).exists(), "the first slice stayed behind"
    assert not (config.HOME / "importing").exists() or not any((config.HOME / "importing").iterdir())
    # 2. run again, as it was: both slices, full names
    p = importer.import_project("demo", "h5ad", [h5], sections=secs, section_column="half", progress=quiet)
    assert p["samples"] == names, p["samples"]
    n = {s: importer._meta(s)["n_cells"] for s in names}
    print("imported", n)
    # 3. the same name again without overwrite: refused before anything is written
    try:
        importer.import_sample(names[0], h5, section_column="half", section=secs[0], progress=quiet)
        raise AssertionError("expected 'already exists'")
    except ValueError as e:
        assert "already exists" in str(e), e
    # 4. overwrite replaces the sample as a whole (a stale file of the earlier import goes)
    (config.SAMPLES / names[0] / "regcheck.json").write_text("{}")
    importer.import_sample(names[0], h5, section_column="half", section=secs[1], overwrite=True, progress=quiet)
    assert importer._meta(names[0])["n_cells"] == n[names[1]]
    assert not (config.SAMPLES / names[0] / "regcheck.json").exists()
    # 5. a name longer than the limit is refused with a message, not shortened
    try:
        importer._check_name("y" * 201, False)
        raise AssertionError("expected a length error")
    except ValueError as e:
        print("refused:", e)
    print("ok")
finally:
    shutil.rmtree(home, ignore_errors=True)
