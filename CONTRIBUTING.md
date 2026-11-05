# Contributing to STAnnotator

Bug reports, ideas and pull requests are welcome.

## Report a bug or suggest a feature

Open an [issue](https://github.com/caokai1073/STAnnotator/issues) with:
- the version (About in the app) and your system (macOS, Windows or Linux);
- the platform of the data (Xenium, CosMx, MERSCOPE, Visium HD, Stereo-seq, h5ad / RDS) and its size;
- what you did, what you expected and what happened;
- the end of the log if the app stopped: `~/stagent/app.log` (Windows: `%LOCALAPPDATA%\STAnnotator\app.log`).

## Set up the source

You need [uv](https://docs.astral.sh/uv/).

```bash
git clone https://github.com/caokai1073/STAnnotator.git && cd STAnnotator
uv venv .venv --python 3.12
uv pip install --python .venv/bin/python -e ".[app]"
.venv/bin/stagent app            # the window; --browser for a browser tab
```

Set `STAGENT_HOME` to a scratch folder while developing, so that your tests do not touch your own samples.
The code layout is described at the end of the [user guide](docs/GUIDE.md#layout).

## Tests

`tests/` holds smoke tests, each a script run directly, e.g.

```bash
STAGENT_HOME=/tmp/stagent-test .venv/bin/python tests/smoke_rds.py
```

`smoke_rds`, `smoke_joint`, `smoke_import_names`, `smoke_knn` and `smoke_file_filters` run on the data in this
repository (`tests/data`, `examples/`) or on synthetic cells; `smoke_engine`, `smoke_align` and `smoke_mcp` need
sections of the SCP2927 IBD Xenium atlas imported (see the user guide). Run the ones that touch your change.

## Pull requests

- One change per pull request, with a short description of what it changes for the user and how you tested it.
- Keep the style of the code around it; no new dependency without a reason.
- Results shown to the user must come from the data: no invented numbers, cell types or regions.

## License of contributions

STAnnotator is licensed under the [PolyForm Noncommercial License 1.0.0](LICENSE). By submitting a contribution
(a pull request, a patch or code in an issue) you confirm that you wrote it or have the right to submit it, and you
grant the STAnnotator authors a perpetual, worldwide, non-exclusive, royalty-free, irrevocable license to use, copy,
modify, sublicense and distribute your contribution under the PolyForm Noncommercial License and under any other
license, including commercial licenses. You keep the copyright of your contribution.
