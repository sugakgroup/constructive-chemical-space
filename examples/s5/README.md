# CCS S5 worked examples

**55 executable notebooks / 84 electron-allocation forms / 168 saved SVGs.**

[Browse by SI record and molecule name](NOTEBOOK_INDEX.md) or open the [notebooks directory](notebooks/). Each notebook defines a literal `Word([Nuc(...), Ele(), Con(...)])` for each Form, executes `replay_state()`, derives its `ElectronAllocation` and `Constitution`, validates the replay against the independently stored SI record, and draws the *replay-generated* objects. All nuclei, including H, are explicit.

## Install and browse

From the repository root:

~~~bash
python -m pip install -e ".[dev]"
python -m pip install jupyterlab
jupyter lab
~~~

Open any file in `examples/s5/notebooks/`. The notebooks use the same implementation in `src/constructive_chemical_space/ccs.py`; there is no duplicate CCS core.

## Source data and saved outputs

| Resource | Purpose |
| --- | --- |
| [`data/manifest.json`](data/manifest.json) | Record index, SI location, Form labels and expected counts |
| [`data/s5_records.json.xz`](data/s5_records.json.xz) | Single complete source-data fixture for all 84 Forms |
| [`_common.py`](_common.py) | Shared example loading, indexed replay checks and optional export |
| [`execution_report.json`](execution_report.json) | Real-kernel run, per-notebook execution counts and SVG/source hashes |
| [`NOTEBOOK_INDEX.md`](NOTEBOOK_INDEX.md) | Direct links to the 55 saved notebooks |

The committed notebooks contain **genuine Jupyter kernel outputs**: 84 Constitution SVGs plus 84 ElectronAllocation SVGs. All 55 notebooks were last regenerated in [run 38391623509](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38391623509), which passed 494 validation tests before and after output commit [`5ef1bfb`](https://github.com/sugakgroup/constructive-chemical-space/commit/5ef1bfb6278331c248d109155472e4a9109f1968).

## Regenerate locally

The following commands *modify* the local Notebook files. `generate_explicit_words.py` first writes source-only notebooks with null execution counts and no outputs; `execute_notebooks.py` then runs real `ipykernel` sessions and saves only actual display data and provenance.

~~~bash
python -m pip install -e ".[dev]"
python -m pip install nbclient nbformat ipykernel
python examples/s5/generate_explicit_words.py
python examples/s5/execute_notebooks.py
python -m pytest -q tests/test_s5_examples.py tests/test_visualization.py
~~~

Alternatively, use the read-only GitHub Actions workflow [`execute_s5_notebooks.yml`](../../.github/workflows/execute_s5_notebooks.yml) to produce a downloadable archive without updating the canonical branch. Generated standalone graphics can be saved through `export_example(...)`; the optional `outputs/` directory is created on demand and ignored by Git.

## Scope and known source correction

The drawings retain direct CCS connections and distinguish all specified electron-allocation Forms. They are **not calculated 3D geometries**, conventional bond-order assignments, Lewis-charge determinations or stereochemical/topological identities. Negative assignment-imbalance signs are vector strokes; residual electrons are SVG circles.

**Record 13:** The validated Forms B/C fixture uses `(C(9), C(14)): (2,2)` and `(C(13), C(14)): (1,1)`; Form C retains `n(C(0)) = 7`, `n(C(7)) = 5`. This correction must also appear in the final manuscript SI. The full source-level audit, original-document hashes and scope of verification are documented in [S5_VALIDATION.md](../../S5_VALIDATION.md).
