# CCS S5 — Worked examples

The 55 numbered examples in Section S5 of *Defining Chemical Space Constructively*
have **one executable Jupyter notebook per SI record**. Together they cover all
**84 fully specified electron-allocation forms**.

Start with [S5_01_H2.ipynb](notebooks/S5_01_H2.ipynb), then browse the
[notebooks/](notebooks/) directory. Each notebook contains **literal**
`word_A = Word([Nuc(...), Ele(), Con(...)])` definitions (one per Form),
then actually calls `replay_state()`, obtains new `allocation_A` and
`constitution_A` objects, and displays **those replay-produced objects**
through `.show()`. The SI source Form is used independently to check that
the explicit Word and resulting objects are correct; it is not what is drawn.

## Installation

From the repository root:

```bash
python -m pip install -e ".[dev]"
python -m pip install jupyterlab
jupyter lab
```

Then open a file in `examples/s5/notebooks/`. The notebook bootstrap works
whether the Jupyter process starts at the repository root or inside that
notebook directory.

## Data and source of truth

- `data/manifest.json` indexes 55 notebooks, SI section/record numbers,
  source page locators, Table S4 expectations and each record's Form labels.
- `data/s5_records.json.xz` stores the complete immutable input dataset for
  84 forms (compressed UTF-8 JSON). This file is the **only record-data fixture**,
  shared by the notebooks and the verification tests.
- `_common.py` contains only convenience loading, verification and export
  utilities. It **does not duplicate or modify** the CCS implementation.
- `outputs/` is for optional locally generated SVG/PNG files and is ignored
  by Git. By contrast, genuine Jupyter `display_data` SVG output **is stored
  directly inside the 55 executed notebooks**; these were produced by
  `nbclient` using real `ipykernel` executions, not synthesized afterward.
- `execution_report.json` stores the original execution run ID, per-notebook
  kernel execution counts, SVG counts and SHA-256 checksums.
- `tests/test_s5_examples.py` verifies notebook completeness and all Forms.

No extra primitive is introduced. All nuclei including H remain explicit.
`W(D)` uses the S5 displayed index order and Con lexicographic ordering.
Index choices do not change intrinsic CCS identity. Named descriptions and
depiction strokes do not assert computed orbital properties or physical
realizability. Catenation, stereochemistry and spatial arrangement remain
outside constitutional identity under SI S1.

## Literal construction words and real execution provenance

The Notebook word-definition source is produced by
`generate_explicit_words.py` from the independently audited S5 dataset:
each Nuc, Ele and Con token is **written literally into the Notebook's code
cell**, rather than being computed at Notebook execution time by calling
`form.to_word()`. The following cell replays **every** Form's literal Word into
a fresh ElectronAllocation and Constitution; it checks the complete nuclear
identity, direct-connection and total-electron signatures under the actual
Replay address-table bijection. This proves that all Form-specific
constitutions are the same, without an exponential general isomorphism search.
Both **each Form's constructed Constitution** and its constructed
ElectronAllocation are drawn using that Form's own replay address order.
Asymmetric electron counts are annotated beside the nucleus contributing each
endpoint allocation (including one-electron or zero-allocation connections).

The output checked into PR #4 was produced by a **real Jupyter ipykernel**
for all 55 notebooks, yielding **168 real SVG images** (84 Form-specific
constitutions and 84 allocation diagrams). Each notebook has six actually
executed code cells. The environment was Python 3.11.17, RDKit 2026.09.1
and nbclient 0.11.0. Source and saved-output SHA-256 fingerprints appear in
`execution_report.json`. The original real-kernel execution log is
[GitHub Actions run 38357651689](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38357651689).

**No output or execution_count field is simulated.** To regenerate, execute
the following two commands together after installing the dependencies:

```bash
python -m pip install -e ".[dev]"
python -m pip install nbclient nbformat ipykernel
python examples/s5/generate_explicit_words.py  # rewrites source-only notebooks
python examples/s5/execute_notebooks.py         # real kernel execution and SVG capture
```

The first script deliberately resets every code cell to *unexecuted* (empty
outputs and null execution counts); the second uses NotebookClient to save
only genuine outputs after all assertions succeed. Do not commit the
source-only intermediate notebooks as though they were executed.

## Working in Python

```python
from examples.s5._common import build_example, verify_example, export_example

example = build_example(10)          # Phenol: Forms A, B, C and D
print(verify_example(example))        # {'A': ..., 'B': ..., ...}
for label, form in example.forms.items():
    word = form.to_word(order=example.order)  # the notebooks instead spell out every token
    state = word.replay_state()
    constitution = state.form.to_constitution()
    print(label, "Replay -> D -> M")
    constitution.show(order=state.address_table, show_indices=True)
    state.form.show(order=state.address_table, show_indices=True,
                    allocation_label_mode="endpoint")

# Explicitly save images; no files are generated by default:
export_example(example, formats=("svg",))
```

The notebooks default to **full, faithful explicit-H** visualization. The
optional `show_hydrogens=False` view drops actual H nuclei and is marked
**lossy**. It is not a complete CCS constitution.

## Verification

```bash
python -m pytest tests/test_s5_examples.py
```

The SI transcription/source audit, its source checksums and the difference
between source-document verification and executable-fixture replay are
recorded in [S5_VALIDATION.md](../../S5_VALIDATION.md). The original PDF/DOCX
are not automatically re-extracted during CI.

In the visualization CI environment, all 55 constitutions and all 84 forms are rendered as well-formed, full-explicit-H SVGs. This is an executable depiction check, not a claim that these formal electronic descriptions correspond to experimentally observed states.
