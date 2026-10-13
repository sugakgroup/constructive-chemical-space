# CCS S5 — Worked examples

The 55 numbered examples in Section S5 of *Defining Chemical Space Constructively*
have **one executable Jupyter notebook per SI record**. Together they cover all
**84 fully specified electron-allocation forms**.

Browse the [complete numbered index](NOTEBOOK_INDEX.md) (all 55 molecule names,
SI sections, Form labels and direct notebook links), start with
[S5_01_H2.ipynb](notebooks/S5_01_H2.ipynb), or open the
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
- `execution_report.json` stores the **latest complete** execution run ID,
  per-notebook kernel execution counts, SVG counts and SHA-256 checksums.
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

The currently committed 55 Notebooks were regenerated using **real Jupyter
ipykernels** on October 13, 2026, following the generic ring-aware layout
improvements. The [complete Actions run 38388645967](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38388645967)
successfully executed all 55 notebooks, 84 electron-allocation forms and
**168 actual SVG outputs** (84 Constitution diagrams and 84 allocation
diagrams). Each notebook contains six genuinely executed code cells. The
environment was Python 3.11.17, RDKit 2026.09.1 and nbclient 0.11.0.
The source/output SHA-256 fingerprints and per-notebook execution counts are
recorded in `execution_report.json`. [Commit a9d27c7](https://github.com/sugakgroup/constructive-chemical-space/commit/a9d27c72a55837aa84c97500d6b6b581a25d8df1)
contains the real-kernel outputs. The test suite passed **492 tests before
and after** this commit in the complete regeneration workflow.

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
    constitution.show(order=state.address_table, show_indices=False)
    state.form.show(order=state.address_table, show_indices=False,
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

## Compact atom-local electron and charge annotations

Allocation pictures now place the CCS-assignment imbalance as `+`, `−`,
`2−`, etc., adjacent to the nucleus. For groups 13–17, the non-core residual
valence count is displayed as `•` glyphs (amine nitrogen has `••`, water
oxygen has four dots). These are explicit electron-allocation bookkeeping
marks, **not inferred spin multiplicities or Lewis formal charges**.
They are drawn in SVG/PNG and have indexed machine-readable SVG metadata.
The default Notebook view suppresses `#0`/`#1` index labels.

All 55 direct-connection graphs are checked for graph planarity, but
**a chemically intelligible ring/cage projection takes priority over
eliminating line crossings**. Even a planar graph may be deliberately drawn
with crossings. Four S5 graphs are mathematically nonplanar (records 21, 42,
43 and 44) and receive a warning that a crossing-free planar embedding is
impossible. Every depiction remains a non-normative 2D layout, not a
calculated 3D structure.

## Record 13 (bis-periazulene) author correction

The separately validated B/C electron-allocation data use
`(C(9),C(14)):(2,2)` and `(C(13),C(14)):(1,1)`. This corrects the earlier
transcription, which had these two values reversed in Forms B/C. Other
connections are inherited from Form A unless specified; Form C retains
`n(C(0))=7`, `n(C(7))=5`. See `S5_VALIDATION.md` for provenance.

### Generic ring-aware 2D depiction

The same graph-only candidate selection applies to every S5 record: RDKit,
CoordGen, and an additional CoordGen candidate with temporarily weakened
high-degree hub connections where needed. The score prioritizes intact
short-ring and cage-face shapes, the position of any vertex connected to a
complete short ring, reasonable connection lengths, and separation between
nonbonded nuclei. Connection-line intersections have only a weak penalty.

The coordinate-generation models are reconstructed from **sorted undirected
nuclear connections**, independently of the electron-allocation drawing
styles. This allows each Form's Constitution and ElectronAllocation to share
the same skeletal coordinates without changing their original `Nuc`,
`Ele`, `Con`, endpoint allocations or actual line/arrow glyphs.

Regression tests confirm that the generic criteria select RDKit for S5-18
Ferrocene and CoordGen for S5-47/S5-48 Perfluorocubane, with no species-
specific or element-specific branches. Actual output depends on the toolkit
version and, for CoordGen, its possible local-minimum variability.

The residual-electron dots are genuine SVG `<circle>` primitives (PNG:
ellipses). Endpoint counts, atomic imbalance signs, and optional indexes
retain machine-readable metadata in the output. Use `layout="auto"` by
default or `layout="rdkit"`/`layout="coordgen"` for explicit comparisons.

### Validation and CI costs

Ordinary source/test pushes run **representative Python 3.11 smoke tests**
only. [Run 38387790562](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38387790562)
passed 53 selected tests, executed 9 genuine Jupyter notebooks for 12 Forms
and produced 24 actual SVG figures. Its preview artifact is temporary.

The explicitly triggered [full-corpus workflow 38388645967](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38388645967)
then regenerated **all 55 notebooks and 168 SVGs**, ran **492 tests twice**,
and committed their authentic outputs to PR #4.
The Python 3.10–3.12 complete matrix remains deferred until explicitly
requested. PR #4 remains Draft and is not merged.

For direct access to each executed molecular notebook, see the
[complete S5 notebook index](NOTEBOOK_INDEX.md).
