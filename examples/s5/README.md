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
[GitHub Actions run 38362379662](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38362379662).

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

All 55 S5 direct-connection graphs are tested for 2D planarity. 51 are planar,
so the renderer enforces **no straight-line connection crossings**.
Four are mathematically nonplanar (records 21, 42, 43 and 44); their diagrams
state that crossing-free 2D representation is impossible. The geometry is a
non-normative layout, not a computed 3D structure.

## Record 13 (bis-periazulene) author correction

The separately validated B/C electron-allocation data use
`(C(9),C(14)):(2,2)` and `(C(13),C(14)):(1,1)`. This corrects the earlier
transcription, which had these two values reversed in Forms B/C. Other
connections are inherited from Form A unless specified; Form C retains
`n(C(0))=7`, `n(C(7))=5`. See `S5_VALIDATION.md` for provenance.

### Latest real-kernel reproducibility check

At the latest audited execution (`38362379662`), all 55 notebooks ran in
actual Jupyter kernels and produced 168 SVGs: one Constitution and one
ElectronAllocation per Form. Every picture uses `show_indices=False`,
with explicit hydrogen nuclei, no generated 3D geometry, and the
crossing-aware 2D layout. The output includes nucleus-local CCS q signs
and residual valence dots, plus warnings where 2D planarity is impossible.
The execution report's captured SVG/source hashes are independently checked
by `tests/test_s5_examples.py`.

### Generic graph-based 2D and automatic drawing-space expansion

There are now **no molecule-, metal-, cluster-, ligand-, or named-SI-record
layout branches**. For every graph the same code compares candidate RDKit,
CoordGen and, when crowded, seeded spring-based arrangements. It selects
using straight-edge intersections, interatomic separation and extreme
connection-length penalties. Planar graphs always use an embedding with
**zero straight-line crossings**, while nonplanar graphs show a visible
warning. No physically computed 3D geometry is implied.

The renderer alone cannot create additional visible spacing simply by
multiplying coordinate positions: RDKit automatically fits them to the
drawing viewport. Instead every S5 notebook calculates the **maximum
connection degree** and enlarges both SVG width and height when crowded,
leaving enough pixels between symbols, connection lines and nuclear
electron annotations. For ferrocene the factor is approximately 1.77
because its metal center has ten incident connections; no special
ferrocene code or shape rules are used.

The residual-electron dots are genuine SVG `<circle>` primitives and
PNG ellipses rather than UTF-8 bullet characters, preventing escaped-byte
glyph problems. Their electron counts remain in machine-readable SVG
metadata. `layout="auto"` is recommended; `"rdkit"` and `"coordgen"`
remain available to compare generic candidate behavior.

### Cost-conscious testing and real output provenance

During PR review, ordinary pushes run **only one representative Python
3.11 workflow** and execute eight SI notebooks with genuine Jupyter
kernels. These example outputs and the real SVG previews are available
as a short-lived GitHub Actions ZIP artifact (5-day retention), not
committed as substitute full-corpus results.

The generic rendering and congestion-scale revision was checked by
[representative run 38373858558](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38373858558):
**47 selected tests passed**, and eight actual notebooks generated
**22 authentic SVGs**. The [corresponding preview artifact](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38373858558/artifacts/11799704062)
includes ferrocene, bis-periazulene, Lu2(Pc)3, Ni@Pb12 and cryptand examples.

The checked-in **55 fully executed notebooks / 168 SVGs** are still
authentic outputs from the earlier [full run 38362379662](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38362379662);
they have **not** been falsely relabeled as images from the current
generic renderer. A final manually dispatched complete kernel run and
full-version CI matrix will refresh those outputs after design approval.
