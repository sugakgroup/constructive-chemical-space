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

### Cluster schematics and portable electron-dot graphics

The optional 2D depiction uses a schematic separated-ligand arrangement for
multi-decker complexes, a metal-in-shell radial layout for Ni@Pb12 and three
visually separated N-to-N bridges for 2.2.2 cryptand hosts. These sketches
retain **all** CCS nuclear incidences but are not computed 3D geometries.
For genuinely nonplanar graphs, some line crossings remain mathematically
unavoidable; the figure retains the explicit warning. Use
`layout="cluster"` to force the specialized view (or `"auto"` to detect
the motif); `layout="coordgen"` / `"rdkit"` remain available to compare
against the generic layout.

Residual valence electron dots are rendered using actual SVG vector
`<circle>` primitives and equivalent PNG ellipses, not the font-dependent
UTF-8 `•` character. Per-nucleus residual counts and origin are retained as
numeric SVG data attributes. This avoids the `\xe2\x80\xa2` text-display
issue while preserving the exact CCS electron bookkeeping.

### CI usage

Ordinary PR #4 branch pushes now run a **single representative smoke**
(Python 3.11, eight S5 records including Lu2(Pc)3, Ni@Pb12, and both cryptand
examples). Their kernel outputs are intentionally ephemeral and are never
committed or called the complete validation. The **full 55-notebook execution
and saved-output refresh** and the full version matrix run only by manually
dispatching the corresponding GitHub Actions workflows at final review.

### Viewing cluster changes before the final full-corpus run

Representative smoke is also available as a short-lived GitHub Actions
artifact (`s5-smoke-clusters-<run id>`) containing actual-kernel SVG
constitution and allocation diagrams for eight selected S5 records,
including the 4 difficult coordination and cluster examples. The artifact
is not merged into the notebooks, expires after five days, and is **not**
claimed to replace the final 55-notebook execution and provenance report.
The all-record notebook regeneration remains a manual final-review action.

### Ligand-first / metal-last placement

For metallocenes and multi-decker complexes the organic ring structures are
laid out first. Metal atoms are then placed last using deterministic candidate
positions outside the crowded ligand nuclei, scored for atom clearance,
direct-connection interference and unnecessary distance. The metals remain
fully connected: no Fe–C, Tb–N or Lu–N incidence is removed. This is a
readability-driven schematic 2D convention, not a proposed 3D geometry.

For cryptands, the three bridging chains and **all explicit H positions**
are placed first. Two hydrogens on a parent atom are fanned out and checked
against nearby atoms; K/Na is added last, away from the ligand. Central Ni
inside Pb12 is deliberately retained as the exception where the metal is
intrinsically the included graph hub. Residual electron dots are kept nearer
their source atom than charge annotations.

### Planarity versus a distant metal: a deliberate limitation

For planar metal-containing graphs such as Fe(Cp)2, U(COT)2 and Tb(Pc)2,
putting the metal entirely outside otherwise regular ligand polygons
creates artificial straight-line crossings. In these cases, the renderer
honors the prior **zero-crossing** requirement by falling back to an
unintersected planar embedding and locally increasing metal clearance while
keeping the ligand vertices fixed. For genuinely nonplanar layered
complexes, the ligand-first / metal-last arrangement is retained and the
unavoidable 2D crossings are identified explicitly. The full CCS
connection set is always preserved.
