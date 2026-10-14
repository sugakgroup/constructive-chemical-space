# Section S5 corpus audit — 55 records, 84 electron-allocation forms

This is a **data-and-executable-witness** verification, not a new CCS formal specification. The normative definitions remain Supporting Information Sections S1–S4, and the displayed-index W(D) convention is from S5.

## Exact source and normalization

Input: the user-provided latest `supporting_information(1).pdf` (72 pages, S5 worked examples on PDF pages 27–67) and `supporting_information(1).docx`. The fixture contains SHA-256 checksums for both originals, the S5 subsection and record number, page/paragraph locators, the temporary displayed nucleus enumeration, the complete undirected connection set, and every full allocation Form A/B/C/D.

The PDF supplies SI S5's n_empty and N_e where Word equation objects do not appear in python-docx plain-text extraction. Form B/C/D changes are expanded relative to Form A, **not** cumulatively from the preceding Form. The sorted [i,j,ai,aj] representation uses i<j with allocations reoriented to match the displayed indices. These integers are **temporary addresses**, not canonical molecular identity fields.

PDF page-header text can divide a single printed connection/allocation token. Such headers were removed before parsing; complete extracted nucleus lists, connection sets, per-form assigned counts, and endpoint allocations then agreed with the separate DOCX extraction. No SI source value has been silently adjusted.

## Original PDF/DOCX extraction audit (2026-10-13; pre-Record-13 correction)

| Check | Result |
| --- | ---: |
| Record coverage against SI Table S4 | 55 / 55 |
| Form coverage against SI Table S4 | 84 / 84 |
| PDF versus DOCX extracted explicit nucleus/connection/allocation entries | 0 mismatches |
| SI Table S4 nucleus/connection/electron counts and charges | 0 mismatches |
| S1.4 endpoint-budget/well-formedness violations | 0 |
| Independent finite Nuc/Ele/Con replay witnesses (not using the package) | 84 / 84 pass |

There are 1,793 nucleus entries and 2,014 connection entries over the 55 example records, and 2,900 explicit endpoint-allocation connection records summed over the 84 forms. Nonzero unassigned-electron counts occur in S5 records 43, 46, and 48. Unusual, asymmetric and zero endpoint allocations are preserved literally.

### Source and fixture fingerprints

- SI PDF SHA-256: `777f7746048a1859aef9c955a0346063082a3ff91384c319dfc90b765bfa31b5`
- SI DOCX SHA-256: `00513db5f6d68fe456911460cbedbf185f77393d01e2642d081079a0c09cb89a`
- Original (pre-author-correction) fixture SHA-256: `1a9c1dcd119159aba4a95ff9941ea73110de18f97583200e33b3578df6aad8dd`
- Corrected fixture SHA-256: `0b10a14089114d4e0aca44a99f45ca68a78bb6bea5232d0ba5304b98ba98e426` (also verified from `manifest.json`)

The original fixture blob was uploaded byte-for-byte (Git blob SHA-1 `f594259a60a87756d223ec0183704ea1f02925ec`); this historical Git blob is not the current author-corrected fixture. The `n_empty=1` cases are records 43, 46, and 48, including the electron-only form on SI PDF page 61.

## Machine-readable fixture

File: `examples/s5/data/s5_records.json.xz`. It contains UTF-8 JSON compressed with the Python-standard-library xz/lzma format to keep the private repository small. Read it with:

```python
import json
import lzma

with lzma.open("examples/s5/data/s5_records.json.xz", "rt", encoding="utf-8") as f:
    corpus = json.load(f)
print(len(corpus["records"]))  # 55
```

- `record["nuclei"]`: element identities in displayed temporary order.
- `record["connections"]`: sorted [i,j] pairs, preserving (0,0) connections.
- `record["expected"]`: Table S4 numbers (`nuclei`, `connections`, `electrons`, `charge`, `forms`).
- `record["forms"]`: each complete form's `assigned` vector, sorted `allocated_connections` as [i,j,ai,aj], and `n_empty`.
- `record["source"]`: original PDF page(s) and DOCX paragraph index (0-based).

## CI

`python -m pytest tests/test_s5_examples.py` tests all 84 forms as CCS value objects, verifies exact S5 token order and replay, the electron-allocation quotient, Table S4 counts and charges, and representative intrinsic relabelings. This test uses only the existing CCS core and Python standard library (no RDKit or visualization dependency). For every form, the replay's temporary address table gives an explicit field-preserving bijection. Generic graph-isomorphism equality is additionally tested on smaller examples instead of repeatedly searching high-symmetry 100+-nucleus graphs; the full fieldwise replay checks establish intrinsic equality under S1 renaming.

The original PDF and DOCX are not committed to the repository. The included hashes and source positions permit an external source audit; CI **does not** independently regenerate the fixture from those documents. Passing these tests establishes consistency and reproducible witnesses for the **specified** 84 forms, not physical chemical validity or exhaustive coverage of all finite CCS words.

**Status**: Source-level fixture checks and post-merge executable CI passed. The final SI manuscript still requires synchronization of the explicitly author-corrected Record 13 Forms B/C; the archived original PDF/DOCX extraction audit covers the pre-correction version.

## Notebook organization and execution

The S5 examples are now exposed as **55 Jupyter notebooks** under
`examples/s5/notebooks/`, indexed by `examples/s5/data/manifest.json`.
A small `examples/s5/_common.py` reads the single compressed data fixture,
constructs the package's existing immutable CCS objects, verifies the exact
display-order `W(D)` and Replay, and supports optional SVG/PNG export. It is
not an alternative CCS core implementation.

Each notebook now contains independently inspectable **literal Nuc/Ele/Con
construction words** for every displayed Form. Its code actually replays
each Word, obtains a new ElectronAllocation and its Constitution, checks
the complete word tokens against the original SI fixture and validates
the constructed objects. The diagrams display **these replay-generated
objects**, not the reference fixture forms.

The actual notebook outputs **are committed**. The latest
[GitHub Actions run 38391623509](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38391623509)
regenerated all 55 notebooks after adopting the ring-aware layout and vector
negative-charge signs. The notebooks ran through real Jupyter ipykernels via
`nbclient.NotebookClient.execute()` (six code cells per Notebook). The
resulting genuine outputs were committed as
[5ef1bfb](https://github.com/sugakgroup/constructive-chemical-space/commit/5ef1bfb6278331c248d109155472e4a9109f1968).
Each Form A/B/C/D has its own complete, literal construction Word, Replay,
ElectronAllocation and Constitution. Its Constitution is compared with Form A
and the SI source under the complete explicit-index bijection, then **both
constructed objects are individually displayed**. The run captured **84
Constitution SVGs and 84 allocation SVGs, totaling 168 genuine images**,
rather than the earlier 139. RDKit SVG labels for unusual/asymmetric
allocations are independently anchored near each contributing nucleus, with
the endpoint atom index and numeric allocation preserved in SVG metadata.
Environment: Python 3.11.17, RDKit 2026.09.1, nbclient 0.11.0.
The run log and `examples/s5/execution_report.json` record the input
commit, execution counts, per-notebook code SHA-256 and SVG-output SHA-256.
`tests/test_s5_examples.py` rehashes these persisted kernel payloads and
compares them against the report. The output fields were never hand-built;
`generate_explicit_words.py` produces code only and deliberately resets
execution counts and output lists before real kernel execution.

Extra exported SVG/PNG files under `examples/s5/outputs/` remain opt-in and
Git-ignored. Complete explicit-H diagrams are the default.
CI checks every Notebook file's structure, correspondence to the manifest,
all 84 Form witnesses, selected notebook loading/replay cells executed as Python,
and **actual SVG rendering of all 84 word-generated constitutions and all 84
word-generated allocations** with explicit H retained. The original SI PDF/DOCX are not re-extracted in CI.

## Author-confirmed correction to Record 13 (bis-periazulene)

The author supplied the explicit Form B and Form C `a changes` and confirmed
the Form C `n changes`. The original machine-readable corpus incorrectly
contained `(9,14):(1,1)` and `(13,14):(2,2)` in both B/C. These two
assignments have now been corrected to `(9,14):(2,2)` and
`(13,14):(1,1)` in each form; the unchanged Form A values at `(10,11)`
and `(12,13)` remain `(2,2)`. Form C still contains `n_0=7`, `n_7=5`.

This is an **explicit author-specified correction**, distinct from the
original PDF–DOCX extraction audit, which applies to the pre-correction
fixture. It must be synchronized into the final submitted SI manuscript if
the source document differs. The new normalized binary fixture fingerprint
is stored in `examples/s5/data/manifest.json` and independently checked
by CI. The correction script is idempotent and checks the old values before
editing, rather than silently rewriting other source fields.

## Display scope: ring-aware 2D drawing and electron annotations

Every notebook depicts the CCS direct-connection graph with **all H nuclei
explicitly represented**. NetworkX classifies 51 S5 graphs as planar and
four (records 21, 42, 43 and 44) as mathematically nonplanar. Planarity is
a graph property; **a planar graph is not required to be drawn without
straight-line crossings**. The generic RDKit/CoordGen candidate score
prioritizes recognizable short rings and projected cages, reasonable
connection lengths, ring-centered vertices incident to an entire ring,
and nonbonded separation. Crossings have only a weak penalty. High-degree
graphs can use a temporary, degree-weakened CoordGen coordinate candidate,
without altering visible connections or the CCS form. A seeded spring
fallback is used only if both principal coordinate methods fail. No
element-, complex- or named-molecule-specific branch and no physical
3D geometry are introduced.

Group-13–17 residual valence-electron dots are drawn as actual SVG
`circle` primitives. Negative assignment-imbalance signs are constructed
from SVG `line` strokes, not Unicode U+2212 or an XML numeric entity
that IPython may decode to a raw UTF-8 minus in persisted notebook output.
Independent scanning of **all 168 genuine SVG outputs** in the latest
corpus confirmed 146 negative-charge annotations accompanied by 146 vector
minus strokes; no raw U+2212, numeric minus entity or visibly escaped
UTF-8 minus-byte sequence remained. Atomic `q_v = Z_v - n_v` and
electron-allocation endpoint counts are annotations, not Lewis formal
charges, physical electronic structure or constitutional identity fields.

The latest full reproduction is
[run 38391623509](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38391623509):
**55 notebooks, 84 forms, 168 actual SVG outputs**, with **494 tests
passing twice** (before and after committing the outputs). The complete
Python 3.10–3.12 core matrix is separately supported in
`.github/workflows/tests.yml`; a future manual real-kernel replay
can be launched through the read-only
`.github/workflows/execute_s5_notebooks.yml` workflow.
