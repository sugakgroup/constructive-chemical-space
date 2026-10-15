# Section S5 corpus audit — 55 records, 84 electron-allocation forms

This is a **data-and-executable-witness** verification, not a new CCS formal specification. The normative definitions remain Supporting Information Sections S1–S4, and the displayed-index W(D) convention is from S5.

## S5 dataset and normalization

The fixture contains all 55 S5 records and 84 electron-allocation forms. Each record includes its SI section and source locator, displayed temporary nucleus enumeration, complete direct connections, electron totals, and fully expanded Forms A/B/C/D.

Form B/C/D changes are interpreted relative to Form A, not cumulatively. Allocated connections are stored in sorted `[i,j,ai,aj]` order with `i<j` and endpoint allocations aligned to their displayed indices. These indices are temporary addresses, not intrinsic CCS identity fields.

## Validation summary

| Check | Result |
| --- | ---: |
| Record coverage against SI Table S4 | 55 / 55 |
| Form coverage against SI Table S4 | 84 / 84 |
| SI Table S4 nucleus/connection/electron counts and charges | 0 mismatches |
| S1.4 endpoint-budget/well-formedness violations | 0 |
| Independent finite Nuc/Ele/Con replay witnesses (not using the package) | 84 / 84 pass |

There are 1,793 nucleus entries and 2,014 connection entries over the 55 example records, and 2,900 explicit endpoint-allocation connection records summed over the 84 forms. Nonzero unassigned-electron counts occur in S5 records 43, 46, and 48. Unusual, asymmetric and zero endpoint allocations are preserved literally.

### Fixture integrity

The normalized S5 fixture has SHA-256
`6478ed51a58d5c0ccb4acdacee9335befedc3f28a4084bc0dcd96a8d9ffe6a43`,
recorded in `examples/s5/data/manifest.json` and verified by CI. The manifest
also retains source-document fingerprints from the original extraction for
traceability; these are not fingerprints for later manuscript revisions.

### Editorial correction to phenol Form C (2026-10-15)

Record 10, Form C has assigned-electron counts C(1):6 and C(3):7, instead of C(1):7 and C(3):6. Thus it represents the para-charge-separated phenol form shown in manuscript Figure 2. Form B retains the ortho assignment C(1):7. In Form C, O(6):7, all endpoint allocations, all connections, unassigned electrons, and the total of 50 electrons are unchanged. The corrected fixture is replayed and the corresponding notebook outputs are genuinely regenerated.

The source-document fingerprints identify the earlier SI extraction and are retained only as historical provenance. They are NOT checksums of the corrected SI PDF; the manuscript SI Form C text must be revised accordingly before submission. For this reason source.preserves_si is set to false in the corrected manifest.

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
- `record["source"]`: SI source page(s) and paragraph locator (0-based).

## CI

`python -m pytest tests/test_s5_examples.py` tests all 84 forms as CCS value objects, verifies exact S5 token order and replay, the electron-allocation quotient, Table S4 counts and charges, and representative intrinsic relabelings. This test uses only the existing CCS core and Python standard library (no RDKit or visualization dependency). For every form, the replay's temporary address table gives an explicit field-preserving bijection. Generic graph-isomorphism equality is additionally tested on smaller examples instead of repeatedly searching high-symmetry 100+-nucleus graphs; the full fieldwise replay checks establish intrinsic equality under S1 renaming.

Manuscript files are not part of the repository, and CI **does not** re-extract or independently compare revised manuscript documents. It checks the committed fixture's integrity and the 84 executable witnesses. Passing these tests does not establish physical chemical validity or exhaustive coverage of all finite CCS words.

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
the complete word tokens against the S5 fixture and validates
the constructed objects. The diagrams display **these replay-generated
objects**, not the reference fixture forms.

The actual notebook outputs **are committed**. The original corpus run was
[38391623509](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38391623509);
the corrected v1.0.0 corpus was rerun in [run 38704002752](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38704002752). The notebooks ran through real Jupyter ipykernels via
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
word-generated allocations** with explicit H retained. CI does not re-extract manuscript files.
