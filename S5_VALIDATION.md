# Section S5 corpus audit — 55 records, 84 electron-allocation forms

This is a **data-and-executable-witness** verification, not a new CCS formal specification. The normative definitions remain Supporting Information Sections S1–S4, and the displayed-index W(D) convention is from S5.

## Exact source and normalization

Input: the user-provided latest `supporting_information(1).pdf` (72 pages, S5 worked examples on PDF pages 27–67) and `supporting_information(1).docx`. The fixture contains SHA-256 checksums for both originals, the S5 subsection and record number, page/paragraph locators, the temporary displayed nucleus enumeration, the complete undirected connection set, and every full allocation Form A/B/C/D.

The PDF supplies SI S5's n_empty and N_e where Word equation objects do not appear in python-docx plain-text extraction. Form B/C/D changes are expanded relative to Form A, **not** cumulatively from the preceding Form. The sorted [i,j,ai,aj] representation uses i<j with allocations reoriented to match the displayed indices. These integers are **temporary addresses**, not canonical molecular identity fields.

PDF page-header text can divide a single printed connection/allocation token. Such headers were removed before parsing; complete extracted nucleus lists, connection sets, per-form assigned counts, and endpoint allocations then agreed with the separate DOCX extraction. No SI source value has been silently adjusted.

## Offline source audit (2026-10-13)

| Check | Result |
| --- | ---: |
| Record coverage against SI Table S4 | 55 / 55 |
| Form coverage against SI Table S4 | 84 / 84 |
| PDF versus DOCX extracted explicit nucleus/connection/allocation entries | 0 mismatches |
| SI Table S4 nucleus/connection/electron counts and charges | 0 mismatches |
| S1.4 endpoint-budget/well-formedness violations | 0 |
| Independent finite Nuc/Ele/Con replay witnesses (not using the package) | 84 / 84 pass |

There are 1,793 nucleus entries and 2,014 connection entries over the 55 example records, and 2,900 explicit endpoint-allocation connection records summed over the 84 forms. Nonzero unassigned-electron counts occur in S5 records 43, 46, and 48. Unusual, asymmetric and zero endpoint allocations are preserved literally.

## Machine-readable fixture

File: `data/s5_records.json.xz`. It contains UTF-8 JSON compressed with the Python-standard-library xz/lzma format to keep the private repository small. Read it with:

```python
import json
import lzma

with lzma.open("data/s5_records.json.xz", "rt", encoding="utf-8") as f:
    corpus = json.load(f)
print(len(corpus["records"]))  # 55
```

- `record["nuclei"]`: element identities in displayed temporary order.
- `record["connections"]`: sorted [i,j] pairs, preserving (0,0) connections.
- `record["expected"]`: Table S4 numbers (`nuclei`, `connections`, `electrons`, `charge`, `forms`).
- `record["forms"]`: each complete form's `assigned` vector, sorted `allocated_connections` as [i,j,ai,aj], and `n_empty`.
- `record["source"]`: original PDF page(s) and DOCX paragraph index (0-based).

## CI

`python -m pytest tests/test_s5_records.py` tests all 84 forms as CCS value objects, verifies exact S5 token order and replay, the electron-allocation quotient, Table S4 counts and charges, and representative intrinsic relabelings. This test uses only the existing CCS core and Python standard library (no RDKit or visualization dependency).

The original PDF and DOCX are not committed to the repository. The included hashes and source positions permit an external source audit; CI **does not** independently regenerate the fixture from those documents. Passing these tests establishes consistency and reproducible witnesses for the **specified** 84 forms, not physical chemical validity or exhaustive coverage of all finite CCS words.

**Status**: source-level checks pass; executable CI status is reported by the PR's latest workflow run. This PR must not change S1–S4 semantics or be merged without explicit approval.
