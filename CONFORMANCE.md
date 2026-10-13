# CCS v1.0 reference-implementation conformance record

Date: 2026-10-13. Basis: *Defining Chemical Space Constructively*, Supporting Information Sections S1-S4 and the Section S5 construction-word convention. The SI is normative; this file is an implementation/verification record, not a replacement specification.

## Requirement-to-test traceability

| SI requirement | Implementation | Regression tests |
| --- | --- | --- |
| S1.2-S1.3 intrinsic, index-independent constitutions and simple direct connections | `NuclearCollection`, `Constitution`, `_constitution_isomorphic` | `test_s1_formal.py` |
| S1.4 nonnegative assignments and endpoint budgets; zero-allocation edges | `ElectronAllocation`, `AllocatedConnection` | `test_s1_formal.py` |
| S1.4 assignment imbalance `q_v = Z(X(v)) - n_v` as auxiliary quantity | `ElectronAllocation.assignment_imbalance` | `test_s1_formal.py`, `test_rdkit_projection.py` |
| S1 intrinsic equality for relabeled, symmetric forms | `_constitution_isomorphic`, `_allocation_isomorphic` | `test_s1_formal.py` (including all 24 relabelings of a 4-nucleus form) |
| S1 value-object well-formedness after construction | immutable `Constitution`, immutable `ElectronAllocation`, defensive read-only assignment mapping | `test_s1_formal.py` |
| S2 admissible Nuc/Ele and total Con no-op replay | `Nuc`, `Ele`, `Con`, `ExecutionState`, `Word` | `test_s2_replay.py` |
| S3.2 exact constructability and all-unassigned lift | `to_word`, `all_unassigned_form` | `test_s3_coverage.py` |
| S4.1-S4.3 group 13-17 eligibility, terminal End, expansion into H core operations | `completion_plan`, `PracticalWord.expand` | `test_s4_practical.py` (Tables S1, S2 and selected S3 checks) |
| S4.5 finite macros, argument admissibility, repeated Gly | `Macro`, `MacroCon`, `Arg` | `test_s4_practical.py` |
| S5 deterministic `W(D)` relative to an explicit display enumeration | `ElectronAllocation.to_word(order=...)` | `test_s3_coverage.py` |
| Non-normative RDKit projection; `n_empty` is not silently localized | `to_rdkit_mol`, `to_smiles` | `test_rdkit_projection.py` |
| Optional S1/S5 visualization; full explicit-H graph, group-13–17 residual valence dots, atom-local q signs, and one generic density-aware 2D layout with zero crossings for planar graphs | `show`, `to_svg`, `to_png`, `save` on intrinsic objects and words | `test_visualization.py` (RDKit+IPython+Pillow CI) |
| S5 all 55 executed notebooks, 84 explicit W(D) words, full Form-specific D and M consistency, and 168 authentic SVG outputs | `examples/s5/data/manifest.json`, `examples/s5/data/s5_records.json.xz`, `ElectronAllocation.to_word(order=...)`, `Word.replay()` | `tests/test_s5_examples.py` |

## Reproducibility boundary

- The intrinsic `Constitution` and `ElectronAllocation` objects do **not** contain nuclear indices or canonical numbering.
- Repeatable token sequences are guaranteed for a *fixed explicit sequence of the same nucleus objects*, passed as `to_word(order=...)`; no guarantee is made for token identity with omitted `order`.
- RDKit atom formal charges, bond types, and SMILES are downstream heuristics, not CCS identity fields. `_CCSAssignmentImbalance` preserves the SI-derived q-value independently; `_CCSAuxCharge` is an older compatibility property.
- `formal_charge(nucleus)` remains a deprecated alias for `assignment_imbalance(nucleus)`, not a claim of chemical formal-charge equivalence.
- The complete S5 **55 records / 84 forms** are indexed by 55 Jupyter notebooks under `examples/s5/notebooks/` and transcribed in `examples/s5/data/s5_records.json.xz`. `tests/test_s5_examples.py` checks each full form, lexicographic explicit-order W(D), and Replay identity. The PDF/DOCX cross-extraction audit and source hashes are recorded in `S5_VALIDATION.md`; CI does not itself regenerate the fixture from the source documents.
- Passing executable regression tests supplies evidence for the specified cases; it is not a formal proof that all possible finite words or all physical systems have been exhaustively tested.

## Running validation

```bash
python -m pip install -e ".[dev]"
python -m pytest
```

During PR #4 development, `.github/workflows/tests.yml` runs one Python 3.11 representative smoke suite on push, including eight real-kernel notebooks. The complete Python 3.10/3.11/3.12 matrix and regeneration of all 55 authentic executed notebooks are separate **manual final-review** workflows. Earlier successful runs must never be presented as proof of the latest HEAD or the unrefreshed notebook outputs.
