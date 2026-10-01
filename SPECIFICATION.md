# CCS v1.0 specification mapping

This document records the implementation contract for the private v1.0 preparation branch/repository. The authoritative formal basis is Supporting Information **Sections S1-S4** of *Defining Chemical Space Constructively*.

## S1 — Scope and formal objects

### Nuclear collection

- `Nucleus` represents one abstract member of the finite carrier `V`.
- Two nuclei with the same `NuclearIdentity` remain distinct carrier elements.
- `NuclearCollection` stores no numerical nuclear index.
- `Constitution.__eq__` and `ElectronAllocation.__eq__` compare intrinsic objects up to a bijective renaming that preserves all retained fields.

### Constitution

`Constitution` implements

`M = (C_X, S, N_e)`

with a finite nuclear collection, simple undirected connection set, and nonnegative total electron count. Connections carry no bond multiplicity or donor direction at the constitutional level.

### Electron allocation

`ElectronAllocation` implements assigned counts, endpoint allocations, and `n_empty`. It enforces the well-formedness inequality at every nucleus. Residual counts and the assignment imbalance `q_v = Z(X(v)) - n_v` are derived, not independently stored.

`ElectronAllocation.to_constitution()` implements the electron-allocation quotient.

## S2 — Primitive construction and total replay

The public core tokens are:

- `Nuc(X, n=None, u=0)`
- `Ele()`
- `Con(i, j, ai, aj)`

`Nuc(X)` means `n = Z(X)`. The optional `u` value is retained only in the active execution state for S4 completion and has no effect on ordinary core replay.

`ExecutionState` contains the intrinsic allocation form plus temporary construction metadata. It is intentionally not part of the top-level stable API.

`Word.replay()` performs deterministic left-to-right replay. `Con` address resolution follows the SI definition for nonnegative and negative relative addresses. Every failed applicability check makes the entire `Con` a no-op.

## S3 — Constructive coverage helpers

`ElectronAllocation.to_word()` implements one exact finite construction witness after choosing a temporary enumeration of the finite carrier.

`Constitution.all_unassigned_form()` implements the all-unassigned reference lift used to prove that every constitution has at least one well-formed form.

These methods are constructive helpers; they do not define unique serializations or preferred electronic descriptions.

## S4 — Practical molecular construction

`PracticalWord` accepts core tokens, finite macro invocations, and `End`.

The first `End`:

1. evaluates the pre-`End` core state,
2. computes the S4 completion quantities,
3. appends explicit neutral H nuclei and `(1,1)` parent-H connections,
4. terminates execution,
5. ignores the entire suffix after that first `End`.

The automatic H-completion domain in S4 is the set of nuclear identities whose underlying elements are in groups 13-17. **Hydrogen itself is excluded.** Nuclei outside the domain, and in-domain nuclei failing the S4 eligibility conditions, receive zero automatically generated H without invalidating the preceding core replay.

`PracticalWord.expand()` returns the explicit-H core `Word`.

### Macros

`Macro` is a finite abbreviation of ordinary core tokens. External macro arguments are nonnegative addresses of nuclei present before invocation. `Arg(k)` represents such an argument inside a macro body; literal nonnegative addresses are not allowed in the body, while negative relative addresses are allowed. An invocation with any invalid external argument contributes zero tokens.

### SearchProfile

`SearchProfile` implements the optional S4 search-policy concept only. It filters proposals and does not alter formal replay semantics.

## RDKit projection is non-normative

RDKit support is outside S1-S4 and is not part of the CCS identity definition. `to_rdkit_mol()` and `to_smiles()` are best-effort downstream projections.

The current bond-type heuristic is:

- `(1,1)` -> single
- `(2,2)` -> double
- `(3,3)` -> triple
- higher equal allocations -> matching RDKit multiplicity when supported
- `(2,0)` / `(0,2)` -> directed dative representation
- `(0,0)` and other unsupported patterns -> unspecified bond

Original CCS allocation data are retained as RDKit properties where possible. A nonzero `n_empty` is not silently assigned to a nucleus.

## Stable v1.0 top-level imports

The supported top-level names are exactly those listed in `constructive_chemical_space.__all__`:

- `NuclearIdentity`
- `Nucleus`
- `NuclearCollection`
- `Connection`
- `AllocatedConnection`
- `Constitution`
- `ElectronAllocation`
- `Nuc`
- `Ele`
- `Con`
- `Word`
- `End`
- `Arg`
- `MacroCon`
- `Macro`
- `PracticalWord`
- `SearchProfile`
- `RDKitConversionError`

Internal helpers, execution-state classes, and module paths are not part of the v1.0 compatibility promise.

## Explicitly excluded from the v1.0 contract

The following are not required by S1-S4 and are not part of the v1.0 reference API:

- automatic inference of `S` from geometry or electronic structure;
- a unique canonical serialization of CCS objects;
- physical ranking of electron-allocation forms;
- stereochemistry, conformation, spin/electronic-state labels, or mechanical topology;
- molecule-generation/search algorithms;
- an automatic unique inverse mapping from RDKit/SMILES to CCS.