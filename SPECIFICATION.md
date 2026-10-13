# CCS v1.0 specification mapping

This document records the implementation contract for the private v1.0 preparation branch/repository. The authoritative formal basis is Supporting Information **Sections S1-S4** of *Defining Chemical Space Constructively*.

## S1 — Scope and formal objects

### Nuclear collection

- `Nucleus` represents one abstract member of the finite carrier `V`.
- Two nuclei with the same `NuclearIdentity` remain distinct carrier elements.
- `NuclearCollection` stores no numerical nuclear index.
- `Constitution.__eq__` and `ElectronAllocation.__eq__` compare intrinsic objects up to a bijective renaming that preserves all retained fields.
- Both are frozen, unhashable value objects; constructor validation establishes well-formedness and ordinary attribute assignment cannot invalidate it. `ElectronAllocation` copies assigned-electron inputs into a read-only mapping, while `assigned_mapping()` returns a copy.

### Constitution

`Constitution` implements

`M = (C_X, S, N_e)`

with a finite nuclear collection, simple undirected connection set, and nonnegative total electron count. Connections carry no bond multiplicity or donor direction at the constitutional level.

### Electron allocation

`ElectronAllocation` implements assigned counts, endpoint allocations, and `n_empty`. It enforces the well-formedness inequality at every nucleus. Residual counts and the S1.4 auxiliary **assignment imbalance** `q_v = Z(X(v)) - n_v` are derived, not independently stored. The public API is `assignment_imbalance(nucleus)`; `formal_charge(nucleus)` is retained only as a deprecated compatibility alias and emits `DeprecationWarning`. Neither name defines an RDKit/Lewis formal charge or oxidation state.

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

`ElectronAllocation.to_word(order=...)` implements the exact S5 word convention relative to one explicit temporary enumeration: Nuc in enumeration order, `n_empty` Ele tokens, and lexicographically ordered Con tokens with `i < j`, with endpoint allocations aligned to ascending addresses. The resulting Replay is intrinsically equal to the input form. Omitting `order` chooses an arbitrary enumeration that is valid for S3.2 but not promised to reproduce identical token sequences across runs.

`Constitution.all_unassigned_form()` implements the all-unassigned reference lift used to prove that every constitution has at least one well-formed form.

These methods are constructive helpers; they do not define unique serializations or preferred electronic descriptions. S5 `W(D)` is deterministic **for a fixed display enumeration**, not a canonical numbering of an intrinsic object.

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

Original CCS allocation data are retained as RDKit properties where possible. The initial RDKit atom formal-charge guess is taken from CCS `assignment_imbalance()` only as a non-normative conversion heuristic, not as a claim of conventional formal-charge equality. The independent CCS quantity is retained on each RDKit atom as `_CCSAssignmentImbalance`; the older `_CCSAuxCharge` property remains as a compatibility alias. Sanitization may retry with different toolkit formal charges. A nonzero `n_empty` is not silently assigned to a nucleus.

## Optional visualization (single implementation module)

Visualization is implemented in `src/constructive_chemical_space/ccs.py`,
without a separate source module and without modifying the SI S1–S4 definition
of Constitution, ElectronAllocation, construction words, or replay. The
package initializer only reexports public names.

The supported display and export methods for `Constitution`,
`ElectronAllocation`, `Word`, and `PracticalWord` are:
`show(**kwargs)` (IPython SVG display), `to_svg(**kwargs)` (SVG string),
`to_png(**kwargs)` (PNG bytes), and `save(path, **kwargs)` (SVG or PNG).
The optional exception `CCSVisualizationError` identifies failures of
optional rendering/display dependencies.

By default, every nucleus (including H) and every direct connection is drawn.
`show_hydrogens=False` deliberately creates a lossy view and the exported
figure warns that H was omitted. In constitution view a stroke represents an
S1 direct connection only, without implying bond multiplicity.

In allocation view, standard graphical strokes are used only as a
non-normative aid for selected endpoint allocations. Nonstandard endpoint
electron counts are positioned separately next to their contributing nuclei,
using the coordinates of the actual drawn connection endpoints
(`allocation_label_mode="endpoint"`, the default). SVG retains explicit
`data-atom` and `data-electrons` metadata for each individual label,
and PNG uses the same endpoint mapping. For backward-compatible central
bond annotations use `allocation_label_mode="bond"`; to suppress labels use
`"none"`. `label_all_allocations=True` additionally shows endpoint
counts on conventional single/double/triple/dative glyphs. Neither labels
nor line styles enter CCS intrinsic identity.
The S1.4 assignment imbalance `q_v=Z_v-n_v` is represented by a compact
atom-local sign (`+`, `−`, `2−`, etc.) and retained as explicit CCS
metadata, **not** as an RDKit/Lewis formal charge. Group-13–17 atoms additionally
display the nonnegative valence residual count
`n_v - Σ_s a_{v,s} - (Z_v-(group_v-10))` as `•` glyphs; this is a formal
accounting count and does **not** claim a particular spin state or a physically
unique lone-pair interpretation. H and nuclei outside that p-block range
receive no residual valence dots by default. The visualization input is immutable
and does not acquire any display indices; caller-provided `order` is used
only for reproducible depiction. Word-based displays use replay order by
default. SVG generation requires RDKit and NetworkX for planarity testing; PNG also
requires Pillow, and `show()` additionally requires IPython. The 2D layout
uses no metal, cluster, ligand or named-example-specific branches: RDKit and
CoordGen candidate drawings are compared with a seeded spring-repulsion
candidate for high-degree/nonplanar graphs using the same intersection,
interatomic-spacing and bond-length criteria. A crossing-free straight-line
planar embedding is mandatory where the direct-connection graph is planar.
Planar CCS graphs must have no unrelated bond-line crossings. Mathematically nonplanar
graphs cannot satisfy that condition and receive a visible disclaimer rather
than pretending the crossings can be removed. No physical 3D model is inferred.

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
- `CCSVisualizationError`

Internal helpers, execution-state classes, and module paths are not part of the v1.0 compatibility promise.

### Stable v1.0 entry points

The principal stable constructors and methods are:

- `Constitution.from_indexed(...)`, `Constitution.all_unassigned_form()`
- `ElectronAllocation.from_indexed(...)`, `ElectronAllocation.to_constitution()`, `ElectronAllocation.to_word(order=None)`, `ElectronAllocation.assignment_imbalance(nucleus)`
- Deprecated compatibility: `ElectronAllocation.formal_charge(nucleus)` (returns `assignment_imbalance()` and emits `DeprecationWarning`)
- `ElectronAllocation.to_rdkit_mol(...)`, `ElectronAllocation.to_smiles(...)`
- `Word.replay()`, `Word.replay_state()`, `Word.to_constitution()`
- `PracticalWord.pre_end_state()`, `PracticalWord.completion_plan()`, `PracticalWord.completion_counts()`
- `PracticalWord.expand()`, `PracticalWord.replay()`, `PracticalWord.to_constitution()`
- `PracticalWord.to_rdkit_mol(...)`, `PracticalWord.to_smiles(...)`
- `SearchProfile.allows(...)`

For `ElectronAllocation`, the public unassigned-electron field/keyword is `num_unassigned_electrons`.

## Scope of S5 example verification

The original PR #1 covers the **S5 deterministic construction-word convention**. PR #4 adds the complete machine-readable corpus of **55 S5 records and 84 electron-allocation forms**, individual explicit-word Jupyter notebooks, replay verification, and source-provenance checks. These are executable witnesses for the specified examples, not an assertion about physical geometry or exhaustive formal chemical space.

## Explicitly excluded from the v1.0 contract

The following are not required by S1-S4 and are not part of the v1.0 reference API:

- automatic inference of `S` from geometry or electronic structure;
- a unique canonical serialization of CCS objects;
- physical ranking of electron-allocation forms;
- stereochemistry, conformation, spin/electronic-state labels, or mechanical topology;
- molecule-generation/search algorithms;
- an automatic unique inverse mapping from RDKit/SMILES to CCS.