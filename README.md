# Constructive Chemical Space

Reference implementation of **Constructive Chemical Space (CCS)**, accompanying the manuscript *Defining Chemical Space Constructively*.

CCS separates three descriptive levels:

- **constitution**: nuclear identities, direct connections, and total electron count;
- **electron allocation**: assigned electron counts, endpoint allocations, and unassigned electrons;
- **construction**: finite ordered operations that build electron-allocation forms.

This repository implements the formal specification in Supporting Information Sections S1-S4. Nuclear indices are not molecular identity fields: integer addresses exist only while a construction is executed or while an indexed convenience input is parsed.

## Installation

For the core implementation and tests:

```bash
python -m pip install -e ".[test]"
```

For optional RDKit projection:

```bash
python -m pip install -e ".[dev]"
```

The package requires Python 3.10 or later.

## Public v1.0 API

The stable import surface is exposed from `constructive_chemical_space`:

```python
from constructive_chemical_space import (
    NuclearIdentity,
    Nucleus,
    NuclearCollection,
    Connection,
    AllocatedConnection,
    Constitution,
    ElectronAllocation,
    Nuc,
    Ele,
    Con,
    Word,
    End,
    Arg,
    MacroCon,
    Macro,
    PracticalWord,
    SearchProfile,
    RDKitConversionError,
)
```

Internal module paths are not part of the v1.0 API contract.

## Core replay

`Nuc(X; n, u)` installs one fresh nucleus. The second argument is the **assigned-electron count**, not a charge. When `n` is omitted, `Nuc(X)` uses the atomic number `Z(X)`. The optional `u` field is construction metadata used only by terminal `End` completion; changing `u` alone does not change ordinary core replay.

```python
from constructive_chemical_space import Nuc, Con, Word

water = Word([
    Nuc("O"),
    Nuc("H"),
    Con(0, 1, 1, 1),
    Nuc("H"),
    Con(0, 2, 1, 1),
])

D = water.replay()
M = D.to_constitution()
```

`Con(i, j; ai, aj)` uses temporary construction addresses. Negative addresses are relative to the current installation order. An inapplicable `Con` is a complete no-op: it does not partially modify the state or overwrite an existing connection.

## Index-independent objects

`Constitution` and `ElectronAllocation` are **immutable intrinsic value objects**. Their nuclear carriers, connections, assigned electron counts, and electron totals cannot be changed through attribute reassignment after construction. Assigned-electron inputs are defensively copied and stored read-only; `assigned_mapping()` returns a separate mutable copy. Convenience constructors may use temporary integer labels, but equality ignores arbitrary renaming of nuclei.

```python
from constructive_chemical_space import Constitution

water_a = Constitution.from_indexed(
    ["O", "H", "H"],
    [(0, 1), (0, 2)],
    10,
)

water_b = Constitution.from_indexed(
    ["H", "O", "H"],
    [(1, 0), (1, 2)],
    10,
)

assert water_a == water_b
```

## Constructive witnesses

A well-formed electron-allocation form can produce one finite witness word:

```python
from constructive_chemical_space import ElectronAllocation

D = ElectronAllocation.from_indexed(
    ["O", "H", "H"],
    [8, 1, 1],
    [(0, 1, 1, 1), (0, 2, 1, 1)],
)

# Any chosen enumeration gives a valid witness.  For reproducible token
# ordering across runs, preserve and pass an explicit sequence of Nucleus
# objects, e.g. a Word.replay_state().address_table from the same form.
w = D.to_word()
assert w.replay() == D
```

For **deterministic S5 construction words**, supply `D.to_word(order=nuclear_order)` where `nuclear_order` contains each nucleus of **that same D** exactly once. The resulting word installs nuclei in that order, emits `Ele` exactly `n_empty` times, and sorts `Con(i,j; ai,aj)` by increasing `(i,j)`, with `i < j` and allocations oriented to those endpoints. Without an explicit order, a valid but potentially different temporary enumeration is chosen. This does **not** define or require an intrinsic canonical numbering.

For example, when a form comes from replay, the accompanying temporary order can be retained explicitly:

```python
state = water.replay_state()
reproducible_word = state.form.to_word(order=state.address_table)
assert reproducible_word.replay() == state.form
```

The S1.4 auxiliary assignment imbalance `q_v = Z(X(v)) - n_v` is available through `D.assignment_imbalance(nucleus)`. The old `D.formal_charge(nucleus)` method remains as a **deprecated compatibility alias**, emitting `DeprecationWarning`; it does not assert a conventional Lewis formal charge or an oxidation state.

A constitution also provides the all-unassigned reference form used in the existence proof:

```python
D0 = water_a.all_unassigned_form()
assert D0.to_constitution() == water_a
```

## Terminal `End` completion

A practical input may contain `End`. The **first** `End` is terminal: it performs hydrogen completion once and all later tokens are ignored.

Under the S4 convention implemented here, automatic hydrogen completion applies only to nuclei whose underlying elements are in periodic-table groups 13-17. Hydrogen itself is excluded from automatic completion.

```python
from constructive_chemical_space import Nuc, Con, End, PracticalWord

ethene = PracticalWord([
    Nuc("C"),
    Nuc("C"),
    Con(0, 1, 2, 2),
    End(),
])

assert ethene.completion_counts() == (2, 2)
D = ethene.replay()
```

`PracticalWord.expand()` returns the explicit-H core `Word` generated by the completion rule.

## Finite macros

Macros are finite abbreviations of ordinary construction tokens; they do not introduce new molecular primitives.

External macro arguments are nonnegative addresses of nuclei that already exist when the macro is invoked. Inside a macro body, connection endpoints may use `Arg(k)` or negative relative addresses.

```python
from constructive_chemical_space import Arg, Macro, MacroCon, Nuc

add_carbon = Macro(
    "add_carbon",
    num_args=1,
    body=[
        Nuc("C"),
        MacroCon(Arg(0), -1, 1, 1),
    ],
)
```

If an external macro argument does not resolve to a pre-existing nucleus, the entire invocation contributes no core tokens.

## Optional search profiles

`SearchProfile` is only a **proposal policy**. It can restrict which candidate operations a search algorithm chooses to propose, but it does not change the construction alphabet, replay semantics, electron-allocation forms, `End`, or CCS constitutions.

## Optional RDKit projection

Electron-allocation forms can be projected to an RDKit `Mol` or, when possible, SMILES:

```python
mol = D.to_rdkit_mol()
smiles = D.to_smiles(remove_hs=True)
```

This projection is deliberately downstream and lossy. CCS endpoint allocations are not bond orders, and CCS assignment imbalance is **not** defined as RDKit/Lewis formal charge. The implementation uses a conservative heuristic for RDKit bond types and initially projects `q_v` into the toolkit formal-charge field solely as a best-effort guess (with a possible dative-charge retry during sanitization). The original CCS quantity is retained independently as the RDKit atom property `_CCSAssignmentImbalance`; the legacy property `_CCSAuxCharge` is retained for compatibility. CCS endpoint-allocation data are also preserved as RDKit bond properties, and `n_empty` is never silently localized. By default, a form with unassigned electrons returns `None` for RDKit/SMILES projection unless `allow_unassigned=True` is requested for a structural `Mol` carrying the count only as metadata.

## Visualization and export

The optional visualization API is implemented **inside the single implementation
module** `src/constructive_chemical_space/ccs.py`. The small package initializer
`__init__.py` only exposes imports; no second CCS core or visualization source
module is maintained.

Install the optional dependencies:

```bash
python -m pip install -e ".[visualization]"
```

`Constitution`, `ElectronAllocation`, `Word`, and `PracticalWord`
provide `.show()`, `.to_svg()`, `.to_png()`, and `.save(path)`.
`.show()` displays an SVG directly in Jupyter/IPython; `.to_svg()`
returns an SVG string, `.to_png()` returns PNG bytes, and
`.save("figure.svg")` / `.save("figure.png")` writes a file.

```python
from constructive_chemical_space import ElectronAllocation

D = ElectronAllocation.from_indexed(
    ["O", "H", "H"], [8, 1, 1],
    [(0, 1, 1, 1), (0, 2, 1, 1)],
)
D.show()                          # Jupyter / IPython
D.show(diagnostics=True, show_indices=True)
D.save("water-allocation.svg")    # same visualization without IPython
D.to_png(width=900, height=600)   # PNG bytes
D.to_constitution().show()        # graph-only view
```

**Faithful defaults:** Every explicit nucleus **including H**, and every direct
connection, is shown by default. In constitution view a plain stroke means a
direct CCS connection, not a conventional bond order. `N_e` is displayed in
the figure footer; optionally use `show_total_charge=True`.

Allocation view uses single/double/triple strokes for `(1,1)`, `(2,2)`,
`(3,3)`, and illustrative arrows for `(2,0)` / `(0,2)`.
All other endpoint allocations, including `(0,0)` and odd values, are
written **with the temporary index of each endpoint**, such as `0:3, 2:1`.
The SI S1.4 auxiliary assignment imbalance `q_v=Z(X(v))-n_v` is shown as
`q=+1` / `q=-1` atom notes where nonzero, **not** as an RDKit/Lewis formal
charge. `diagnostics=True` additionally displays `n`, residual electron
count `r`, and `q`.

`show_hydrogens=False` is an **explicit opt-in, lossy simplified depiction**:
it omits H nuclei *and their connections* from the picture, including bridging
or isotope-resolved hydrogen. Such images are marked `H omitted - lossy view`.
They must not be treated as complete CCS constitutions or as faithful molecular
identity diagrams.

For reproducible indexed diagrams use `order=(...)` with the exact nucleus
objects of the source form, or show a `Word` / `PracticalWord` (which preserves
replay address order). Intrinsic CCS objects do not have canonical nuclear
numbering. The visualization never mutates the CCS object or changes replay.

RDKit depiction is downstream and non-normative. In particular, line styles,
positions, and arrow glyphs are graphical interpretations only, not CCS identity
fields or physical bonding assignments.

## Tests

The tests are organized by the formal specification:

```text
tests/test_s1_formal.py
tests/test_s2_replay.py
tests/test_s3_coverage.py
tests/test_s4_practical.py
tests/test_rdkit_projection.py
tests/test_visualization.py
```

The core suite also checks intrinsic isomorphism, immutability and defensive copying, deterministic S5 word generation for explicit nuclear order, multi-residue glycine macros, End eligibility boundaries, and valid/invalid macro arguments. Complete machine-readable transcription and independent cross-checking of all **55 S5 records / 84 forms** is intentionally deferred to a separate data-validation change after v1.0 API freeze; it is not claimed by this implementation PR.

Run them with:

```bash
python -m pytest
```

## Scope

The reference implementation does **not** infer a unique CCS connection set from geometry or electronic structure, rank electron-allocation forms by physical relevance, or treat RDKit conversion as a lossless CCS serialization. Stereochemistry, conformational geometry, spin/electronic-state labels, and mechanical topology are outside the current CCS constitutional identity.

See [`SPECIFICATION.md`](SPECIFICATION.md) for the implementation-to-SI mapping and the v1.0 API boundary, and [`CONFORMANCE.md`](CONFORMANCE.md) for a requirement-to-test matrix and explicit limitations of the current evidence.

## License

MIT License. See [`LICENSE`](LICENSE).