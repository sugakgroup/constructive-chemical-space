# Constructive Chemical Space (CCS)

Reference implementation for *Defining Chemical Space Constructively*.

CCS distinguishes **constitution**, **electron allocation**, and **construction**:

- Constitution: nuclear identities, direct connections and total electron count, `M = (C_X, S, N_e)`.
- Electron allocation: electrons assigned to nuclei, allocation at each connection endpoint, and unassigned electrons.
- Construction: finite words over `Nuc`, `Ele` and `Con`, with optional practical operations.

The authoritative formal definitions are in the manuscript's Supporting Information (SI Sections S1–S4). Worked examples and the fixed-index `W(D)` convention are in SI Section S5. Neither nucleus indices nor construction histories are intrinsic molecular identity fields.

## Install

Python **3.10 or later** is supported. From a checkout of this repository:

~~~bash
python -m pip install -e .
~~~

For drawing, RDKit conversion and local testing:

~~~bash
python -m pip install -e ".[dev]"
~~~

The core package has no mandatory third-party runtime dependencies. Visualization is optional; see the `visualization` extra in [pyproject.toml](pyproject.toml).

## Quick start

Construct a form from an explicit word, then obtain its constitution:

~~~python
from constructive_chemical_space import Nuc, Con, Word

word = Word([
    Nuc("O"), Nuc("H"), Nuc("H"),
    Con(0, 1, 1, 1),
    Con(0, 2, 1, 1),
])
form = word.replay()                     # ElectronAllocation
constitution = form.to_constitution()   # Constitution
assert constitution.num_electrons == 10
~~~

`Nuc(X, n=None, u=0)` installs a nucleus with `n` assigned electrons (default `Z(X)`), **not** a formal charge. `Con(i, j, ai, aj)` uses temporary nucleus addresses and contributes a direct connection if applicable; an inapplicable `Con` is a complete no-op. `Ele()` adds an unassigned electron.

Constitutional equality is index-independent:

~~~python
from constructive_chemical_space import Constitution

a = Constitution.from_indexed(["O", "H", "H"], [(0, 1), (0, 2)], 10)
b = Constitution.from_indexed(["H", "O", "H"], [(1, 0), (1, 2)], 10)
assert a == b
~~~

For any well-formed electron-allocation form `D`, `D.to_word(order=...)` constructs a replay witness. Supply an explicit sequence of the *same nuclei* when token-order reproducibility matters; omitting `order` gives a valid but not canonically indexed witness. `PracticalWord` supports finite macros and a terminal `End` operation for deterministic explicit-hydrogen completion. The [specification mapping](SPECIFICATION.md) records the exact v1.0 API and operation contracts.

## Molecular representations and visualization

CCS connection and allocation data are not conventional bond orders. The derived assignment imbalance `q_v = Z(X(v)) - n_v` is **not** automatically a Lewis formal charge or an oxidation state. RDKit conversion (`to_rdkit_mol`, `to_smiles`) is an explicitly non-normative, potentially lossy projection.

`Constitution`, `ElectronAllocation`, `Word` and `PracticalWord` support `show()`, `to_svg()`, `to_png()` and `save(path)` when visualization dependencies are installed.

~~~python
form.show()                                    # IPython/Jupyter
form.save("water-allocation.svg")
constitution.save("water-constitution.svg")
~~~

Diagrams retain **all nuclei, including H**, and every direct connection by default. Unusual endpoint allocations are labeled beside their own nuclei; group-13–17 residual valence electrons appear as vector dots, and negative assignment imbalances as vector strokes. Layout uses generic RDKit/CoordGen candidates and graph-based ring/cage quality criteria. Crossing-free diagrams are not enforced when this would degrade recognizable 2D chemical shapes. No 3D geometry, stereochemistry or physical electronic state is inferred. Omitting H from a depiction is an explicit **lossy** option.

## Reproducible SI S5 examples

The repository includes **55 individually executed Jupyter notebooks** covering **84 electron-allocation forms**, each containing explicit `Word([Nuc(...), Ele(), Con(...)])` source, actual Replay and tests against the S5 record. Every Form has an output for its Constitution and its electron allocation: **168 genuine saved SVG figures**.

- [Index of all 55 S5 notebooks](examples/s5/NOTEBOOK_INDEX.md)
- [Execution and regeneration instructions](examples/s5/README.md)
- [Source transcription, correction and validation audit](S5_VALIDATION.md)
- [Machine-readable source manifest](examples/s5/data/manifest.json) and [execution provenance](examples/s5/execution_report.json)

The last full execution was [GitHub Actions run 38391623509](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38391623509) (Python 3.11.17, RDKit 2026.09.1, nbclient 0.11.0). The saved notebook outputs and their SHA-256 checksums are preserved in this repository; they have not been fabricated. An optional, **manual read-only** workflow regenerates the full corpus without overwriting the committed notebooks.

**Editorial caveat:** The Record 13 (bis-periazulene) Forms B/C data include an author-confirmed correction. The final SI manuscript must independently incorporate that correction; the original PDF/DOCX source audit predates it. See [S5_VALIDATION.md](S5_VALIDATION.md).

## Testing and scope

~~~bash
python -m pytest -q
~~~

The test suite covers `tests/test_s1_formal.py` through `test_s4_practical.py`, RDKit conversion, visualization, all 84 S5 witnesses and saved-notebook provenance. CI runs S1–S4 core tests with Python **3.10, 3.11, 3.12**, plus the full regression suite on Python 3.11. The audited post-merge implementation passed [550 tests](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38393815432).

CCS does not infer the connection graph from physical structure, rank electron-allocation forms, or define stereochemical, conformational, spin or topological identities. The [specification mapping](SPECIFICATION.md) states the supported API and boundaries.

## License

MIT. See [LICENSE](LICENSE).
