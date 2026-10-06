# constructive-chemical-space
Reference implementation of Constructive Chemical Space (CCS): molecular constitutions, electron-allocation forms, finite construction replay, hydrogen completion, and construction macros.

## Visualization

The optional `ccs_visualize.py` layer renders CCS objects with RDKit without
changing CCS identity or replay semantics.

```python
from ccs_s1_s4 import ElectronAllocation, Nuc, Con, Word
from ccs_visualize import draw, save_drawing

d = ElectronAllocation.from_indexed(
    ["C", "C"],
    [6, 5],
    [(0, 1, 1, 0)],
)

save_drawing(
    d,
    "one_electron_bond.svg",
    show_indices=True,
    show_hydrogens=False,
)
```

`draw(...)` accepts a `Constitution`, `ElectronAllocation`, `Word`, or
`PracticalWord`. Words are replayed before drawing; practical words are
expanded first.

Depiction rules are deliberately conservative:

- Constitution view: direct connections are drawn as graph edges only; no bond
  order or local charge is inferred. `N_e` is shown in the footer.
- Allocation view: derived local charges `q_v = Z(X(v)) - n_v` are displayed.
- `(1,1)`, `(2,2)`, and `(3,3)` use single, double, and triple bond glyphs.
- `(2,0)` and `(0,2)` use a directed dative-bond glyph.
- Every other allocation is drawn as a connection with a literal label, e.g.
  `(0,0)`, `1e (1,0)`, `(3,2)`, `(5,5)`, `(0,6)`, or `(0,9)`.
- `n_empty` is shown in the footer when nonzero.
- `diagnostics=True` adds `n`, residual-electron count, and derived charge to
  atom annotations.

SVG and PNG output are supported. `show_indices=True` displays temporary
construction/display addresses; these indices are not intrinsic CCS identity.

