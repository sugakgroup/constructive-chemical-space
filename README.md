# constructive-chemical-space
Reference implementation of Constructive Chemical Space (CCS): molecular constitutions, electron-allocation forms, finite construction replay, hydrogen completion, and construction macros.

## Visualization

The core CCS objects have an optional RDKit-backed `.show()` method for
Jupyter/IPython use.

```python
from ccs_s1_s4 import ElectronAllocation

d = ElectronAllocation.from_indexed(
    ["C", "C"],
    [6, 5],
    [(0, 1, 1, 0)],
)

d.show(show_indices=True)
```

The same method is available on `Constitution`, `Word`, and
`PracticalWord`. Words are replayed before display; practical words are
expanded first.

Depiction is deliberately conservative. Constitution view shows only the
direct-connection graph and `N_e`. Allocation view displays the derived local
charges `q_v = Z(X(v)) - n_v`; `(1,1)`, `(2,2)`, and `(3,3)` use
single/double/triple glyphs, `(2,0)` and `(0,2)` use a dative arrow, and
all other endpoint allocations are written literally on the connection
(e.g. `1e (1,0)`, `(3,2)`, `(5,5)`, `(0,6)`, `(0,9)`).

Hydrogens are hidden automatically when heavy atoms are present, but all-H
systems remain visible. Use `show_hydrogens=True` to force explicit H display,
`diagnostics=True` for `n`/residual/charge annotations, and
`show_total_charge=True` to add total charge to the footer.
