from unittest.mock import patch

from rdkit import Chem

from ccs_s1_s4 import (
    Constitution,
    ElectronAllocation,
    Nuc,
    Con,
    End,
    PracticalWord,
    Word,
    _build_show_mol,
    _render_show_svg,
    _show_source,
)


def test_constitution_show_is_graph_only():
    m = Constitution.from_indexed(["C", "C"], [(0, 1)], 12)
    source, order, kind = _show_source(m, view="constitution")
    mol = _build_show_mol(
        source, order,
        source_kind=kind,
        show_hydrogens=True,
        show_indices=False,
        diagnostics=False,
        label_all_allocations=False,
    )
    assert mol.GetBondWithIdx(0).GetBondType() == Chem.BondType.SINGLE
    assert all(a.GetFormalCharge() == 0 for a in mol.GetAtoms())


def test_allocation_show_charge_and_one_electron_label():
    d = ElectronAllocation.from_indexed(
        ["C", "C"], [6, 5], [(0, 1, 1, 0)]
    )
    source, order, kind = _show_source(d, view="allocation")
    mol = _build_show_mol(
        source, order,
        source_kind=kind,
        show_hydrogens=True,
        show_indices=False,
        diagnostics=False,
        label_all_allocations=False,
    )
    assert sorted(a.GetFormalCharge() for a in mol.GetAtoms()) == [0, 1]
    bond = mol.GetBondWithIdx(0)
    assert bond.GetBondType() == Chem.BondType.ZERO
    assert "1e" in bond.GetProp("bondNote")


def test_standard_allocations_use_standard_glyphs():
    for k, expected in [
        (1, Chem.BondType.SINGLE),
        (2, Chem.BondType.DOUBLE),
        (3, Chem.BondType.TRIPLE),
    ]:
        d = ElectronAllocation.from_indexed(
            ["C", "C"], [6, 6], [(0, 1, k, k)]
        )
        source, order, kind = _show_source(d, view="allocation")
        mol = _build_show_mol(
            source, order,
            source_kind=kind,
            show_hydrogens=True,
            show_indices=False,
            diagnostics=False,
            label_all_allocations=False,
        )
        assert mol.GetBondWithIdx(0).GetBondType() == expected


def test_unusual_allocations_are_literal():
    for ai, aj in [(0, 0), (3, 2), (5, 5), (0, 6), (0, 9)]:
        d = ElectronAllocation.from_indexed(
            ["U", "O"], [92, 10], [(0, 1, ai, aj)]
        )
        source, order, kind = _show_source(d, view="allocation")
        mol = _build_show_mol(
            source, order,
            source_kind=kind,
            show_hydrogens=True,
            show_indices=False,
            diagnostics=False,
            label_all_allocations=False,
        )
        bond = mol.GetBondWithIdx(0)
        assert bond.GetBondType() == Chem.BondType.ZERO
        assert bond.HasProp("bondNote")


def test_word_show_preserves_address_order():
    w = Word([Nuc("O"), Nuc("H"), Con(0, 1, 1, 1)])
    source, order, kind = _show_source(w, view="allocation")
    assert order == w.replay_state().address_table


def test_practical_word_show_expands_end():
    p = PracticalWord([Nuc("C"), End()])
    source, order, kind = _show_source(p, view="allocation")
    assert len(order) == 5


def test_footer_contains_electron_count_and_unassigned_count():
    d = ElectronAllocation.from_indexed(
        ["C"], [6], [], num_unassigned_electron=1
    )
    svg = _render_show_svg(d, width=500, height=350)
    assert "N_e = 7" in svg
    assert "n_empty = 1" in svg


def test_show_method_calls_ipython_display():
    d = ElectronAllocation.from_indexed(["H", "H"], [1, 1], [(0, 1, 1, 1)])
    with patch("IPython.display.display") as mocked:
        d.show(width=500, height=350)
        assert mocked.called
