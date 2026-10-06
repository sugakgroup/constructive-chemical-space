from pathlib import Path

from rdkit import Chem

from ccs_s1_s4 import Constitution, ElectronAllocation, Nuc, Con, End, PracticalWord, Word
from ccs_visualize import (
    draw_allocation,
    draw_constitution,
    draw_word,
    save_drawing,
    to_depiction_mol,
)


def test_constitution_view_is_graph_only_and_has_no_local_charge():
    m = Constitution.from_indexed(["C", "C"], [(0, 1)], 12)
    mol = to_depiction_mol(m, order=tuple(m.nuclei))
    assert mol.GetNumAtoms() == 2
    assert mol.GetNumBonds() == 1
    assert mol.GetBondWithIdx(0).GetBondType() == Chem.BondType.SINGLE
    assert all(atom.GetFormalCharge() == 0 for atom in mol.GetAtoms())


def test_allocation_view_shows_ccs_derived_formal_charge_and_literal_one_e_bond():
    d = ElectronAllocation.from_indexed(
        ["C", "C"],
        [6, 5],
        [(0, 1, 1, 0)],
    )
    mol = to_depiction_mol(d)
    assert sorted(atom.GetFormalCharge() for atom in mol.GetAtoms()) == [0, 1]
    bond = mol.GetBondWithIdx(0)
    assert bond.GetBondType() == Chem.BondType.ZERO
    assert bond.GetProp("bondNote") in {"1e (1,0)", "1e (0,1)"}


def test_standard_symmetric_allocations_use_conventional_glyphs():
    for k, expected in [
        (1, Chem.BondType.SINGLE),
        (2, Chem.BondType.DOUBLE),
        (3, Chem.BondType.TRIPLE),
    ]:
        d = ElectronAllocation.from_indexed(
            ["C", "C"],
            [6, 6],
            [(0, 1, k, k)],
        )
        mol = to_depiction_mol(d)
        assert mol.GetBondWithIdx(0).GetBondType() == expected
        assert not mol.GetBondWithIdx(0).HasProp("bondNote")


def test_high_or_unusual_allocations_are_literal_not_bond_order_inference():
    for pair in [(0, 0), (3, 2), (5, 5), (0, 6)]:
        ai, aj = pair
        d = ElectronAllocation.from_indexed(
            ["U", "O"],
            [92, 10],
            [(0, 1, ai, aj)],
        )
        mol = to_depiction_mol(d)
        bond = mol.GetBondWithIdx(0)
        assert bond.GetBondType() == Chem.BondType.ZERO
        assert bond.HasProp("bondNote")


def test_dative_allocation_is_directed_donor_to_acceptor():
    d = ElectronAllocation.from_indexed(
        ["N", "B"],
        [7, 5],
        [(0, 1, 2, 0)],
    )
    mol = to_depiction_mol(d, order=tuple(d.nuclei))
    bond = mol.GetBondWithIdx(0)
    assert bond.GetBondType() == Chem.BondType.DATIVE


def test_word_replay_preserves_construction_addresses_for_display_indices():
    w = Word([
        Nuc("O"),
        Nuc("H"),
        Con(0, 1, 1, 1),
    ])
    mol = to_depiction_mol(w, show_indices=True)
    assert mol.GetAtomWithIdx(0).GetProp("atomNote").startswith("0")
    assert mol.GetAtomWithIdx(1).GetProp("atomNote").startswith("1")


def test_practical_word_is_expanded_before_drawing():
    p = PracticalWord([Nuc("C"), End()])
    mol = to_depiction_mol(p, show_hydrogens=True)
    assert mol.GetNumAtoms() == 5  # CH4 with explicit H completion


def test_footer_contains_total_electron_count_and_n_empty_when_nonzero():
    d = ElectronAllocation.from_indexed(
        ["C"],
        [6],
        [],
        num_unassigned_electron=1,
    )
    svg = draw_allocation(d, fmt="svg", width=500, height=350)
    assert "N_e = 7" in svg
    assert "n_empty = 1" in svg


def test_constitution_footer_can_show_derived_total_charge():
    m = Constitution.from_indexed(["H", "H"], [(0, 1)], 1)
    svg = draw_constitution(
        m,
        fmt="svg",
        show_total_charge=True,
        width=500,
        height=350,
    )
    assert "N_e = 1" in svg
    assert "Q = +1" in svg


def test_save_drawing_svg_and_png(tmp_path: Path):
    d = ElectronAllocation.from_indexed(
        ["C", "C"],
        [6, 6],
        [(0, 1, 2, 2)],
    )
    svg_path = save_drawing(d, tmp_path / "ethene.svg")
    png_path = save_drawing(d, tmp_path / "ethene.png")
    assert svg_path.read_text(encoding="utf-8").startswith("<?xml")
    assert png_path.read_bytes().startswith(b"\x89PNG")


def test_draw_word_alias_runs():
    w = Word([Nuc("H"), Nuc("H"), Con(0, 1, 1, 1)])
    svg = draw_word(w, fmt="svg", width=500, height=350)
    assert "N_e = 2" in svg