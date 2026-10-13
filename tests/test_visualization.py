"""Optional visualization regression tests using installed CCS package objects."""
from io import BytesIO
import xml.etree.ElementTree as ET
from unittest.mock import patch

import pytest

pytest.importorskip("rdkit")
from rdkit import Chem
from PIL import Image

from constructive_chemical_space import (
    AllocatedConnection, Con, Constitution, Ele, ElectronAllocation,
    End, Nuc, Nucleus, PracticalWord, Word,
)
from constructive_chemical_space.ccs import (
    _build_ccs_visual_mol, _ccs_footer, _visual_source,
)


def _mol(obj, **kwargs):
    source, order, kind = _visual_source(obj, view=kwargs.pop("view", "auto"),
                                        order=kwargs.pop("order", None))
    return _build_ccs_visual_mol(
        source, order, source_kind=kind,
        show_hydrogens=kwargs.pop("show_hydrogens", True),
        show_indices=kwargs.pop("show_indices", False),
        diagnostics=kwargs.pop("diagnostics", False),
        label_all_allocations=kwargs.pop("label_all_allocations", False),
    )


def test_all_h_nuclei_connections_retained_in_constitution_default():
    m = Constitution.from_indexed(
        ["O", "H", "H"], [(0, 1), (0, 2)], 10
    )
    mol = _mol(m)
    assert mol.GetNumAtoms() == 3
    assert mol.GetNumBonds() == 2
    assert mol.GetProp("_CCSView") == "constitution"
    assert not mol.GetBoolProp("_CCSLossyHydrogenOmission")
    assert all(b.GetBondType() == Chem.BondType.SINGLE for b in mol.GetBonds())
    assert "N_e = 10" in m.to_svg()
    assert m.num_electrons == 10


def test_bridging_h_bifluoride_is_not_removed_by_default():
    d = ElectronAllocation.from_indexed(
        ["F", "H", "F"], [9, 2, 9],
        [(0, 1, 1, 1), (1, 2, 1, 1)],
    )
    for view in ("constitution", "allocation"):
        mol = _mol(d, view=view)
        assert (mol.GetNumAtoms(), mol.GetNumBonds()) == (3, 2)
        h = next(atom for atom in mol.GetAtoms() if atom.GetSymbol() == "H")
        assert h.GetDegree() == 2
    assert not d.to_svg().count("lossy view")


def test_hidden_h_is_explicitly_lossy_and_retains_full_system_footer():
    m = Constitution.from_indexed(["F", "H", "F"], [(0, 1), (1, 2)], 20)
    mol = _mol(m, show_hydrogens=False)
    assert mol.GetNumAtoms() == 2
    assert mol.GetNumBonds() == 0
    assert mol.GetBoolProp("_CCSLossyHydrogenOmission")
    svg = m.to_svg(show_hydrogens=False)
    assert "H omitted - lossy view" in svg
    assert "N_e = 20" in svg
    assert "H omitted - lossy view" not in m.to_svg()


def test_isotope_resolved_h_is_shown_explicitly():
    d = ElectronAllocation.from_indexed(
        ["2H", "O", "H"], [1, 8, 1],
        [(0, 1, 1, 1), (1, 2, 1, 1)],
    )
    mol = _mol(d)
    assert sorted((a.GetSymbol(), a.GetIsotope()) for a in mol.GetAtoms()) == [
        ("H", 0), ("H", 2), ("O", 0)
    ]
    assert mol.GetNumBonds() == 2


def test_assignment_imbalance_annotation_not_rdkit_formal_charge():
    d = ElectronAllocation.from_indexed(
        ["O", "H"], [7, 1], [(0, 1, 1, 1)]
    )
    mol = _mol(d)
    o = next(a for a in mol.GetAtoms() if a.GetSymbol() == "O")
    assert o.GetFormalCharge() == 0
    assert o.GetIntProp("_CCSAssignmentImbalance") == 1
    assert "q=+1" in o.GetProp("atomNote")
    assert "q=+1" in d.to_svg(show_indices=True)
    assert d.to_constitution().total_charge == 1


def test_diagnostics_show_primary_counts_and_residuals():
    d = ElectronAllocation.from_indexed(
        ["N", "B"], [7, 5], [(0, 1, 2, 0)]
    )
    mol = _mol(d, diagnostics=True)
    assert all(a.HasProp("atomNote") for a in mol.GetAtoms())
    assert all("n=" in a.GetProp("atomNote") for a in mol.GetAtoms())


@pytest.mark.parametrize(
    ("a", "b", "bond_type"), [
        (1, 1, Chem.BondType.SINGLE),
        (2, 2, Chem.BondType.DOUBLE),
        (3, 3, Chem.BondType.TRIPLE),
        (2, 0, Chem.BondType.DATIVE),
        (0, 2, Chem.BondType.DATIVE),
    ],
)
def test_standard_allocation_glyphs(a, b, bond_type):
    d = ElectronAllocation.from_indexed(
        ["U", "O"], [92, 10], [(0, 1, a, b)]
    )
    mol = _mol(d)
    assert mol.GetBondWithIdx(0).GetBondType() == bond_type


@pytest.mark.parametrize(("a", "b"), [(0, 0), (1, 0), (3, 2), (5, 5), (0, 6), (0, 9)])
def test_unusual_allocations_have_oriented_indexed_labels(a, b):
    u, o = Nucleus("U"), Nucleus("O")
    d = ElectronAllocation(
        [u, o], {u: 92, o: 10},
        [AllocatedConnection(u, o, a, b)]
    )
    mol = _mol(d, order=(o, u))
    bond = mol.GetBondWithIdx(0)
    assert bond.GetBondType() == Chem.BondType.ZERO
    assert bond.GetProp("bondNote") == f"0:{b}, 1:{a}"


def test_zero_allocation_edge_not_deleted():
    d = ElectronAllocation.from_indexed(
        ["C", "C"], [6, 6], [(0, 1, 0, 0)]
    )
    mol = _mol(d)
    assert mol.GetNumBonds() == 1
    assert mol.GetBondWithIdx(0).HasProp("bondNote")


def test_word_and_practical_word_preserve_replay_order():
    w = Word([Nuc("O"), Nuc("H"), Con(0, 1, 1, 1)])
    state = w.replay_state()
    # Replay creates new abstract Nucleus objects; compare temporary order
    # by nuclear identities rather than Python object identities.
    assert tuple(v.identity for v in _visual_source(w)[1]) == tuple(
        v.identity for v in state.address_table
    )
    p = PracticalWord([Nuc("C"), End()])
    src, order, kind = _visual_source(p)
    assert len(src.nuclei) == len(order) == 5
    assert kind == "allocation"


def test_disconnected_and_electron_only_can_be_rendered():
    d = Word([Ele()]).replay()
    mol = _mol(d)
    assert mol.GetNumAtoms() == 0
    svg = d.to_svg()
    assert "N_e = 1" in svg
    assert "n_empty = 1" in svg
    m = Constitution.from_indexed(["Na", "Cl"], [], 28)
    mol2 = _mol(m)
    assert (mol2.GetNumAtoms(), mol2.GetNumBonds()) == (2, 0)


def test_exotic_multicenter_allocations_have_literal_annotations():
    # A representative high-electron metal/cage allocation.
    d = ElectronAllocation.from_indexed(
        ["Ni", "Pb", "Pb"], [28, 82, 82],
        [(0, 1, 0, 7), (0, 2, 0, 8), (1, 2, 0, 0)],
    )
    mol = _mol(d)
    assert (mol.GetNumAtoms(), mol.GetNumBonds()) == (3, 3)
    assert all(b.GetBondType() == Chem.BondType.ZERO for b in mol.GetBonds())
    assert all(b.HasProp("bondNote") for b in mol.GetBonds())


def test_svg_png_export_and_save(tmp_path):
    p = PracticalWord([Nuc("C"), End()])
    svg = p.to_svg(width=600, height=450)
    ET.fromstring(svg)  # Well-formed SVG.
    assert "N_e = 10" in svg
    png = p.to_png(width=600, height=450)
    with Image.open(BytesIO(png)) as image:
        assert image.format == "PNG"
        assert image.size == (600, 450)
    svg_path = tmp_path / "methane.svg"
    png_path = tmp_path / "methane.png"
    p.save(svg_path, width=600, height=450)
    p.save(png_path, width=600, height=450)
    assert svg_path.read_text(encoding="utf-8") == svg
    assert p.to_svg(width=600, height=450) == svg
    assert png_path.read_bytes().startswith(b"\x89PNG")
    with pytest.raises(ValueError):
        p.save(tmp_path / "not-supported.pdf")


def test_show_calls_ipython_display_without_changing_core_data():
    pytest.importorskip("IPython")
    d = ElectronAllocation.from_indexed(["H", "H"], [1, 1], [(0, 1, 1, 1)])
    baseline = d.to_word().replay()
    with patch("IPython.display.display") as mocked:
        assert d.show(width=600, height=450) is None
        assert mocked.called
    assert d == baseline


def test_explicit_order_is_preserved_and_validated():
    a, b = Nucleus("C"), Nucleus("N")
    d = ElectronAllocation(
        [a, b], {a: 6, b: 7},
        [AllocatedConnection(a, b, 1, 1)]
    )
    assert _visual_source(d, order=(b, a))[1] == (b, a)
    with pytest.raises(ValueError):
        d.to_svg(order=(a, a))
