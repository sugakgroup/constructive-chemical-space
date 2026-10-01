import pytest

pytest.importorskip("rdkit")

from constructive_chemical_space import Con, ElectronAllocation, End, Nuc, PracticalWord


def test_methane_compact_smiles():
    assert PracticalWord([Nuc("C"), End()]).to_smiles(remove_hs=True) == "C"


def test_ethene_compact_smiles():
    p = PracticalWord([Nuc("C"), Nuc("C"), Con(0, 1, 2, 2), End()])
    assert p.to_smiles(remove_hs=True) == "C=C"


def test_acetylene_compact_smiles():
    p = PracticalWord([Nuc("C"), Nuc("C"), Con(0, 1, 3, 3), End()])
    assert p.to_smiles(remove_hs=True) == "C#C"


def test_zero_allocation_preserves_connection_as_unspecified_bond():
    d = ElectronAllocation.from_indexed(["C", "C"], [6, 6], [(0, 1, 0, 0)])
    smiles = d.to_smiles()
    assert smiles is not None
    assert "~" in smiles


def test_unassigned_electrons_are_not_silently_localized():
    d = ElectronAllocation.from_indexed(["C"], [6], [], num_unassigned_electrons=1)
    assert d.to_rdkit_mol() is None
    assert d.to_smiles() is None


def test_unassigned_count_can_be_carried_as_mol_property_when_requested():
    d = ElectronAllocation.from_indexed(["C"], [6], [], num_unassigned_electrons=1)
    mol = d.to_rdkit_mol(allow_unassigned=True, sanitize=False)
    assert mol.GetIntProp("_CCSUnassignedElectrons") == 1