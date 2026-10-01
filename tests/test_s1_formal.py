import pytest

from constructive_chemical_space import (
    AllocatedConnection,
    Connection,
    Constitution,
    ElectronAllocation,
    Nucleus,
)


def test_nuclei_are_abstract_carrier_elements():
    h1 = Nucleus("H")
    h2 = Nucleus("H")
    assert h1 is not h2
    assert h1 != h2
    assert len({h1, h2}) == 2


def test_constitution_equality_ignores_temporary_enumeration():
    a = Constitution.from_indexed(["O", "H", "H"], [(0, 1), (0, 2)], 10)
    b = Constitution.from_indexed(["H", "O", "H"], [(1, 0), (1, 2)], 10)
    assert a == b


def test_constitution_rejects_negative_convenience_indices():
    with pytest.raises(ValueError):
        Constitution.from_indexed(["H", "H"], [(-1, 0)], 2)


def test_connection_is_simple_and_undirected():
    h1, h2 = Nucleus("H"), Nucleus("H")
    assert Connection(h1, h2) == Connection(h2, h1)
    with pytest.raises(ValueError):
        Connection(h1, h1)


def test_well_formedness_rejects_overallocation():
    h1, h2 = Nucleus("H"), Nucleus("H")
    with pytest.raises(ValueError):
        ElectronAllocation(
            [h1, h2],
            {h1: 1, h2: 1},
            [AllocatedConnection(h1, h2, 2, 0)],
            0,
        )


def test_zero_allocation_connection_is_retained():
    d = ElectronAllocation.from_indexed(
        ["C", "C"],
        [6, 6],
        [(0, 1, 0, 0)],
    )
    assert len(d.connections) == 1
    assert len(d.to_constitution().connections) == 1


def test_electron_allocation_quotient_retains_total_electron_count():
    d = ElectronAllocation.from_indexed(
        ["H", "H"],
        [1, 1],
        [(0, 1, 1, 1)],
        num_unassigned_electron=3,
    )
    assert d.num_electrons == 5
    assert d.num_unassigned_electrons == 3
    assert d.to_constitution().num_electrons == 5