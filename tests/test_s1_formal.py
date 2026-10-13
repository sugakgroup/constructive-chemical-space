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
        num_unassigned_electrons=3,
    )
    assert d.num_electrons == 5
    assert d.num_unassigned_electrons == 3
    assert d.to_constitution().num_electrons == 5

def test_intrinsic_form_equality_under_arbitrary_relabeling():
    # Same graph, isotope labels and endpoint allocation after a permutation.
    a = ElectronAllocation.from_indexed(
        ["13C", "N", "H", "H"],
        [6, 7, 1, 1],
        [(0, 1, 2, 0), (1, 2, 1, 1), (0, 3, 1, 1)],
        num_unassigned_electrons=2,
    )
    b = ElectronAllocation.from_indexed(
        ["H", "H", "N", "13C"],
        [1, 1, 7, 6],
        [(3, 2, 2, 0), (0, 2, 1, 1), (1, 3, 1, 1)],
        num_unassigned_electrons=2,
    )
    assert a == b
    assert a.to_constitution() == b.to_constitution()


def test_symmetric_form_isomorphism_preserves_oriented_endpoint_allocations():
    a = ElectronAllocation.from_indexed(
        ["C", "C", "C"], [6, 6, 6],
        [(0, 1, 1, 2), (1, 2, 2, 1)],
    )
    b = ElectronAllocation.from_indexed(
        ["C", "C", "C"], [6, 6, 6],
        [(2, 1, 1, 2), (1, 0, 2, 1)],
    )
    nonisomorphic = ElectronAllocation.from_indexed(
        ["C", "C", "C"], [6, 6, 6],
        [(0, 1, 1, 2), (1, 2, 1, 2)],
    )
    assert a == b
    assert a != nonisomorphic
    assert a.to_constitution() == nonisomorphic.to_constitution()


def test_assignment_imbalance_is_not_a_lewis_charge_api():
    d = ElectronAllocation.from_indexed(
        ["N", "B"], [7, 5], [(0, 1, 2, 0)]
    )
    n = next(v for v in d.nuclei if v.identity.symbol == "N")
    assert d.assignment_imbalance(n) == 0
    with pytest.warns(DeprecationWarning, match="assignment_imbalance"):
        assert d.formal_charge(n) == 0


def test_constitution_immutable_value_object():
    from dataclasses import FrozenInstanceError

    m = Constitution.from_indexed(["H", "H"], [(0, 1)], 2)
    with pytest.raises(FrozenInstanceError):
        m.num_electrons = -1
    with pytest.raises(FrozenInstanceError):
        m.connections = frozenset()
    with pytest.raises(FrozenInstanceError):
        m.collection = None
    assert m.num_electrons == 2
    assert len(m.connections) == 1


def test_form_immutable_with_defensively_copied_assignments():
    from dataclasses import FrozenInstanceError

    h1, h2 = Nucleus("H"), Nucleus("H")
    assignments = {h1: 1, h2: 1}
    d = ElectronAllocation(
        [h1, h2], assignments, [AllocatedConnection(h1, h2, 1, 1)]
    )
    assignments[h1] = -7  # Caller-owned mutable input cannot affect the form.
    assert d.assigned_electrons(h1) == 1
    with pytest.raises(TypeError):
        d._assigned[h1] = -1
    copy = d.assigned_mapping()
    copy[h1] = -1
    assert d.assigned_electrons(h1) == 1
    with pytest.raises(FrozenInstanceError):
        d.num_unassigned_electrons = -1
    with pytest.raises(FrozenInstanceError):
        d.allocated_connections = frozenset()
    with pytest.raises(FrozenInstanceError):
        d._assigned = {}
    assert d.to_word(order=(h1, h2)).replay() == d

def test_intrinsic_equality_exhaustive_four_nucleus_relabelings():
    from itertools import permutations

    identities = ["C", "N", "H", "H"]
    assignments = [6, 7, 1, 1]
    edges = [(0, 1, 2, 0), (1, 2, 1, 1), (0, 3, 1, 1)]
    original = ElectronAllocation.from_indexed(
        identities, assignments, edges, num_unassigned_electrons=1
    )
    for permutation in permutations(range(4)):
        inverse = {old: new for new, old in enumerate(permutation)}
        relabeled = ElectronAllocation.from_indexed(
            [identities[old] for old in permutation],
            [assignments[old] for old in permutation],
            [(inverse[i], inverse[j], ai, aj) for i, j, ai, aj in edges],
            num_unassigned_electrons=1,
        )
        assert relabeled == original
        assert relabeled.to_constitution() == original.to_constitution()


def test_total_charge_includes_unassigned_electron_contribution():
    d = ElectronAllocation.from_indexed(
        ["O", "H"], [8, 1], [(0, 1, 1, 1)],
        num_unassigned_electrons=2,
    )
    assert sum(d.assignment_imbalance(v) for v in d.nuclei) == 0
    assert d.to_constitution().total_charge == -2
