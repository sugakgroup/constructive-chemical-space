from constructive_chemical_space import Constitution, ElectronAllocation


def test_exact_constructability_witness_replays_target():
    d = ElectronAllocation.from_indexed(
        ["O", "H", "H"],
        [8, 1, 1],
        [(0, 1, 1, 1), (0, 2, 1, 1)],
    )
    assert d.to_word().replay() == d


def test_exact_constructability_with_unassigned_electrons_and_zero_allocations():
    d = ElectronAllocation.from_indexed(
        ["B", "B"],
        [5, 5],
        [(0, 1, 0, 0)],
        num_unassigned_electrons=2,
    )
    assert d.to_word().replay() == d


def test_all_unassigned_reference_lift_maps_back_to_constitution():
    m = Constitution.from_indexed(["O", "H", "H"], [(0, 1), (0, 2)], 10)
    d = m.all_unassigned_form()
    assert d.to_constitution() == m
    assert d.num_unassigned_electrons == 10
    assert all(d.assigned_electrons(v) == 0 for v in d.nuclei)
    assert all(e.allocation1 == 0 and e.allocation2 == 0 for e in d.allocated_connections)


def test_empty_carrier_reference_lift():
    m = Constitution.from_indexed([], [], 4)
    d = m.all_unassigned_form()
    assert d.to_constitution() == m
    assert d.num_unassigned_electrons == 4