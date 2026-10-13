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

def test_s5_word_has_lexicographic_connections_and_oriented_endpoints():
    from constructive_chemical_space import AllocatedConnection, Nucleus, Nuc, Con, Ele

    n, o, h = Nucleus("N"), Nucleus("O"), Nucleus("H")
    d = ElectronAllocation(
        [n, o, h], {n: 7, o: 8, h: 1},
        [AllocatedConnection(h, n, 1, 2), AllocatedConnection(n, o, 1, 3)],
        num_unassigned_electrons=2,
    )
    order = (o, h, n)
    expected = (
        Nuc("O", 8), Nuc("H", 1), Nuc("N", 7),
        Ele(), Ele(),
        Con(0, 2, 3, 1), Con(1, 2, 1, 2),
    )
    assert d.to_word(order=order).tokens == expected
    assert d.to_word(order=order).tokens == d.to_word(order=order).tokens
    assert d.to_word(order=order).replay() == d
    assert d.to_word(order=tuple(reversed(order))).replay() == d
    assert d.to_word().replay() == d


def test_word_order_requires_an_exact_enumeration():
    import pytest
    from constructive_chemical_space import Nucleus

    a, b = Nucleus("H"), Nucleus("H")
    d = ElectronAllocation(
        [a, b], {a: 1, b: 1},
        [],
    )
    with pytest.raises(ValueError):
        d.to_word(order=(a, a))
    with pytest.raises(ValueError):
        d.to_word(order=(a,))
    with pytest.raises(ValueError):
        d.to_word(order=(a, b, Nucleus("H")))


def test_s5_small_record_hydrogen_zero_connection_is_replayed_exactly():
    # S5.1 hydrogen forms A and B share a constitution but not a form.
    a = ElectronAllocation.from_indexed(
        ["H", "H"], [1, 1], [(0, 1, 1, 1)]
    )
    b = ElectronAllocation.from_indexed(
        ["H", "H"], [1, 1], [(0, 1, 0, 0)]
    )
    assert a.to_word().replay() == a
    assert b.to_word().replay() == b
    assert a != b
    assert a.to_constitution() == b.to_constitution()
