from constructive_chemical_space import Con, Ele, Nuc, Word


def test_nuc_default_assigns_atomic_number():
    state = Word([Nuc("O")]).replay_state()
    o = state.address_table[0]
    assert state.form.assigned_electrons(o) == 8


def test_explicit_zero_assignment_is_not_neutral_shorthand():
    state = Word([Nuc("H", 0)]).replay_state()
    h = state.address_table[0]
    assert state.form.assigned_electrons(h) == 0


def test_reservation_u_does_not_change_core_replay():
    assert Word([Nuc("C", 6, 0)]).replay() == Word([Nuc("C", 6, 4)]).replay()


def test_unassigned_electron_installation():
    d = Word([Ele()]).replay()
    assert len(d.nuclei) == 0
    assert d.num_unassigned_electrons == 1
    assert d.num_electrons == 1


def test_interleaved_water_construction():
    w = Word([
        Nuc("O"),
        Nuc("H"),
        Con(0, 1, 1, 1),
        Nuc("H"),
        Con(0, 2, 1, 1),
    ])
    state = w.replay_state()
    o, h1, h2 = state.address_table
    assert state.form.residual_electrons(o) == 6
    assert state.form.residual_electrons(h1) == 0
    assert state.form.residual_electrons(h2) == 0


def test_noop_and_relative_address_behavior():
    tested = Word([
        Con(0, 1, 1, 1),
        Nuc("H"),
        Nuc("H"),
        Con(-2, -1, 1, 1),
        Con(0, 1, 0, 0),
    ]).replay()
    reference = Word([Nuc("H"), Nuc("H"), Con(0, 1, 1, 1)]).replay()
    assert tested == reference


def test_two_addresses_resolving_to_same_nucleus_make_con_noop():
    tested = Word([
        Nuc("H"),
        Nuc("H"),
        Con(1, -1, 0, 0),
    ]).replay()
    assert len(tested.connections) == 0


def test_insufficient_residual_budget_makes_con_noop():
    tested = Word([
        Nuc("H"),
        Nuc("H"),
        Con(0, 1, 2, 0),
    ]).replay()
    assert len(tested.connections) == 0

def test_con_noop_returns_entire_original_execution_state():
    base = Word([Nuc("H"), Nuc("H")]).replay_state()
    assert base.apply(Con(0, 1, 2, 0)) is base
    assert base.apply(Con(1, -1, 0, 0)) is base
    assert base.apply(Con(0, 9, 0, 0)) is base
    with_edge = base.apply(Con(0, 1, 0, 0))
    assert with_edge is not base
    assert with_edge.apply(Con(0, 1, 1, 1)) is with_edge


def test_nuc_and_ele_remain_executable_after_inapplicable_connection():
    d = Word([Con(0, 1, 1, 1), Nuc("H"), Ele(), Nuc("H")]).replay()
    assert len(d.nuclei) == 2
    assert d.num_unassigned_electrons == 1
    assert len(d.connections) == 0
