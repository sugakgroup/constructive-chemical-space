from constructive_chemical_space import Arg, Con, End, Macro, MacroCon, Nuc, PracticalWord


def element_count(m, symbol):
    return sum(v.identity.symbol == symbol for v in m.nuclei)


def test_hydrogen_is_outside_automatic_completion_domain():
    m = PracticalWord([Nuc("H"), End()]).to_constitution()
    assert element_count(m, "H") == 1
    assert m.num_electrons == 1


def test_isolated_carbon_completes_to_ch4_composition():
    p = PracticalWord([Nuc("C"), End()])
    assert p.completion_counts() == (4,)
    m = p.to_constitution()
    assert element_count(m, "C") == 1
    assert element_count(m, "H") == 4


def test_ethene_completion():
    p = PracticalWord([Nuc("C"), Nuc("C"), Con(0, 1, 2, 2), End()])
    assert p.completion_counts() == (2, 2)
    m = p.to_constitution()
    assert element_count(m, "C") == 2
    assert element_count(m, "H") == 4


def test_donor_allocation_ammonia_borane_completion():
    p = PracticalWord([Nuc("N"), Nuc("B"), Con(0, 1, 2, 0), End()])
    assert p.completion_counts() == (3, 3)
    assert element_count(p.to_constitution(), "H") == 6


def test_first_end_is_terminal_and_suffix_is_ignored():
    m = PracticalWord([Nuc("C"), End(), Nuc("O"), End()]).to_constitution()
    assert element_count(m, "C") == 1
    assert element_count(m, "O") == 0
    assert element_count(m, "H") == 4


def test_ineligible_parent_receives_no_h_without_invalidating_prefix():
    p = PracticalWord([Nuc("C", 0), End()])
    info = p.completion_plan()[0]
    assert info.in_completion_domain
    assert not info.eligible
    assert info.h == 0
    assert p.replay().num_electrons == 0


def test_reservation_suppresses_completion_but_not_core_form():
    assert PracticalWord([Nuc("C", 6, 0), End()]).completion_counts() == (4,)
    assert PracticalWord([Nuc("C", 6, 4), End()]).completion_counts() == (0,)


def test_invalid_macro_argument_makes_entire_invocation_noop():
    add_c = Macro(
        "add_c",
        num_args=1,
        body=[Nuc("C"), MacroCon(Arg(0), -1, 1, 1)],
    )
    m = PracticalWord([Nuc("C"), add_c(5), End()]).to_constitution()
    assert element_count(m, "C") == 1


def test_gly_macro_matches_s4_5_one_residue_completion():
    gly = Macro(
        "Gly",
        body=[
            Nuc("N"),
            Con(-2, -1, 1, 1),
            Nuc("C"),
            Con(-2, -1, 1, 1),
            Nuc("O"),
            Nuc("C"),
            Con(-3, -1, 1, 1),
            Con(-2, -1, 2, 2),
        ],
    )
    m = PracticalWord([
        gly(),
        Nuc("O"),
        Con(-2, -1, 1, 1),
        End(),
    ]).to_constitution()
    assert {x: element_count(m, x) for x in ("C", "H", "N", "O")} == {
        "C": 2,
        "H": 5,
        "N": 1,
        "O": 2,
    }
    assert m.num_electrons == 40