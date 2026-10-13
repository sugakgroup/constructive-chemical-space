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

def test_gly_macro_multiresidue_compositions_follow_s4_5():
    gly = Macro(
        "Gly",
        body=[
            Nuc("N"), Con(-2, -1, 1, 1),
            Nuc("C"), Con(-2, -1, 1, 1),
            Nuc("O"), Nuc("C"),
            Con(-3, -1, 1, 1), Con(-2, -1, 2, 2),
        ],
    )
    for r in (2, 3, 4):
        p = PracticalWord(
            [gly() for _ in range(r)] +
            [Nuc("O"), Con(-2, -1, 1, 1), End()]
        )
        m = p.to_constitution()
        assert {x: element_count(m, x) for x in ("C", "H", "N", "O")} == {
            "C": 2 * r, "H": 3 * r + 2, "N": r, "O": r + 1,
        }
        assert m.num_electrons == 30 * r + 10


def test_completion_domain_boundaries_and_isotopes():
    # Group 12 and 18 are excluded; 13-17 are eligible; H is excluded.
    for symbol, expected in {
        "Mg": 0, "B": 3, "C": 4, "N": 3,
        "O": 2, "F": 1, "Ne": 0, "H": 0,
    }.items():
        assert PracticalWord([Nuc(symbol), End()]).completion_counts() == (expected,)
    p = PracticalWord([Nuc("13C"), End()])
    assert p.completion_counts() == (4,)
    assert element_count(p.to_constitution(), "H") == 4


def test_end_eligibility_b_u_and_e_boundaries():
    # For C, c=2. b=0 and b=8 are eligible; -1 and 9 are not.
    for n, eligible in ((1, False), (2, True), (10, True), (11, False)):
        info = PracticalWord([Nuc("C", n), End()]).completion_plan()[0]
        assert info.eligible is eligible
        if not eligible:
            assert info.h == 0

    assert PracticalWord([Nuc("C", 6, 4), End()]).completion_plan()[0].eligible
    too_large_u = PracticalWord([Nuc("C", 6, 5), End()]).completion_plan()[0]
    assert not too_large_u.eligible and too_large_u.h == 0

    negative_e = PracticalWord([
        Nuc("C", 6, 1), Nuc("H"), Con(0, 1, 4, 0), End(),
    ]).completion_plan()[0]
    assert negative_e.E == -1
    assert not negative_e.eligible and negative_e.h == 0


def test_macro_external_arguments_valid_and_invalid():
    import pytest

    attach = Macro(
        "attach", num_args=1,
        body=[Nuc("H"), MacroCon(Arg(0), -1, 1, 1)],
    )
    valid = PracticalWord([Nuc("O"), attach(0)]).to_constitution()
    assert len(valid.nuclei) == 2 and len(valid.connections) == 1

    out_of_range = PracticalWord([Nuc("O"), attach(1)]).to_constitution()
    assert len(out_of_range.nuclei) == 1
    assert len(out_of_range.connections) == 0

    with pytest.raises(ValueError):
        attach(-1)
    with pytest.raises(ValueError):
        attach("0")
    with pytest.raises(ValueError):
        attach()
    with pytest.raises(ValueError):
        attach.expand((-1,), 1)
    with pytest.raises(ValueError):
        Macro("invalid", num_args=1,
              body=[Nuc("H"), MacroCon(Arg(1), -1, 1, 1)])
    with pytest.raises(ValueError):
        MacroCon(0, -1, 1, 1)


def test_first_end_excludes_subsequent_macro_invocations_and_core_tokens():
    gly = Macro("one", [Nuc("C")])
    p = PracticalWord([Nuc("H"), End(), gly(), Nuc("O"), End()])
    assert p.expand().tokens == (Nuc("H"),)

def test_s4_table_s1_all_reservation_values():
    # Every entry of SI Table S1, including the non-monotonic group 15-17 rows.
    expected = {
        "B": (3, 2, 1, 0, 0, 0, 0, 0),
        "C": (4, 3, 2, 1, 0, 0, 0, 0),
        "N": (3, 2, 1, 0, 1, 0, 0, 0),
        "O": (2, 1, 0, 1, 0, 1, 0, 0),
        "F": (1, 0, 1, 0, 1, 0, 1, 0),
    }
    for symbol, values in expected.items():
        for u, expected_h in enumerate(values):
            actual_h = PracticalWord([Nuc(symbol, u=u), End()]).completion_counts()
            assert actual_h == (expected_h,), (symbol, u, expected_h, actual_h)


def test_s4_table_s2_all_assigned_electron_offsets():
    expected = {
        "B": (2, 3, 4),
        "C": (3, 4, 3),
        "N": (4, 3, 2),
        "O": (3, 2, 1),
        "F": (2, 1, 0),
    }
    atomic_number = {"B": 5, "C": 6, "N": 7, "O": 8, "F": 9}
    for symbol, values in expected.items():
        for offset, expected_h in zip((-1, 0, 1), values):
            p = PracticalWord([Nuc(symbol, atomic_number[symbol] + offset), End()])
            assert p.completion_counts() == (expected_h,)


def test_s4_table_s3_expanded_position_branch():
    # P bonded to two substituents (1,1) and one O (2,2):
    # b=5, B=T=4, E=1, L=max(0,1)=1, h=1.
    p = PracticalWord([
        Nuc("P"), Nuc("H"), Nuc("H"), Nuc("O"),
        Con(0, 1, 1, 1), Con(0, 2, 1, 1),
        Con(0, 3, 2, 2), End(),
    ])
    info = p.completion_plan()[0]
    assert (info.b, info.B, info.T, info.E, info.L, info.h) == (5, 4, 4, 1, 1, 1)


def test_zero_allocation_edge_occupies_one_completion_position():
    p = PracticalWord([
        Nuc("C"), Nuc("H"), Con(0, 1, 0, 0), End(),
    ])
    info = p.completion_plan()[0]
    assert (info.B, info.T, info.E, info.L, info.h) == (0, 1, 4, 3, 2)
