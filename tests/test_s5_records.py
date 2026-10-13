"""Machine-verifiable SI S5 example records (55 records, 84 forms)."""
import json
import lzma
from pathlib import Path

import pytest

from constructive_chemical_space import (
    AllocatedConnection, Con, Connection, Constitution, Ele,
    ElectronAllocation, Nuc, Nucleus,
)

FILE = Path(__file__).resolve().parents[1] / "data" / "s5_records.json.xz"
with lzma.open(FILE, "rt", encoding="utf-8") as stream:
    DATA = json.load(stream)
RECORDS = DATA["records"]
CASES = [(r["number"], f["form"]) for r in RECORDS for f in r["forms"]]


def prepare(record):
    nuclei = tuple(Nucleus(x) for x in record["nuclei"])
    connections = tuple(tuple(edge) for edge in record["connections"])
    expected = record["expected"]
    m = Constitution(
        nuclei,
        [Connection(nuclei[i], nuclei[j]) for i, j in connections],
        expected["electrons"],
    )
    return nuclei, connections, m


def allocation_form(entry, nuclei):
    return ElectronAllocation(
        nuclei,
        {v: n for v, n in zip(nuclei, entry["assigned"])},
        [AllocatedConnection(nuclei[i], nuclei[j], a, b)
         for i, j, a, b in entry["allocated_connections"]],
        entry["n_empty"],
    )


def test_corpus_exact_size_and_traceability():
    assert DATA["schema"] == "ccs-s5-1.0"
    assert len(RECORDS) == 55
    assert len(CASES) == 84
    assert [r["number"] for r in RECORDS] == list(range(1, 56))
    assert all(len(DATA["source"][name]) == 64 for name in
               ("pdf_sha256", "docx_sha256"))
    assert sum(r["expected"]["forms"] for r in RECORDS) == 84


@pytest.mark.parametrize("number", range(1, 56))
def test_table_s4_and_complete_connection_sets(number):
    r = RECORDS[number - 1]
    nuclei, connections, m = prepare(r)
    exp = r["expected"]
    assert len(nuclei) == exp["nuclei"]
    assert len(connections) == exp["connections"]
    assert len(set(connections)) == len(connections)
    assert all(0 <= i < j < len(nuclei) for i, j in connections)
    assert m.num_electrons == exp["electrons"]
    assert m.total_charge == exp["charge"]
    assert len(r["forms"]) == exp["forms"]
    assert [x["form"] for x in r["forms"]] == [
        chr(ord("A") + i) for i in range(exp["forms"])
    ]
    assert r["source"]["pdf_pages"]
    assert r["source"]["docx_paragraph"] >= 233


@pytest.mark.parametrize("number,letter", CASES,
                         ids=[f"S5-{n:02d}-{letter}" for n, letter in CASES])
def test_every_s5_form_has_exact_deterministic_replay(number, letter):
    r = RECORDS[number - 1]
    f = next(item for item in r["forms"] if item["form"] == letter)
    nuclei, connections, m = prepare(r)
    d = allocation_form(f, nuclei)

    assert d.num_electrons == m.num_electrons
    assert d.to_constitution() == m
    assert len(d.allocated_connections) == len(connections)
    assert sum(d.assignment_imbalance(v) for v in nuclei) - f["n_empty"] == m.total_charge

    # S5 W(D) uses exactly the listed temporary nucleus order.
    word = d.to_word(order=nuclei)
    assert word.tokens == d.to_word(order=nuclei).tokens
    assert len(word.tokens) == len(nuclei) + f["n_empty"] + len(connections)
    assert word.tokens[:len(nuclei)] == tuple(
        Nuc(v.identity, n) for v, n in zip(nuclei, f["assigned"])
    )
    assert word.tokens[len(nuclei):len(nuclei) + f["n_empty"]] == (
        Ele(),
    ) * f["n_empty"]
    assert word.tokens[len(nuclei) + f["n_empty"]:] == tuple(
        Con(i, j, ai, aj) for i, j, ai, aj in f["allocated_connections"]
    )
    assert list(connections) == sorted(connections)
    assert len({tuple(a[:2]) for a in f["allocated_connections"]}) == len(connections)

    state = word.replay_state()
    assert tuple(v.identity for v in state.address_table) == tuple(v.identity for v in nuclei)
    assert state.form == d  # Intrinsic equality, not address equality.
    assert state.form.to_constitution() == m


@pytest.mark.parametrize("number", [1, 2, 9, 10, 14, 17, 22, 28, 39, 42, 46])
def test_relabeling_preserves_intrinsic_identity(number):
    r = RECORDS[number - 1]
    nuclei, _, m = prepare(r)
    perm = tuple(reversed(range(len(nuclei))))
    new_address = {old: new for new, old in enumerate(perm)}
    for f in r["forms"]:
        original = allocation_form(f, nuclei)
        relabeled = ElectronAllocation.from_indexed(
            [r["nuclei"][i] for i in perm],
            [f["assigned"][i] for i in perm],
            [(new_address[i], new_address[j], ai, aj)
             for i, j, ai, aj in f["allocated_connections"]],
            num_unassigned_electrons=f["n_empty"],
        )
        assert relabeled == original
        assert relabeled.to_constitution() == m


def test_unassigned_electron_and_zero_allocation_boundaries():
    assert [(r["number"], f["form"]) for r in RECORDS for f in r["forms"]
            if f["n_empty"] > 0] == [(43, "A"), (46, "A"), (48, "A")]
    assert RECORDS[45]["nuclei"] == []
    assert RECORDS[45]["connections"] == []
    assert RECORDS[45]["forms"][0]["n_empty"] == 1
    assert RECORDS[44]["connections"] == []  # disconnected Na/Cl
    assert RECORDS[0]["forms"][1]["allocated_connections"] == [[0, 1, 0, 0]]


def test_cross_record_identity_boundaries():
    examples = {r["number"]: r for r in RECORDS}
    assert prepare(examples[49])[2] != prepare(examples[50])[2]
    assert prepare(examples[51])[2] != prepare(examples[52])[2]
    assert prepare(examples[47])[2] != prepare(examples[48])[2]
    assert prepare(examples[39])[2] != prepare(examples[40])[2]
    # Spatial catenation and stereo labels add no S1/S5 identity field.
    assert all(examples[n]["expected"]["forms"] == 1 for n in (53, 54, 55))
