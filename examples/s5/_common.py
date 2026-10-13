"""Example-only helpers for SI S5 notebooks. No CCS semantics are defined here."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping, Sequence
import hashlib
import json
import lzma

from constructive_chemical_space import (
    AllocatedConnection,
    Con,
    Connection,
    Constitution,
    Ele,
    ElectronAllocation,
    Nuc,
    Nucleus,
)

S5_ROOT = Path(__file__).resolve().parent
MANIFEST_PATH = S5_ROOT / "data" / "manifest.json"
OUTPUT_DIR = S5_ROOT / "outputs"


@lru_cache(maxsize=1)
def _manifest() -> dict[str, Any]:
    with MANIFEST_PATH.open("r", encoding="utf-8") as stream:
        result = json.load(stream)
    if result.get("schema") != "ccs-s5-notebooks-1.0":
        raise ValueError("Unsupported S5 notebook manifest schema")
    return result


def load_manifest() -> dict[str, Any]:
    """Read the displayed S5 record index and source metadata."""
    return deepcopy(_manifest())


@lru_cache(maxsize=1)
def _records() -> dict[str, Any]:
    index = _manifest()
    location = S5_ROOT / "data" / index["data_file"]
    digest = hashlib.sha256(location.read_bytes()).hexdigest()
    if digest != index["data_sha256"]:
        raise ValueError("The S5 data file checksum does not match the manifest")
    with lzma.open(location, "rt", encoding="utf-8") as stream:
        result = json.load(stream)
    if result.get("schema") != "ccs-s5-1.0":
        raise ValueError("Unsupported S5 record data schema")
    return result


def load_record(number: int) -> dict[str, Any]:
    """Return the full S5 record, including every Form A/B/C/D."""
    if type(number) is not int or not 1 <= number <= 55:
        raise ValueError("S5 record number must be an integer from 1 through 55")
    record = _records()["records"][number - 1]
    if record["number"] != number:
        raise ValueError("S5 fixture order is inconsistent")
    return deepcopy(record)


@dataclass(frozen=True)
class S5Example:
    """A record with a shared, temporary nucleus order and all its CCS Forms."""

    number: int
    name: str
    section: str
    order: tuple[Nucleus, ...]
    constitution: Constitution
    forms: Mapping[str, ElectronAllocation]
    record: Mapping[str, Any]


def build_example(number: int) -> S5Example:
    """Instantiate one S5 record using the existing CCS v1.0 classes."""
    record = load_record(number)
    index = _manifest()["records"][number - 1]
    if (index["number"] != number or index["name"] != record["name"]
            or index["expected"] != record["expected"]):
        raise ValueError("S5 manifest and full record are inconsistent")

    nuclei = tuple(Nucleus(identity) for identity in record["nuclei"])
    constitution = Constitution(
        nuclei,
        [Connection(nuclei[i], nuclei[j]) for i, j in record["connections"]],
        record["expected"]["electrons"],
    )
    forms = {}
    for form in record["forms"]:
        forms[form["form"]] = ElectronAllocation(
            nuclei,
            {nucleus: count for nucleus, count in zip(nuclei, form["assigned"])},
            [AllocatedConnection(nuclei[i], nuclei[j], a, b)
             for i, j, a, b in form["allocated_connections"]],
            form["n_empty"],
        )
    return S5Example(
        number=number, name=record["name"], section=record["section"],
        order=nuclei, constitution=constitution,
        forms=MappingProxyType(forms), record=record,
    )


def summarize_example(example: S5Example) -> str:
    expected = example.record["expected"]
    letters = ", ".join(example.forms)
    return (
        f"{example.section}, Record {example.number:02d}: {example.name}\n"
        f"|V|={expected['nuclei']}, |S|={expected['connections']}, "
        f"N_e={expected['electrons']}, Q={expected['charge']:+d}; "
        f"Forms: {letters}"
    )


def verify_example(example: S5Example) -> dict[str, int]:
    """Assert Table S4 and complete W(D) Replay with an explicit nucleus bijection.

    Using the replay address table avoids factorial generic isomorphism searches
    on large symmetric S5 structures without changing intrinsic identity.
    """
    record = example.record
    expected = record["expected"]
    order = example.order
    assert len(order) == expected["nuclei"]
    assert len(example.constitution.connections) == expected["connections"]
    assert example.constitution.num_electrons == expected["electrons"]
    assert example.constitution.total_charge == expected["charge"]
    assert len(example.forms) == expected["forms"]
    expected_edges = {tuple(x) for x in record["connections"]}
    result = {}

    for raw in record["forms"]:
        letter = raw["form"]
        form = example.forms[letter]
        assert form.num_electrons == expected["electrons"]
        assert sum(form.assignment_imbalance(v) for v in order) - raw["n_empty"] == expected["charge"]
        quotient = form.to_constitution()
        assert quotient.nuclei == example.constitution.nuclei
        assert quotient.connections == example.constitution.connections
        assert quotient.num_electrons == example.constitution.num_electrons
        assert {tuple(item[:2]) for item in raw["allocated_connections"]} == expected_edges

        word = form.to_word(order=order)
        assert word.tokens == form.to_word(order=order).tokens
        assert word.tokens[:len(order)] == tuple(
            Nuc(v.identity, count) for v, count in zip(order, raw["assigned"])
        )
        e_end = len(order) + raw["n_empty"]
        assert word.tokens[len(order):e_end] == (Ele(),) * raw["n_empty"]
        assert word.tokens[e_end:] == tuple(
            Con(i, j, a, b) for i, j, a, b in raw["allocated_connections"]
        )
        assert len(word.tokens) == len(order) + raw["n_empty"] + len(expected_edges)

        state = word.replay_state()
        assert len(state.address_table) == len(order)
        assert all(state.address_table[i].identity == v.identity for i, v in enumerate(order))
        assert state.form.num_unassigned_electrons == raw["n_empty"]
        assert len(state.form.connections) == len(expected_edges)
        assert state.form.num_electrons == expected["electrons"]
        for i, count in enumerate(raw["assigned"]):
            assert state.form.assigned_electrons(state.address_table[i]) == count
        for i, j, a, b in raw["allocated_connections"]:
            edge = state.form.allocation_on(state.address_table[i], state.address_table[j])
            assert edge is not None
            assert edge.allocation_for(state.address_table[i]) == a
            assert edge.allocation_for(state.address_table[j]) == b
        # The full matched address correspondence is the explicit S1 bijection.
        if len(order) <= 15:
            assert state.form == form
        result[letter] = len(word.tokens)
    return result



def indexed_constitution_signature(
    constitution: Constitution, order: Sequence[Nucleus]
):
    """Exact S1 constitution under a specified, explicitly witnessed bijection.

    The ordering is external to molecular identity. Comparing these complete
    signatures proves identity under this supplied address correspondence,
    even for highly symmetric large molecular graphs.
    """
    address = {nucleus: i for i, nucleus in enumerate(order)}
    if (len(order) != len(constitution.nuclei)
            or len(address) != len(constitution.nuclei)
            or set(address) != set(constitution.nuclei)):
        raise ValueError("order must enumerate the constitution exactly")
    edges = tuple(sorted(
        tuple(sorted((address[e.endpoint1], address[e.endpoint2])))
        for e in constitution.connections
    ))
    return (
        tuple(nucleus.identity for nucleus in order),
        edges,
        constitution.num_electrons,
    )


def indexed_allocation_signature(
    form: ElectronAllocation, order: Sequence[Nucleus]
):
    """Complete CCS D fields, oriented relative to external index order."""
    address = {nucleus: i for i, nucleus in enumerate(order)}
    if (len(order) != len(form.nuclei)
            or len(address) != len(form.nuclei)
            or set(address) != set(form.nuclei)):
        raise ValueError("order must enumerate the electron allocation exactly")
    endpoints = []
    for edge in form.allocated_connections:
        i, j = sorted((address[edge.endpoint1], address[edge.endpoint2]))
        endpoints.append((
            i, j,
            edge.allocation_for(order[i]),
            edge.allocation_for(order[j]),
        ))
    return (
        tuple(nucleus.identity for nucleus in order),
        tuple(form.assigned_electrons(v) for v in order),
        tuple(sorted(endpoints)),
        form.num_unassigned_electrons,
    )



def export_example(
    example: S5Example,
    *,
    output_dir: str | Path | None = None,
    formats: Sequence[str] = ("svg",),
    include_constitution: bool = True,
    **draw_kwargs: Any,
) -> list[Path]:
    """Save constitution and every allocation in S5/outputs/ when requested."""
    out = OUTPUT_DIR if output_dir is None else Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    filename = _manifest()["records"][example.number - 1]["notebook"]
    stem = Path(filename).stem
    options = {"order": example.order, "show_indices": True,
               "show_total_charge": True, **draw_kwargs}
    saved: list[Path] = []
    for fmt in formats:
        if fmt not in {"svg", "png"}:
            raise ValueError("Supported S5 figure formats are svg and png")
        if include_constitution:
            target = out / f"{stem}_Constitution.{fmt}"
            example.constitution.save(target, **options)
            saved.append(target)
        for letter, form in example.forms.items():
            target = out / f"{stem}_Form_{letter}.{fmt}"
            form.save(target, **options)
            saved.append(target)
    return saved
