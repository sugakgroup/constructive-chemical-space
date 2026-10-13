#!/usr/bin/env python3
"""Generate *source-only* S5 notebooks with literal Nuc/Ele/Con construction words.

Never generates any Jupyter output or execution count. NotebookClient using
a real Jupyter kernel must subsequently execute the resulting notebooks.
"""
from __future__ import annotations

import json
import lzma
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
S5 = ROOT / "examples" / "s5"
SOURCE = S5 / "data" / "s5_records.json.xz"
MANIFEST = S5 / "data" / "manifest.json"
NOTEBOOKS = S5 / "notebooks"


def cell(kind: str, cell_id: str, source: str) -> dict:
    result = {
        "cell_type": kind, "id": cell_id, "metadata": {},
        "source": source.splitlines(keepends=True),
    }
    if kind == "code":
        # Explicitly UNEXECUTED. This script cannot invent Jupyter output.
        result.update(execution_count=None, outputs=[])
    return result


def word_code(record: dict) -> str:
    """Literal Word([...]) for every completely specified SI Form."""
    lines = [
        "from constructive_chemical_space import Word, Nuc, Ele, Con",
        "",
        "# Literal SI S5 W(D): Nuc in displayed order, n_empty times Ele,",
        "# then lexicographically ordered Con(i, j; ai, aj).",
        "# Explicit input tokens, never auto-generated from an allocation at runtime.",
        "",
    ]
    for form in record["forms"]:
        letter = form["form"]
        lines += [f"# Form {letter}", f"word_{letter} = Word(["]
        for identity, n in zip(record["nuclei"], form["assigned"]):
            lines.append(f"    Nuc({json.dumps(identity, ensure_ascii=False)}, {n}),")
        for _ in range(form["n_empty"]):
            lines.append("    Ele(),")
        for i, j, ai, aj in form["allocated_connections"]:
            lines.append(f"    Con({i}, {j}, {ai}, {aj}),")
        lines += ["])", ""]
    return "\n".join(lines).rstrip() + "\n"


def replay_code(record: dict) -> str:
    lines = [
        "# These are real Replay results from the literal Words above,",
        "# NOT the previously loaded reference input objects.",
        "",
    ]
    for form in record["forms"]:
        letter = form["form"]
        lines += [
            f"# Form {letter}: Word -> Replay -> ElectronAllocation -> Constitution",
            f"state_{letter} = word_{letter}.replay_state()",
            f"allocation_{letter} = state_{letter}.form",
            f"constitution_{letter} = allocation_{letter}.to_constitution()",
            f"reference_{letter} = example.forms[{json.dumps(letter)}]",
            f"assert word_{letter}.tokens == reference_{letter}.to_word(order=example.order).tokens",
            f"assert allocation_{letter}.to_word(order=state_{letter}.address_table).tokens == word_{letter}.tokens",
            f"assert tuple(v.identity for v in state_{letter}.address_table) == tuple(v.identity for v in example.order)",
            f"assert allocation_{letter}.num_electrons == reference_{letter}.num_electrons",
            f"assert allocation_{letter}.num_unassigned_electrons == reference_{letter}.num_unassigned_electrons",
            f"assert constitution_{letter}.num_electrons == example.constitution.num_electrons",
            f"assert constitution_{letter}.total_charge == example.constitution.total_charge",
            f"assert len(constitution_{letter}.connections) == len(example.constitution.connections)",
            f"signature_{letter} = indexed_constitution_signature(constitution_{letter}, state_{letter}.address_table)",
            f"assert signature_{letter} == reference_constitution_signature",
            f"assert indexed_allocation_signature(allocation_{letter}, state_{letter}.address_table) == indexed_allocation_signature(reference_{letter}, example.order)",
            "if len(example.order) <= 15:",
            f"    assert allocation_{letter} == reference_{letter}",
            f"    assert constitution_{letter} == example.constitution",
            f"print('Form {letter}: explicit Word -> Replay -> D -> M verified; '",
            f"      'N_e =', allocation_{letter}.num_electrons,",
            f"      ', tokens =', len(word_{letter}.tokens))",
            "",
        ]
    # Exact equality of all resulting constitutions under the displayed
    # index bijections, including the large/high-symmetry S5 examples.
    labels = [item["form"] for item in record["forms"]]
    for letter in labels[1:]:
        lines += [
            f"assert signature_{letter} == signature_{labels[0]}",
        ]
    lines.append(
        f"print('All {len(labels)} word-generated constitutions are identical "
        f"under the S5 displayed-index bijection.')"
    )
    return "\n".join(lines).rstrip() + "\n"


def constitution_code(record: dict) -> str:
    lines = [
        "# Display the constitution from EACH explicit word, including B/C/D.",
        "# Index labels are stored in metadata but omitted from clean figures.",
        "# Both drawing layers use an explicitly crossing-aware 2D layout.",
        "size = len(example.order)",
        "# Increasing coordinates alone is canceled by RDKit auto-fit.",
        "# Broaden the actual canvas proportionally to graph congestion.",
        "degrees = [0] * size",
        "for i, j in example.record['connections']:",
        "    degrees[i] += 1",
        "    degrees[j] += 1",
        "crowding = min(1.85, 1.0 + 0.11 * max(0, max(degrees, default=0) - 3))",
        "base_width = 2000 if size > 100 else (1650 if size > 55 else (1250 if size > 20 else 950))",
        "base_height = 1450 if size > 100 else (1200 if size > 55 else (900 if size > 20 else 700))",
        "draw_options = dict(show_hydrogens=True, show_indices=False,",
        "                    show_total_charge=True, allocation_label_mode='endpoint',",
        "                    layout='auto',",
        "                    width=round(base_width * crowding),",
        "                    height=round(base_height * crowding))",
    ]
    for form in record["forms"]:
        letter = form["form"]
        lines += [
            f"print('Constitution generated from word_{letter}')",
            f"constitution_{letter}.show(order=state_{letter}.address_table, **draw_options)",
        ]
    return "\n".join(lines) + "\n"


def allocation_code(record: dict) -> str:
    lines = [
        "# Each image below visualizes a newly REPLAYED allocation object.",
        "# Reference fixture forms are NOT used to generate these pictures.",
    ]
    for form in record["forms"]:
        letter = form["form"]
        lines += [
            f"print('Electron allocation: Form {letter} (generated by word_{letter})')",
            f"allocation_{letter}.show(order=state_{letter}.address_table, **draw_options)",
        ]
    return "\n".join(lines) + "\n"


def export_code(record: dict, notebook_name: str) -> str:
    stem = Path(notebook_name).stem
    lines = [
        "# Optional SVG export of every Form's replay-generated D AND M.",
        "EXPORT = False",
        "if EXPORT:",
        "    output_dir = repo_root / 'examples' / 's5' / 'outputs'",
        "    output_dir.mkdir(parents=True, exist_ok=True)",
    ]
    for form in record["forms"]:
        letter = form["form"]
        lines += [
            f"    constitution_{letter}.save(output_dir / '{stem}_Constitution_{letter}.svg',",
            f"                               order=state_{letter}.address_table, **draw_options)",
            f"    allocation_{letter}.save(output_dir / '{stem}_Form_{letter}.svg',",
            f"                           order=state_{letter}.address_table, **draw_options)",
        ]
    lines.append("    print('Saved actual images under:', output_dir)")
    return "\n".join(lines) + "\n"


def create_notebook(index: dict, record: dict) -> dict:
    n = record["number"]
    prefix = f"s5-{n:02d}-"
    exp = record["expected"]
    title = (
        f"# {record['section']} · Record {n:02d} — {record['name']}\n\n"
        f"**SI:** page(s) {', '.join(map(str, record['source']['pdf_pages']))}; "
        f"Table S4 record {n}.  **Forms:** "
        f"{', '.join(f['form'] for f in record['forms'])}.\n\n"
        f"**Constitution:** |V| = {exp['nuclei']}, "
        f"|S| = {exp['connections']}, "
        f"N_e = {exp['electrons']}, Q = {exp['charge']:+d}.\n\n"
        "This notebook explicitly defines each construction word, executes "
        "its Replay, constructs the resulting electron allocation and "
        "constitution, checks both against the source, and draws the results. "
        "Nuclear indices are temporary addresses only. All H atoms and "
        "connections are shown by default.\n"
    )
    setup = "\n".join([
        "from pathlib import Path",
        "import sys",
        "",
        "cwd = Path.cwd().resolve()",
        "repo_root = next((p for p in (cwd, *cwd.parents)",
        "                  if (p / 'examples' / 's5' / '_common.py').is_file()), None)",
        "if repo_root is None:",
        "    raise FileNotFoundError('Open within the CCS repository tree')",
        "if str(repo_root) not in sys.path:",
        "    sys.path.insert(0, str(repo_root))",
        "from examples.s5._common import (",
        "    build_example, summarize_example, verify_example,",
        "    indexed_constitution_signature, indexed_allocation_signature,",
        ")",
        "",
        f"example = build_example({n})",
        "print(summarize_example(example))",
        "replay_lengths = verify_example(example)",
        "reference_constitution_signature = indexed_constitution_signature(example.constitution, example.order)",
        "print('SI source record verified:', replay_lengths)",
        "",
    ])
    return {
        "cells": [
            cell("markdown", prefix + "intro", title),
            cell("code", prefix + "load", setup),
            cell("markdown", prefix + "words-title",
                 "## 1. Explicit construction words\n\n"
                 "Each Form has a concrete Word([Nuc(...), Ele(), Con(...)]) "
                 "declaration with every individual token written out. "
                 "These words do not come from a runtime call to to_word().\n"),
            cell("code", prefix + "words", word_code(record)),
            cell("markdown", prefix + "replay-title",
                 "## 2. Replay Word, construct D and M, verify the outputs\n\n"
                 "Each state_X is a real Replay state; allocation_X is the "
                 "resulting intrinsic electron-allocation form, and "
                 "constitution_X its CCS constitutional quotient. "
                 "The source records are used only as verification references.\n"),
            cell("code", prefix + "replay", replay_code(record)),
            cell("markdown", prefix + "view-title",
                 "## 3. Visualize the generated CCS objects\n\n"
                 "Each Form A/B/C/... displays its own Replay-generated "
                 "Constitution and its own ElectronAllocation. The complete "
                 "constitutional signatures are proven equal under the S5 "
                 "index bijection. Nucleus labels omit temporary #indices "
                 "by default, but retain them in the SVG metadata. Residual "
                 "valence electrons (groups 13-17) appear as dots, and CCS "
                 "assignment imbalance as signs near the nucleus. Planar "
                 "graphs use crossing-free 2D connections; genuinely nonplanar "
                 "graphs are explicitly marked. All H remain explicit.\n"),
            cell("code", prefix + "constitution", constitution_code(record)),
            cell("code", prefix + "allocations", allocation_code(record)),
            cell("markdown", prefix + "export-title",
                 "## 4. Optional export of replayed results\n\n"
                 "EXPORT=False by default; change to True to save SVGs under "
                 "the ignored examples/s5/outputs directory.\n"),
            cell("code", prefix + "export", export_code(record, index["notebook"])),
        ],
        "metadata": {
            "kernelspec": {"display_name": "Python 3",
                           "language": "python", "name": "python3"},
            "language_info": {"name": "python"},
            "ccs_s5_record": n,
            "ccs_s5_name": record["name"],
            "ccs_s5_section": record["section"],
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def main() -> int:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    with lzma.open(SOURCE, "rt", encoding="utf-8") as f:
        data = json.load(f)
    records = data["records"]
    if len(records) != 55 or sum(len(r["forms"]) for r in records) != 84:
        raise RuntimeError("Expected exactly 55 records and 84 forms")
    for index, record in zip(manifest["records"], records, strict=True):
        if index["number"] != record["number"]:
            raise ValueError("Manifest and dataset are inconsistent")
        path = NOTEBOOKS / index["notebook"]
        if not path.is_file():
            raise FileNotFoundError(path)
        # Source-only notebook with null execution_count and empty outputs.
        # The real-kernel execution script runs in a separate workflow step.
        path.write_text(
            json.dumps(create_notebook(index, record), ensure_ascii=False,
                       indent=1) + "\n", encoding="utf-8"
        )
        print(
            f"SOURCE ONLY {record['number']:02d}/55 {path.name}: "
            f"{len(record['forms'])} literal Word definitions",
            flush=True,
        )
    print("Generated 55 UNEXECUTED notebooks. Real Jupyter execution required.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
