#!/usr/bin/env python3
"""Execute the 55 SI S5 notebooks in *real* Jupyter Python kernels.

Only NotebookClient/Jupyter kernel output is recorded. This script never writes
or reconstructs display_data, execute_result, stream, or execution_count values
by hand. Errors stop the run and prevent a successful completion report.

Run from any directory:
    python examples/s5/execute_notebooks.py

Dependencies: pip install -e ".[dev]" && pip install nbclient nbformat ipykernel
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
import xml.etree.ElementTree as ET

import nbclient
from nbclient import NotebookClient
import nbformat
from rdkit import rdBase


ROOT = Path(__file__).resolve().parents[2]
S5 = ROOT / "examples" / "s5"
MANIFEST = S5 / "data" / "manifest.json"
NOTEBOOKS = S5 / "notebooks"
REPORT = S5 / "execution_report.json"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git_commit() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None


def code_digest(notebook) -> str:
    cells = "\n\n".join(
        "".join(cell.source) for cell in notebook.cells
        if cell.cell_type == "code"
    )
    return sha256(cells.encode("utf-8"))


def svg_from_kernel_output(output) -> str | None:
    if output.output_type not in ("display_data", "execute_result"):
        return None
    svg = output.get("data", {}).get("image/svg+xml")
    if svg is None:
        return None
    return "".join(svg) if isinstance(svg, list) else str(svg)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--records", nargs="+", type=int, default=None,
        help="Representative record numbers; executes real kernels but does "
             "not overwrite the authoritative 55-notebook execution report.",
    )
    parser.add_argument("--report-path", type=Path, default=None)
    args = parser.parse_args(argv)

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if manifest["schema"] != "ccs-s5-notebooks-1.0":
        raise RuntimeError("Unexpected S5 manifest schema")
    all_entries = manifest["records"]
    if len(all_entries) != 55 or sum(len(x["forms"]) for x in all_entries) != 84:
        raise RuntimeError("S5 manifest must contain exactly 55 records / 84 forms")
    full_corpus = args.records is None
    if full_corpus:
        entries = all_entries
    else:
        chosen = args.records
        if (not chosen or len(chosen) != len(set(chosen))
                or any(not 1 <= num <= 55 for num in chosen)):
            raise ValueError("Representative --records must be distinct numbers 1..55")
        selected = set(chosen)
        entries = [row for row in all_entries if row["number"] in selected]
        if len(entries) != len(chosen):
            raise RuntimeError("Requested example record not found in manifest")

    records_file = S5 / "data" / manifest["data_file"]
    if sha256(records_file.read_bytes()) != manifest["data_sha256"]:
        raise RuntimeError("S5 source data checksum does not match the manifest")

    started_at = datetime.now(timezone.utc)
    input_commit = git_commit()
    results = []
    print(
        f"ACTUAL JUPYTER EXECUTION: {len(entries)} notebooks, "
        f"{sum(len(x['forms']) for x in entries)} forms; "
        f"Python={platform.python_version()}, RDKit={rdBase.rdkitVersion}, "
        f"nbclient={nbclient.__version__}; input_commit={input_commit}",
        flush=True,
    )

    for entry in entries:
        path = NOTEBOOKS / entry["notebook"]
        nb = nbformat.read(path, as_version=4)
        if nb.metadata.get("ccs_s5_record") != entry["number"]:
            raise RuntimeError(f"Wrong S5 record metadata in {path}")
        if len([c for c in nb.cells if c.cell_type == "code"]) != 6:
            raise RuntimeError(f"Unexpected code-cell count in {path}")

        # Verify that full, literal Word declarations exist in notebook source;
        # an indirect form.to_word() call is not an explicit Word definition.
        code_sources = ["".join(c.source) for c in nb.cells
                        if c.cell_type == "code"]
        for label in entry["forms"]:
            if f"word_{label} = Word([" not in code_sources[1]:
                raise RuntimeError(
                    f"Missing literal construction word for Form {label}: {path}"
                )
        if ".replay_state()" not in code_sources[2] or ".to_constitution()" not in code_sources[2]:
            raise RuntimeError(f"Missing explicit Replay-to-Constitution step in {path}")

        source_digest = code_digest(nb)
        # NotebookClient will clear old outputs and populate the real outputs
        # and execution counts from an ipykernel. We do not construct them.
        started = time.monotonic()
        client = NotebookClient(
            nb,
            timeout=300,
            kernel_name="python3",
            allow_errors=False,
            record_timing=True,
            resources={"metadata": {"path": str(ROOT)}},
        )
        executed = client.execute()
        elapsed = time.monotonic() - started
        nbformat.validate(executed)

        code_cells = [c for c in executed.cells if c.cell_type == "code"]
        if any(not isinstance(c.execution_count, int) or c.execution_count < 1
               for c in code_cells):
            raise RuntimeError(f"Kernel did not execute every code cell: {path}")

        svg_payloads = [
            svg_from_kernel_output(output)
            for cell in code_cells
            for output in cell.get("outputs", [])
        ]
        svg_payloads = [payload for payload in svg_payloads if payload is not None]
        for svg in svg_payloads:
            ET.fromstring(svg)  # Validate each *actual kernel-produced* SVG.
        # A/C/D etc EACH displays its own replay-generated M AND D.
        expected_svg_count = 2 * len(entry["forms"])
        def real_svgs(cell):
            return sum(svg_from_kernel_output(out) is not None
                       for out in cell.get("outputs", []))
        if (real_svgs(code_cells[3]) != len(entry["forms"])
                or real_svgs(code_cells[4]) != len(entry["forms"])):
            raise RuntimeError(
                f"{path}: both constitution and allocation cells must have "
                "one genuine SVG per Form"
            )
        if len(svg_payloads) != expected_svg_count:
            raise RuntimeError(
                f"{path}: expected {expected_svg_count} real SVG outputs, "
                f"received {len(svg_payloads)}"
            )

        output_hash = sha256("\n".join(svg_payloads).encode("utf-8"))
        finished_at = datetime.now(timezone.utc)
        executed.metadata["ccs_execution"] = {
            "status": "completed",
            "method": "nbclient real ipykernel execution",
            "started_utc": started_at.isoformat() if not results else results[-1]["finished_utc"],
            "finished_utc": finished_at.isoformat(),
            "source_code_sha256": source_digest,
            "svg_outputs_sha256": output_hash,
            "actual_svg_outputs": len(svg_payloads),
            "actual_code_cell_execution_counts": [c.execution_count for c in code_cells],
            "python_version": platform.python_version(),
            "rdkit_version": rdBase.rdkitVersion,
            "nbclient_version": nbclient.__version__,
            "input_git_commit": input_commit,
            "github_actions_run_id": os.environ.get("GITHUB_RUN_ID"),
        }
        # Write only after a complete and successful kernel run; no crafted
        # outputs or unexecuted cells are ever committed.
        temporary = path.with_suffix(".ipynb.running")
        nbformat.write(executed, temporary)
        temporary.replace(path)
        size = path.stat().st_size

        results.append({
            "number": entry["number"],
            "name": entry["name"],
            "notebook": entry["notebook"],
            "form_count": len(entry["forms"]),
            "code_cell_count": len(code_cells),
            "kernel_execution_counts": [c.execution_count for c in code_cells],
            "real_svg_outputs": len(svg_payloads),
            "real_svg_sha256": output_hash,
            "source_code_sha256": source_digest,
            "notebook_bytes": size,
            "runtime_seconds": round(elapsed, 3),
            "finished_utc": finished_at.isoformat(),
            "result": "passed",
        })
        print(
            f"EXECUTED {entry['number']:02d}/55 {entry['notebook']}: "
            f"{len(code_cells)} kernel cells, {len(svg_payloads)} real SVGs, "
            f"{size:,} bytes, {elapsed:.2f}s",
            flush=True,
        )

    report = {
        "schema": ("ccs-s5-real-kernel-execution-1.0" if full_corpus
                   else "ccs-s5-representative-kernel-smoke-1.0"),
        "status": "completed" if full_corpus else "representative_completed",
        "generation_method": "nbclient NotebookClient.execute() with ipykernel",
        "notebook_count": len(results),
        "form_count": sum(x["form_count"] for x in results),
        "actual_svg_count": sum(x["real_svg_outputs"] for x in results),
        "python_version": platform.python_version(),
        "rdkit_version": rdBase.rdkitVersion,
        "nbclient_version": nbclient.__version__,
        "started_utc": started_at.isoformat(),
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "input_git_commit": input_commit,
        "github_actions_run_id": os.environ.get("GITHUB_RUN_ID"),
        "source_fixture_sha256": manifest["data_sha256"],
        "notebooks": results,
    }
    destination = (REPORT if full_corpus else
                   (args.report_path or S5 / "smoke_execution_report.json"))
    if not full_corpus and destination.resolve() == REPORT.resolve():
        raise ValueError("Representative smoke must not overwrite the full report")
    destination.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        f"{'FULL CORPUS' if full_corpus else 'REPRESENTATIVE ONLY'}: "
        f"executed {report['notebook_count']} real Jupyter kernels "
        f"with {report['actual_svg_count']} captured real SVG outputs. "
        f"Execution report: {destination}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
