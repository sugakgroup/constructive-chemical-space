#!/usr/bin/env python3
"""Idempotent, author-specified correction for S5 Record 13 Forms B and C.

The prior transcript accidentally moved the (2,2) allocation from (9,14)
onto (13,14). The author's explicit SI specification instead retains the
correct original (1,1) on (13,14) and puts (2,2) on (9,14).

Only these two allocation fields in Forms B/C are changed. Form A, all
nuclear identities, assigned counts, connections, n_empty and every other
record remain unchanged. This tool never fabricates notebook outputs.
"""
from __future__ import annotations

import hashlib
import json
import lzma
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "data" / "manifest.json"


def main() -> None:
    index = json.loads(MANIFEST.read_text(encoding="utf-8"))
    file = ROOT / "data" / index["data_file"]
    old_bytes = file.read_bytes()
    if hashlib.sha256(old_bytes).hexdigest() != index["data_sha256"]:
        raise RuntimeError("Refusing to edit an unverified S5 dataset")

    data = json.loads(lzma.decompress(old_bytes).decode("utf-8"))
    record = data["records"][12]
    if record["number"] != 13 or record["expected"]["forms"] != 3:
        raise RuntimeError("Unexpected Record 13 fixture structure")

    changed = 0
    for form in record["forms"]:
        if form["form"] == "A":
            continue
        if form["form"] not in {"B", "C"}:
            raise RuntimeError("Unexpected Record 13 form")
        if form["form"] == "C" and (
            form["assigned"][0] != 7 or form["assigned"][7] != 5
        ):
            raise RuntimeError("Form C has unexpected assigned-electron counts")

        mapping = {(i, j): (a, b) for i, j, a, b
                   in form["allocated_connections"]}
        if len(mapping) != len(form["allocated_connections"]):
            raise RuntimeError("Repeated connection in Form " + form["form"])

        desired = {(9, 14): (2, 2), (13, 14): (1, 1)}
        old = {(9, 14): (1, 1), (13, 14): (2, 2)}
        if all(mapping[k] == desired[k] for k in desired):
            continue
        if not all(mapping[k] == old[k] for k in desired):
            raise RuntimeError(
                "Unexpected allocations: refusing to overwrite manual edits"
            )
        form["allocated_connections"] = [
            [i, j, *(desired.get((i, j), (a, b)))]
            for i, j, a, b in form["allocated_connections"]
        ]
        changed += 1

    if changed == 0:
        print("Record 13 already matches the author-specified Form B/C data")
        return
    if changed != 2:
        raise RuntimeError("Partial Record 13 update refused")

    normalized = json.dumps(
        data, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    compressed = lzma.compress(normalized, preset=6)
    file.write_bytes(compressed)
    index["data_sha256"] = hashlib.sha256(compressed).hexdigest()
    MANIFEST.write_text(
        json.dumps(index, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        "Corrected Record 13 Forms B/C only; "
        f"new fixture SHA-256 {index['data_sha256']}"
    )


if __name__ == "__main__":
    main()
