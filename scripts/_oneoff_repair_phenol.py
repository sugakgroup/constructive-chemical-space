#!/usr/bin/env python3
"""One-time S5 record-10 data correction and provenance update."""
from pathlib import Path
import hashlib
import json
import lzma
import os

root = Path(".")
data_path = root / "examples/s5/data/s5_records.json.xz"
manifest_path = root / "examples/s5/data/manifest.json"
with lzma.open(data_path, "rt", encoding="utf-8") as f:
    corpus = json.load(f)
assert len(corpus["records"]) == 55
record = next(r for r in corpus["records"] if r["number"] == 10)
assert record["name"] == "Phenol"
forms = {f["form"]: f for f in record["forms"]}
assert set(forms) == {"A", "B", "C", "D"}
assert (forms["B"]["assigned"][1], forms["B"]["assigned"][3], forms["B"]["assigned"][6]) == (7, 6, 7)
assert (forms["C"]["assigned"][1], forms["C"]["assigned"][3], forms["C"]["assigned"][6]) == (7, 6, 7)
assert record["nuclei"][1] == record["nuclei"][3] == "C" and record["nuclei"][6] == "O"
forms["C"]["assigned"][1] = 6
forms["C"]["assigned"][3] = 7
for f in record["forms"]:
    assert sum(f["assigned"]) + f["n_empty"] == 50
    spent = [0] * len(record["nuclei"])
    for i, j, ai, aj in f["allocated_connections"]:
        spent[i] += ai
        spent[j] += aj
    assert all(x <= y for x, y in zip(spent, f["assigned"]))
assert forms["B"]["assigned"] != forms["C"]["assigned"]
with lzma.open(data_path, "wt", encoding="utf-8") as f:
    json.dump(corpus, f, ensure_ascii=False, separators=(",", ":"))
new_sha = hashlib.sha256(data_path.read_bytes()).hexdigest()

manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
old_sha = "0b10a14089114d4e0aca44a99f45ca68a78bb6bea5232d0ba5304b98ba98e426"
assert manifest["data_sha256"] == old_sha
manifest["data_sha256"] = new_sha
manifest["source"]["preserves_si"] = False
manifest["source_corrections"] = [{
    "date": "2026-10-15",
    "record": 10,
    "form": "C",
    "reason": "Match para-charge-separated Form C in main-text Figure 2; SI text requires matching correction",
    "assigned_changes": {"C(1)": {"before": 7, "after": 6}, "C(3)": {"before": 6, "after": 7}},
    "unchanged": ["n_empty", "connections", "endpoint_allocations", "total_electrons"]
}]
manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

run_id = os.environ["GITHUB_RUN_ID"]
url = f"https://github.com/sugakgroup/constructive-chemical-space/actions/runs/{run_id}"
validation_path = root / "S5_VALIDATION.md"
validation = validation_path.read_text(encoding="utf-8")
assert old_sha in validation
validation = validation.replace(old_sha, new_sha)
note = (
    "### Editorial correction to phenol Form C (2026-10-15)\n\n"
    "Record 10, Form C has assigned-electron counts C(1):6 and C(3):7, "
    "instead of C(1):7 and C(3):6. Thus it represents the para-charge-separated "
    "phenol form shown in manuscript Figure 2. Form B retains the ortho "
    "assignment C(1):7. In Form C, O(6):7, all endpoint allocations, all "
    "connections, unassigned electrons, and the total of 50 electrons are "
    "unchanged. The corrected fixture is replayed and the corresponding "
    "notebook outputs are genuinely regenerated.\n\n"
    "The source-document fingerprints identify the earlier SI extraction and "
    "are retained only as historical provenance. They are NOT checksums of "
    "the corrected SI PDF; the manuscript SI Form C text must be revised "
    "accordingly before submission. For this reason source.preserves_si "
    "is set to false in the corrected manifest.\n\n"
)
assert "## Machine-readable fixture" in validation
validation = validation.replace("## Machine-readable fixture", note + "## Machine-readable fixture")
validation = validation.replace(
    "The actual notebook outputs **are committed**. The latest\n"
    "[GitHub Actions run 38391623509](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38391623509)\n"
    "regenerated all 55 notebooks after adopting the ring-aware layout and vector\n"
    "negative-charge signs.",
    "The actual notebook outputs **are committed**. The original corpus run was\n"
    "[38391623509](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38391623509);\n"
    f"the corrected v1.0.0 corpus was rerun in [run {run_id}]({url})."
)
validation_path.write_text(validation, encoding="utf-8")

readme_path = root / "README.md"
readme = readme_path.read_text(encoding="utf-8")
old_line = (
    "The last full execution was [GitHub Actions run 38391623509]"
    "(https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38391623509)"
    " (Python 3.11.17, RDKit 2026.09.1, nbclient 0.11.0)."
)
assert old_line in readme
readme = readme.replace(old_line,
    f"The corrected v1.0.0 corpus was re-executed in [GitHub Actions run {run_id}]({url}).")
readme = readme.replace(
    "The audited post-merge implementation passed [550 tests]"
    "(https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38393815432).",
    f"The corrected v1.0.0 corpus was re-tested in [run {run_id}]({url})."
)
readme_path.write_text(readme, encoding="utf-8")

s5readme_path = root / "examples/s5/README.md"
s5readme = s5readme_path.read_text(encoding="utf-8")
old_line = (
    "All 55 notebooks were last regenerated in "
    "[run 38391623509](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38391623509),"
    " which passed 494 validation tests before and after output commit "
    "[5ef1bfb](https://github.com/sugakgroup/constructive-chemical-space/commit/5ef1bfb6278331c248d109155472e4a9109f1968)."
)
# Use a less formatting-sensitive replacement, preserving its historical link.
marker = "All 55 notebooks were last regenerated in "
start = s5readme.find(marker)
end = s5readme.find("\n\n", start)
assert start >= 0 and end > start
s5readme = (s5readme[:start] +
    f"All 55 corrected v1.0.0 notebooks were regenerated in [run {run_id}]({url}). "
    "The original pre-correction execution was "
    "[run 38391623509](https://github.com/sugakgroup/constructive-chemical-space/actions/runs/38391623509)." +
    s5readme[end:])
s5readme_path.write_text(s5readme, encoding="utf-8")

test_path = root / "tests/test_s5_phenol_form_c.py"
assert not test_path.exists()
test_path.write_text(
    '"""Regression: phenol Form C represents para-charge separation in Figure 2."""\n\n'
    'from examples.s5._common import build_example, load_record, verify_example\n\n\n'
    'def test_record10_ortho_and_para_electron_assignments():\n'
    '    record = load_record(10)\n'
    '    forms = {f["form"]: f for f in record["forms"]}\n'
    '    assert record["nuclei"][1] == record["nuclei"][3] == "C"\n'
    '    assert record["nuclei"][6] == "O"\n'
    '    assert (forms["B"]["assigned"][1], forms["B"]["assigned"][3], forms["B"]["assigned"][6]) == (7, 6, 7)\n'
    '    assert (forms["C"]["assigned"][1], forms["C"]["assigned"][3], forms["C"]["assigned"][6]) == (6, 7, 7)\n'
    '    assert forms["B"]["assigned"] != forms["C"]["assigned"]\n'
    '    assert all(sum(f["assigned"]) + f["n_empty"] == 50 for f in forms.values())\n'
    '    assert verify_example(build_example(10)) == {"A": 26, "B": 26, "C": 26, "D": 26}\n',
    encoding="utf-8"
)
print("Phenol Form C fixed. Dataset SHA256:", new_sha)
