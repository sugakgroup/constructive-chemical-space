"""End-to-end validation of every SI S5 example notebook and its CCS Form(s)."""

from pathlib import Path
import json
import re

import pytest

from examples.s5._common import (
    S5_ROOT, build_example, load_manifest, load_record,
    summarize_example, verify_example,
)
from constructive_chemical_space import ElectronAllocation

MANIFEST = load_manifest()
ENTRIES = MANIFEST["records"]
FORM_CASES = [(e["number"], label) for e in ENTRIES for label in e["forms"]]


def test_manifest_complete_and_single_source_of_truth():
    assert MANIFEST["schema"] == "ccs-s5-notebooks-1.0"
    assert MANIFEST["record_count"] == 55
    assert MANIFEST["form_count"] == 84
    assert [e["number"] for e in ENTRIES] == list(range(1, 56))
    assert sum(len(e["forms"]) for e in ENTRIES) == 84
    assert len({e["notebook"] for e in ENTRIES}) == 55
    assert len(MANIFEST["source"]["pdf_sha256"]) == 64
    assert len(MANIFEST["source"]["docx_sha256"]) == 64
    assert MANIFEST["data_sha256"] == "1a9c1dcd119159aba4a95ff9941ea73110de18f97583200e33b3578df6aad8dd"
    assert not (S5_ROOT.parents[1] / "data" / "s5_records.json.xz").exists()


@pytest.mark.parametrize("number", range(1, 56))
def test_each_notebook_structure_and_manifest_matches_source(number):
    entry = ENTRIES[number - 1]
    original = load_record(number)
    assert entry["name"] == original["name"]
    assert entry["section"] == original["section"]
    assert entry["expected"] == original["expected"]
    assert entry["forms"] == [x["form"] for x in original["forms"]]
    assert entry["source"] == original["source"]

    name = entry["notebook"]
    assert re.fullmatch(rf"S5_{number:02d}_[A-Za-z0-9_]+\.ipynb", name)
    notebook_file = S5_ROOT / "notebooks" / name
    assert notebook_file.is_file()
    obj = json.loads(notebook_file.read_text(encoding="utf-8"))
    assert (obj["nbformat"], obj["nbformat_minor"]) == (4, 5)
    assert obj["metadata"]["ccs_s5_record"] == number
    assert obj["cells"][0]["cell_type"] == "markdown"
    assert entry["name"] in "".join(obj["cells"][0]["source"])
    source_cells = ["".join(cell["source"]) for cell in obj["cells"]
                    if cell["cell_type"] == "code"]
    assert len(source_cells) == 6
    assert f"build_example({number})" in source_cells[0]
    assert "verify_example(example)" in source_cells[0]
    assert "from constructive_chemical_space import Word, Nuc, Ele, Con" in source_cells[1]
    assert "word_A = Word([" in source_cells[1]
    assert "to_word(" not in source_cells[1]  # no runtime auto-generated words
    assert "Nuc(" in source_cells[1] or number == 46  # electron-only boundary
    assert "Con(" in source_cells[1] or entry["expected"]["connections"] == 0
    for letter in entry["forms"]:
        assert f"word_{letter} = Word([" in source_cells[1]
        assert f"state_{letter} = word_{letter}.replay_state()" in source_cells[2]
        assert f"allocation_{letter} = state_{letter}.form" in source_cells[2]
        assert f"constitution_{letter} = allocation_{letter}.to_constitution()" in source_cells[2]
        assert f"allocation_{letter}.show(" in source_cells[4]
    for letter in entry["forms"]:
        assert f"constitution_{letter}.show(" in source_cells[3]
        assert f"constitution_{letter}.save(" in source_cells[5]
        assert f"signature_{letter} == reference_constitution_signature" in source_cells[2]
    assert "allocation_label_mode='endpoint'" in source_cells[3]
    assert "EXPORT = False" in source_cells[5]
    # During the transition to executed notebooks, both clean source files
    # and fully kernel-executed files are valid; never fabricate output.
    provenance = obj["metadata"].get("ccs_execution")
    display_count = 0
    for cell in obj["cells"]:
        if cell["cell_type"] != "code":
            continue
        compile("".join(cell["source"]), str(notebook_file), "exec")
        if provenance is None:
            assert cell["execution_count"] is None
            assert cell["outputs"] == []
            continue
        assert isinstance(cell["execution_count"], int)
        assert cell["execution_count"] > 0
        for output in cell["outputs"]:
            if output["output_type"] in {"display_data", "execute_result"}:
                if "image/svg+xml" in output.get("data", {}):
                    display_count += 1
    if provenance is not None:
        assert provenance["status"] == "completed"
        assert provenance["method"] == "nbclient real ipykernel execution"
        assert len(provenance["source_code_sha256"]) == 64
        assert len(provenance["svg_outputs_sha256"]) == 64
        assert display_count == provenance["actual_svg_outputs"]
        assert display_count == 2 * entry["expected"]["forms"]


@pytest.mark.parametrize("number", range(1, 56))
def test_s5_example_full_compatibility_and_counts(number):
    example = build_example(number)
    original = load_record(number)
    expected = original["expected"]
    assert len(example.order) == expected["nuclei"]
    assert len(example.constitution.connections) == expected["connections"]
    assert example.constitution.num_electrons == expected["electrons"]
    assert example.constitution.total_charge == expected["charge"]
    assert len(example.forms) == expected["forms"]
    assert "Forms:" in summarize_example(example)
    word_lengths = verify_example(example)
    assert list(word_lengths) == [x["form"] for x in original["forms"]]


@pytest.mark.parametrize("number,letter", FORM_CASES,
                         ids=[f"S5-{n:02d}-{letter}" for n, letter in FORM_CASES])
def test_exact_WD_word_replay_for_all_84_forms(number, letter):
    example = build_example(number)
    form = example.forms[letter]
    assert isinstance(form, ElectronAllocation)
    word = form.to_word(order=example.order)
    assert len(word.tokens) >= len(example.order)
    assert word.tokens == form.to_word(order=example.order).tokens
    # Completeness including all endpoint allocations is checked in verify_example.
    assert form.num_electrons == example.constitution.num_electrons


@pytest.mark.parametrize("number", [1, 2, 9, 10, 14, 17, 34, 35, 46])
def test_s5_intrinsic_equality_under_nuclear_permutation(number):
    original = load_record(number)
    example = build_example(number)
    permutation = tuple(reversed(range(len(example.order))))
    mapping = {old: new for new, old in enumerate(permutation)}
    for f in original["forms"]:
        d = example.forms[f["form"]]
        changed = ElectronAllocation.from_indexed(
            [original["nuclei"][i] for i in permutation],
            [f["assigned"][i] for i in permutation],
            [(mapping[i], mapping[j], a, b)
             for i, j, a, b in f["allocated_connections"]],
            num_unassigned_electrons=f["n_empty"],
        )
        assert changed == d
        assert changed.to_constitution() == example.constitution


def test_s5_critical_boundary_examples():
    assert [(r["number"], f["form"]) for r in
            (load_record(i) for i in range(1, 56)) for f in r["forms"]
            if f["n_empty"] > 0] == [(43, "A"), (46, "A"), (48, "A")]
    assert load_record(1)["forms"][1]["allocated_connections"] == [[0, 1, 0, 0]]
    assert load_record(45)["connections"] == []
    assert load_record(46)["nuclei"] == []
    assert load_record(46)["forms"][0]["n_empty"] == 1
    for a, b in [(39, 40), (47, 48), (49, 50), (51, 52)]:
        assert build_example(a).constitution != build_example(b).constitution


def test_example_notebook_execution_setup_without_rendering():
    # Execute source loading and complete replay-check cells from representative
    # notebooks, without requiring notebook display or a graphical backend.
    for number in (1, 2, 17, 42, 46, 55):
        entry = ENTRIES[number - 1]
        notebook = json.loads(
            (S5_ROOT / "notebooks" / entry["notebook"]).read_text(encoding="utf-8")
        )
        code = ["".join(c["source"]) for c in notebook["cells"]
                if c["cell_type"] == "code"]
        namespace = {}
        for source in code[:3]:
            exec(compile(source, entry["notebook"], "exec"), namespace)
        assert namespace["example"].number == number
        assert namespace["replay_lengths"] == verify_example(namespace["example"])
        for letter in entry["forms"]:
            state = namespace[f"state_{letter}"]
            allocation = namespace[f"allocation_{letter}"]
            constitution = namespace[f"constitution_{letter}"]
            assert state.form is allocation
            assert allocation.to_constitution() == constitution or len(state.address_table) > 15
            assert constitution.num_electrons == entry["expected"]["electrons"]


def test_render_and_export_small_examples(tmp_path):
    pytest.importorskip("rdkit")
    example = build_example(17)  # F-H-F, the explicit-H bridge must survive.
    output = example.constitution.to_svg(order=example.order)
    assert "N_e = 20" in output
    assert "lossy view" not in output
    from examples.s5._common import export_example
    saved = export_example(example, output_dir=tmp_path, formats=("svg",))
    assert len(saved) == 1 + len(example.forms)
    assert all(path.read_text().startswith("<?xml") for path in saved)


def test_invalid_record_number_is_rejected():
    for invalid in (-1, 0, 56, True, 1.5, "1"):
        with pytest.raises(ValueError):
            load_record(invalid)


@pytest.mark.parametrize("number", range(1, 56),
                         ids=[f"S5-{i:02d}" for i in range(1, 56)])
def test_every_s5_notebook_can_render_all_forms(number):
    """Integration check of every explicit-H figure used by the 55 notebooks."""
    pytest.importorskip("rdkit")
    import xml.etree.ElementTree as ET

    example = build_example(number)
    width = 1500 if len(example.order) > 60 else 950
    height = 1000 if len(example.order) > 60 else 650
    objects = [example.constitution, *example.forms.values()]
    for obj in objects:
        svg = obj.to_svg(
            order=example.order, show_hydrogens=True,
            show_indices=True, show_total_charge=True,
            width=width, height=height,
        )
        ET.fromstring(svg)  # Well-formed SVG, not just a nonempty string.
        assert f"N_e = {example.constitution.num_electrons}" in svg
        assert "lossy view" not in svg


def test_notebook_execution_report_when_available():
    """Never accept an incomplete run report as evidence of execution."""
    report_path = S5_ROOT / "execution_report.json"
    if not report_path.exists():
        return  # Notebooks have not yet been executed; no fake report.
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["schema"] == "ccs-s5-real-kernel-execution-1.0"
    assert report["status"] == "completed"
    assert report["notebook_count"] == 55
    assert report["form_count"] == 84
    assert report["actual_svg_count"] == 168
    assert len(report["notebooks"]) == 55
    assert [x["number"] for x in report["notebooks"]] == list(range(1, 56))
    for entry in report["notebooks"]:
        notebook = json.loads(
            (S5_ROOT / "notebooks" / entry["notebook"]).read_text(encoding="utf-8")
        )
        assert notebook["metadata"]["ccs_execution"]["svg_outputs_sha256"] == entry["real_svg_sha256"]
        assert notebook["metadata"]["ccs_execution"]["actual_svg_outputs"] == entry["real_svg_outputs"]
        assert notebook["metadata"]["ccs_execution"]["actual_code_cell_execution_counts"] == entry["kernel_execution_counts"]


@pytest.mark.parametrize("number", range(1, 56))
def test_all_literal_notebook_words_replay_into_full_ccs_objects(number):
    """Execute literal Word definitions, not an indirect to_word() witness."""
    entry = ENTRIES[number - 1]
    path = S5_ROOT / "notebooks" / entry["notebook"]
    nb = json.loads(path.read_text(encoding="utf-8"))
    code = ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]
    scope = {}
    for block in code[:3]:  # data, explicit words, explicit replay/assertions
        exec(compile(block, str(path), "exec"), scope)
    example = scope["example"]
    assert len(example.forms) == len(entry["forms"])
    for label in entry["forms"]:
        word = scope[f"word_{label}"]
        state = scope[f"state_{label}"]
        form = scope[f"allocation_{label}"]
        constitution = scope[f"constitution_{label}"]
        assert state.form is form
        assert form.to_word(order=state.address_table).tokens == word.tokens
        assert constitution.num_electrons == example.constitution.num_electrons
        assert constitution.total_charge == example.constitution.total_charge
        assert len(constitution.nuclei) == len(example.constitution.nuclei)
        assert len(constitution.connections) == len(example.constitution.connections)
        from examples.s5._common import (
            indexed_constitution_signature, indexed_allocation_signature,
        )
        # Exact construction/replay equality for every S5 Form, including
        # large symmetric graphs where general graph isomorphism is expensive.
        assert indexed_constitution_signature(
            constitution, state.address_table
        ) == indexed_constitution_signature(
            example.constitution, example.order
        )
        assert indexed_allocation_signature(
            form, state.address_table
        ) == indexed_allocation_signature(
            example.forms[label], example.order
        )


@pytest.mark.parametrize("number", range(1, 56))
def test_saved_notebook_output_hashes_match_kernel_report(number):
    """Re-hash REAL captured SVGs and source code, not fabricated outputs."""
    import hashlib
    import xml.etree.ElementTree as ET

    report_path = S5_ROOT / "execution_report.json"
    if not report_path.is_file():
        pytest.skip("No real-kernel execution report has been committed yet")

    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["schema"] == "ccs-s5-real-kernel-execution-1.0"
    assert report["status"] == "completed"
    assert report["generation_method"] == "nbclient NotebookClient.execute() with ipykernel"
    assert report["notebook_count"] == 55
    assert report["form_count"] == 84
    assert report["actual_svg_count"] == 168

    entry = ENTRIES[number - 1]
    measured = report["notebooks"][number - 1]
    assert measured["number"] == number
    assert measured["notebook"] == entry["notebook"]
    assert measured["result"] == "passed"
    notebook_path = S5_ROOT / "notebooks" / entry["notebook"]
    nb = json.loads(notebook_path.read_text(encoding="utf-8"))
    meta = nb["metadata"]["ccs_execution"]
    assert meta["status"] == "completed"
    assert meta["method"] == "nbclient real ipykernel execution"
    assert meta["github_actions_run_id"] == report["github_actions_run_id"]
    assert meta["input_git_commit"] == report["input_git_commit"]

    cells = [c for c in nb["cells"] if c["cell_type"] == "code"]
    assert len(cells) == 6
    execution_counts = [c["execution_count"] for c in cells]
    assert execution_counts == [1, 2, 3, 4, 5, 6]
    assert execution_counts == measured["kernel_execution_counts"]
    source_text = "\n\n".join("".join(c["source"]) for c in cells)
    source_digest = hashlib.sha256(source_text.encode("utf-8")).hexdigest()
    assert source_digest == meta["source_code_sha256"]
    assert source_digest == measured["source_code_sha256"]

    payloads = []
    for c in cells:
        for output in c["outputs"]:
            assert output["output_type"] != "error"
            if output["output_type"] not in {"display_data", "execute_result"}:
                continue
            svg = output.get("data", {}).get("image/svg+xml")
            if svg is None:
                continue
            svg = "".join(svg) if isinstance(svg, list) else svg
            ET.fromstring(svg)
            payloads.append(svg)
    assert len(payloads) == 2 * len(entry["forms"])
    data_digest = hashlib.sha256("\n".join(payloads).encode("utf-8")).hexdigest()
    assert data_digest == meta["svg_outputs_sha256"]
    assert data_digest == measured["real_svg_sha256"]
