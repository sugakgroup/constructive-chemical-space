"""Regression: phenol Form C represents para-charge separation in Figure 2."""

from examples.s5._common import build_example, load_record, verify_example


def test_record10_ortho_and_para_electron_assignments():
    record = load_record(10)
    forms = {f["form"]: f for f in record["forms"]}
    assert record["nuclei"][1] == record["nuclei"][3] == "C"
    assert record["nuclei"][6] == "O"
    assert (forms["B"]["assigned"][1], forms["B"]["assigned"][3], forms["B"]["assigned"][6]) == (7, 6, 7)
    assert (forms["C"]["assigned"][1], forms["C"]["assigned"][3], forms["C"]["assigned"][6]) == (6, 7, 7)
    assert forms["B"]["assigned"] != forms["C"]["assigned"]
    assert all(sum(f["assigned"]) + f["n_empty"] == 50 for f in forms.values())
    assert verify_example(build_example(10)) == {"A": 26, "B": 26, "C": 26, "D": 26}
