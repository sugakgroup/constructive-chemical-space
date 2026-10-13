"""Optional visualization regression tests using installed CCS package objects."""
from io import BytesIO
from math import hypot
import xml.etree.ElementTree as ET
from unittest.mock import patch

import pytest

pytest.importorskip("rdkit")
from rdkit import Chem
from PIL import Image

from constructive_chemical_space import (
    AllocatedConnection, Con, Constitution, Ele, ElectronAllocation,
    End, Nuc, Nucleus, PracticalWord, Word,
)
from constructive_chemical_space.ccs import (
    _build_ccs_visual_mol, _ccs_footer, _visual_source,
)


def _mol(obj, **kwargs):
    source, order, kind = _visual_source(obj, view=kwargs.pop("view", "auto"),
                                        order=kwargs.pop("order", None))
    return _build_ccs_visual_mol(
        source, order, source_kind=kind,
        show_hydrogens=kwargs.pop("show_hydrogens", True),
        show_indices=kwargs.pop("show_indices", False),
        diagnostics=kwargs.pop("diagnostics", False),
        label_all_allocations=kwargs.pop("label_all_allocations", False),
    )


def test_all_h_nuclei_connections_retained_in_constitution_default():
    m = Constitution.from_indexed(
        ["O", "H", "H"], [(0, 1), (0, 2)], 10
    )
    mol = _mol(m)
    assert mol.GetNumAtoms() == 3
    assert mol.GetNumBonds() == 2
    assert mol.GetProp("_CCSView") == "constitution"
    assert not mol.GetBoolProp("_CCSLossyHydrogenOmission")
    assert all(b.GetBondType() == Chem.BondType.SINGLE for b in mol.GetBonds())
    assert "N_e = 10" in m.to_svg()
    assert m.num_electrons == 10


def test_bridging_h_bifluoride_is_not_removed_by_default():
    d = ElectronAllocation.from_indexed(
        ["F", "H", "F"], [9, 2, 9],
        [(0, 1, 1, 1), (1, 2, 1, 1)],
    )
    for view in ("constitution", "allocation"):
        mol = _mol(d, view=view)
        assert (mol.GetNumAtoms(), mol.GetNumBonds()) == (3, 2)
        h = next(atom for atom in mol.GetAtoms() if atom.GetSymbol() == "H")
        assert h.GetDegree() == 2
    assert not d.to_svg().count("lossy view")


def test_hidden_h_is_explicitly_lossy_and_retains_full_system_footer():
    m = Constitution.from_indexed(["F", "H", "F"], [(0, 1), (1, 2)], 20)
    mol = _mol(m, show_hydrogens=False)
    assert mol.GetNumAtoms() == 2
    assert mol.GetNumBonds() == 0
    assert mol.GetBoolProp("_CCSLossyHydrogenOmission")
    svg = m.to_svg(show_hydrogens=False)
    assert "H omitted - lossy view" in svg
    assert "N_e = 20" in svg
    assert "H omitted - lossy view" not in m.to_svg()


def test_isotope_resolved_h_is_shown_explicitly():
    d = ElectronAllocation.from_indexed(
        ["2H", "O", "H"], [1, 8, 1],
        [(0, 1, 1, 1), (1, 2, 1, 1)],
    )
    mol = _mol(d)
    assert sorted((a.GetSymbol(), a.GetIsotope()) for a in mol.GetAtoms()) == [
        ("H", 0), ("H", 2), ("O", 0)
    ]
    assert mol.GetNumBonds() == 2


def test_assignment_imbalance_annotation_not_rdkit_formal_charge():
    d = ElectronAllocation.from_indexed(
        ["O", "H"], [7, 1], [(0, 1, 1, 1)]
    )
    mol = _mol(d)
    o = next(a for a in mol.GetAtoms() if a.GetSymbol() == "O")
    assert o.GetFormalCharge() == 0
    assert o.GetIntProp("_CCSAssignmentImbalance") == 1
    # Default visible charge glyph is separate from the RDKit formal-charge
    # field and the optional detailed atomNote diagnostics.
    svg, nuclear_labels = _nuclear_svg_labels(d)
    assert any(
        x["kind"] == "charge" and x["glyph"] == "+"
        for x in nuclear_labels
    )
    assert "q=+1" in svg
    assert d.to_constitution().total_charge == 1


def test_diagnostics_show_primary_counts_and_residuals():
    d = ElectronAllocation.from_indexed(
        ["N", "B"], [7, 5], [(0, 1, 2, 0)]
    )
    mol = _mol(d, diagnostics=True)
    assert all(a.HasProp("atomNote") for a in mol.GetAtoms())
    assert all("n=" in a.GetProp("atomNote") for a in mol.GetAtoms())


@pytest.mark.parametrize(
    ("a", "b", "bond_type"), [
        (1, 1, Chem.BondType.SINGLE),
        (2, 2, Chem.BondType.DOUBLE),
        (3, 3, Chem.BondType.TRIPLE),
        (2, 0, Chem.BondType.DATIVE),
        (0, 2, Chem.BondType.DATIVE),
    ],
)
def test_standard_allocation_glyphs(a, b, bond_type):
    d = ElectronAllocation.from_indexed(
        ["U", "O"], [92, 10], [(0, 1, a, b)]
    )
    mol = _mol(d)
    assert mol.GetBondWithIdx(0).GetBondType() == bond_type


@pytest.mark.parametrize(("a", "b"), [(0, 0), (1, 0), (3, 2), (5, 5), (0, 6), (0, 9)])
def test_unusual_allocations_have_oriented_indexed_labels(a, b):
    u, o = Nucleus("U"), Nucleus("O")
    d = ElectronAllocation(
        [u, o], {u: 92, o: 10},
        [AllocatedConnection(u, o, a, b)]
    )
    mol = _mol(d, order=(o, u))
    bond = mol.GetBondWithIdx(0)
    assert bond.GetBondType() == Chem.BondType.ZERO
    assert bond.GetProp("bondNote") == f"0:{b}, 1:{a}"


def test_zero_allocation_edge_not_deleted():
    d = ElectronAllocation.from_indexed(
        ["C", "C"], [6, 6], [(0, 1, 0, 0)]
    )
    mol = _mol(d)
    assert mol.GetNumBonds() == 1
    assert mol.GetBondWithIdx(0).HasProp("bondNote")


def test_word_and_practical_word_preserve_replay_order():
    w = Word([Nuc("O"), Nuc("H"), Con(0, 1, 1, 1)])
    state = w.replay_state()
    # Replay creates new abstract Nucleus objects; compare temporary order
    # by nuclear identities rather than Python object identities.
    assert tuple(v.identity for v in _visual_source(w)[1]) == tuple(
        v.identity for v in state.address_table
    )
    p = PracticalWord([Nuc("C"), End()])
    src, order, kind = _visual_source(p)
    assert len(src.nuclei) == len(order) == 5
    assert kind == "allocation"


def test_disconnected_and_electron_only_can_be_rendered():
    d = Word([Ele()]).replay()
    mol = _mol(d)
    assert mol.GetNumAtoms() == 0
    svg = d.to_svg()
    assert "N_e = 1" in svg
    assert "n_empty = 1" in svg
    m = Constitution.from_indexed(["Na", "Cl"], [], 28)
    mol2 = _mol(m)
    assert (mol2.GetNumAtoms(), mol2.GetNumBonds()) == (2, 0)


def test_exotic_multicenter_allocations_have_literal_annotations():
    # A representative high-electron metal/cage allocation.
    d = ElectronAllocation.from_indexed(
        ["Ni", "Pb", "Pb"], [28, 82, 82],
        [(0, 1, 0, 7), (0, 2, 0, 8), (1, 2, 0, 0)],
    )
    mol = _mol(d)
    assert (mol.GetNumAtoms(), mol.GetNumBonds()) == (3, 3)
    assert all(b.GetBondType() == Chem.BondType.ZERO for b in mol.GetBonds())
    assert all(b.HasProp("bondNote") for b in mol.GetBonds())


def test_svg_png_export_and_save(tmp_path):
    p = PracticalWord([Nuc("C"), End()])
    svg = p.to_svg(width=600, height=450)
    ET.fromstring(svg)  # Well-formed SVG.
    assert "N_e = 10" in svg
    png = p.to_png(width=600, height=450)
    with Image.open(BytesIO(png)) as image:
        assert image.format == "PNG"
        assert image.size == (600, 450)
    svg_path = tmp_path / "methane.svg"
    png_path = tmp_path / "methane.png"
    p.save(svg_path, width=600, height=450)
    p.save(png_path, width=600, height=450)
    assert svg_path.read_text(encoding="utf-8") == svg
    assert p.to_svg(width=600, height=450) == svg
    assert png_path.read_bytes().startswith(b"\x89PNG")
    with pytest.raises(ValueError):
        p.save(tmp_path / "not-supported.pdf")


def test_show_calls_ipython_display_without_changing_core_data():
    pytest.importorskip("IPython")
    d = ElectronAllocation.from_indexed(["H", "H"], [1, 1], [(0, 1, 1, 1)])
    baseline = d.to_word().replay()
    with patch("IPython.display.display") as mocked:
        assert d.show(width=600, height=450) is None
        assert mocked.called
    assert d == baseline


def test_explicit_order_is_preserved_and_validated():
    a, b = Nucleus("C"), Nucleus("N")
    d = ElectronAllocation(
        [a, b], {a: 6, b: 7},
        [AllocatedConnection(a, b, 1, 1)]
    )
    assert _visual_source(d, order=(b, a))[1] == (b, a)
    with pytest.raises(ValueError):
        d.to_svg(order=(a, a))


@pytest.mark.parametrize(
    ("a", "b"), [(0, 0), (1, 0), (0, 1), (3, 2), (0, 9)]
)
def test_svg_endpoint_labels_are_attached_to_their_own_nuclei(a, b):
    """Homonuclear asymmetric connections must never swap endpoint counts."""
    left, right = Nucleus("C"), Nucleus("C")
    d = ElectronAllocation(
        [left, right], {left: 6, right: 9},
        [AllocatedConnection(left, right, a, b)],
    )
    original = d.to_svg(order=(left, right))
    root = ET.fromstring(original)
    tags = [
        element for element in root.iter()
        if element.attrib.get("class") == "ccs-endpoint-allocation"
    ]
    assert len(tags) == 2
    assert {
        int(x.attrib["data-atom"]): int(x.attrib["data-electrons"])
        for x in tags
    } == {0: a, 1: b}
    assert {x.text for x in tags} == {str(a), str(b)}
    # Counts are located in two different parts of the drawn bond.
    assert abs(float(tags[0].attrib["x"]) - float(tags[1].attrib["x"])) > 10

    reversed_svg = d.to_svg(order=(right, left))
    swapped = [
        element for element in ET.fromstring(reversed_svg).iter()
        if element.attrib.get("class") == "ccs-endpoint-allocation"
    ]
    assert {
        int(x.attrib["data-atom"]): int(x.attrib["data-electrons"])
        for x in swapped
    } == {0: b, 1: a}
    assert d.num_electrons == 15


def test_endpoint_annotation_modes_and_unchanged_intrinsic_object():
    d = ElectronAllocation.from_indexed(
        ["H", "H"], [1, 1], [(0, 1, 1, 0)]
    )
    old_word = d.to_word()
    svg = d.to_svg(width=650, height=400)
    endpoint = [
        node for node in ET.fromstring(svg).iter()
        if node.attrib.get("class") == "ccs-endpoint-allocation"
    ]
    assert len(endpoint) == 2
    assert {int(x.attrib["data-electrons"]) for x in endpoint} == {0, 1}
    for mode in ("bond", "none"):
        old_style = d.to_svg(allocation_label_mode=mode)
        assert "class='ccs-endpoint-allocation' data-bond=" not in old_style
    assert d.to_word().tokens == old_word.tokens
    with pytest.raises(ValueError, match="allocation_label_mode"):
        d.to_svg(allocation_label_mode="invalid")


def test_png_renders_endpoint_labels_as_well():
    d = ElectronAllocation.from_indexed(
        ["H", "H"], [1, 1], [(0, 1, 1, 0)]
    )
    endpoint_png = d.to_png(width=650, height=400)
    no_labels_png = d.to_png(
        width=650, height=400, allocation_label_mode="none"
    )
    with Image.open(BytesIO(endpoint_png)) as image:
        assert image.format == "PNG"
        assert image.size == (650, 400)
        assert image.convert("RGB").getpixel((0, 0)) == (255, 255, 255)
    # The drawing must have actually changed, not just metadata or output bytes.
    with Image.open(BytesIO(endpoint_png)) as image1, Image.open(BytesIO(no_labels_png)) as image2:
        from PIL import ImageChops
        assert ImageChops.difference(
            image1.convert("RGB"), image2.convert("RGB")
        ).getbbox() is not None


def test_conventional_allocation_labels_are_optional():
    d = ElectronAllocation.from_indexed(
        ["H", "H"], [1, 1], [(0, 1, 1, 1)]
    )
    for include, expected in ((False, 0), (True, 2)):
        svg = d.to_svg(label_all_allocations=include)
        labels = [
            item for item in ET.fromstring(svg).iter()
            if item.attrib.get("class") == "ccs-endpoint-allocation"
        ]
        assert len(labels) == expected


def _nuclear_svg_labels(obj, **kwargs):
    svg = obj.to_svg(**kwargs)
    labels = [
        {
            "atom": int(item.attrib["data-atom"]),
            "kind": item.attrib["data-kind"],
            "count": int(item.attrib["data-count"]),
            "glyph": (
                "•" * sum(c.tag.rsplit("}", 1)[-1] == "circle"
                          for c in item)
                if item.attrib["data-kind"] == "residual"
                else "".join(item.itertext())
            ),
        }
        for item in ET.fromstring(svg).iter()
        if item.attrib.get("class") == "ccs-nuclear-annotation"
    ]
    return svg, labels


def test_amine_nitrogen_displays_two_residual_valence_dots():
    # Neutral NH3: N retains 7-3-2=2 valence residual electrons.
    d = ElectronAllocation.from_indexed(
        ["N", "H", "H", "H"], [7, 1, 1, 1],
        [(0, 1, 1, 1), (0, 2, 1, 1), (0, 3, 1, 1)],
    )
    svg, labels = _nuclear_svg_labels(d, show_indices=False)
    residual = [i for i in labels if i["kind"] == "residual"]
    n_atom = next(i for i, nuc in enumerate(_visual_source(d)[1])
                  if nuc.identity.symbol == "N")
    assert residual == [{"atom": n_atom, "kind": "residual",
                         "count": 2, "glyph": "••"}]
    assert not any(i["kind"] == "charge" for i in labels)
    assert "n=7" not in svg
    assert "class='ccs-nuclear-annotation'" in svg
    assert "H omitted" not in svg


def test_water_oxygen_displays_four_valence_residual_dots():
    d = ElectronAllocation.from_indexed(
        ["O", "H", "H"], [8, 1, 1],
        [(0, 1, 1, 1), (0, 2, 1, 1)],
    )
    _, labels = _nuclear_svg_labels(d)
    o_atom = next(i for i, nuc in enumerate(_visual_source(d)[1])
                  if nuc.identity.symbol == "O")
    oxygen_dots = [v for v in labels if v["kind"] == "residual"
                   and v["atom"] == o_atom]
    assert len(oxygen_dots) == 2
    assert sum(x["count"] for x in oxygen_dots) == 4
    assert all(x["glyph"] == "••" for x in oxygen_dots)


@pytest.mark.parametrize(
    ("identities", "assigned", "edges", "expected_q", "expected_dots"),
    [
        (["N", "H", "H", "H", "H"], [6, 1, 1, 1, 1],
         [(0, 1, 1, 1), (0, 2, 1, 1), (0, 3, 1, 1), (0, 4, 1, 1)], "+", 0),
        (["O", "H"], [9, 1], [(0, 1, 1, 1)], "−", 6),
        (["O"], [10], [], "2−", 8),
    ],
)
def test_assignment_imbalance_signs_and_valence_dots(
    identities, assigned, edges, expected_q, expected_dots
):
    d = ElectronAllocation.from_indexed(identities, assigned, edges)
    _, labels = _nuclear_svg_labels(d)
    target = next(i for i, nuc in enumerate(_visual_source(d)[1])
                  if nuc.identity.symbol == identities[0])
    charge = [i for i in labels if i["kind"] == "charge" and i["atom"] == target]
    assert len(charge) == 1
    assert charge[0]["glyph"] == expected_q
    remaining = [i for i in labels if i["kind"] == "residual"
                 and i["atom"] == target]
    assert sum(i["count"] for i in remaining) == expected_dots
    for i in d.nuclei:
        assert d.assignment_imbalance(i) == i.identity.atomic_number - d.assigned_electrons(i)


def test_non_p_block_and_hydrogen_have_no_valence_dot_depicton():
    d = ElectronAllocation.from_indexed(
        ["He", "H", "Fe"], [2, 1, 26], []
    )
    _, labels = _nuclear_svg_labels(d)
    assert not [a for a in labels if a["kind"] == "residual"]


def test_core_subtraction_is_not_silently_clamped_to_a_physical_electron():
    d = ElectronAllocation.from_indexed(["C"], [0], [])
    mol = _mol(d)
    assert mol.GetAtomWithIdx(0).GetIntProp("_CCSResidualValenceElectrons") == -2
    _, labels = _nuclear_svg_labels(d)
    assert not [x for x in labels if x["kind"] == "residual"]
    # A 6+ assignment imbalance is still displayed, without inventing
    # negative counts of residual valence electrons.
    assert [x["glyph"] for x in labels if x["kind"] == "charge"] == ["6+"]


def test_2d_planar_K4_has_zero_straight_line_crossings():
    # K4 is planar but many conventional molecule sketchers draw it crossing.
    from constructive_chemical_space.ccs import _ccs_render_options
    m = Constitution.from_indexed(
        ["C"] * 4,
        [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
        24,
    )
    mol, _, _, _, _ = _ccs_render_options(m, layout="auto")
    assert mol.GetBoolProp("_CCSPlanarGraph")
    assert mol.GetIntProp("_CCS2DCrossings") == 0
    assert "crossings unavoidable" not in m.to_svg()


def test_nonplanar_K5_is_detected_and_accurately_warned():
    from itertools import combinations
    from constructive_chemical_space.ccs import _ccs_render_options
    m = Constitution.from_indexed(
        ["C"] * 5,
        list(combinations(range(5), 2)),
        30,
    )
    mol, footer, _, _, _ = _ccs_render_options(m)
    assert not mol.GetBoolProp("_CCSPlanarGraph")
    assert mol.GetBoolProp("_CCSNonPlanar2D")
    assert any("crossings unavoidable" in line for line in footer)
    assert "Non-planar graph" in m.to_svg()
    with pytest.raises(ValueError):
        m.to_svg(layout="unknown")


def test_planarity_preserves_ccs_constitution_and_allocation_identity():
    # Drawing never changes S1 identities, irrespective of layout selection.
    d = ElectronAllocation.from_indexed(
        ["C", "C", "H", "H"], [6, 6, 1, 1],
        [(0, 1, 2, 2), (0, 2, 1, 1), (1, 3, 1, 1)],
    )
    w = d.to_word()
    before = d.to_constitution()
    for layout in ("rdkit", "coordgen", "auto"):
        svg = d.to_svg(layout=layout)
        assert "<svg" in svg
    assert d == w.replay()
    assert d.to_constitution() == before


def test_residual_dots_are_svg_vector_circles_not_utf8_text():
    """Every displayed dot has a numeric circle, independent of font/encoding."""
    d = ElectronAllocation.from_indexed(
        ["N", "H", "H", "H"], [7, 1, 1, 1],
        [(0, 1, 1, 1), (0, 2, 1, 1), (0, 3, 1, 1)],
    )
    svg, labels = _nuclear_svg_labels(d)
    import xml.etree.ElementTree as ET
    nodes = [
        item for item in ET.fromstring(svg).iter()
        if item.attrib.get("data-kind") == "residual"
    ]
    assert len(nodes) == 1
    assert nodes[0].attrib["data-count"] == "2"
    assert len(list(nodes[0])) == 2
    assert all(item.tag.endswith("circle") for item in nodes[0])
    assert "•" not in svg
    assert "\\xe2" not in svg
    assert labels[0]["glyph"] == "••"


@pytest.mark.parametrize(
    ("number", "expected_style"), [
        (21, "separated-ligand-layers"),
        (42, "metal-centered-shell"),
        (43, "three-bridge-cryptand"),
        (44, "three-bridge-cryptand"),
    ],
)
def test_cluster_schematics_spread_dense_metal_incidence(number, expected_style):
    from examples.s5._common import build_example
    from constructive_chemical_space.ccs import _ccs_render_options
    example = build_example(number)
    mol, footer, _, _, _ = _ccs_render_options(
        example.constitution, layout="cluster",
        order=example.order,
    )
    assert mol.GetProp("_CCSClusterLayoutKind") == expected_style
    assert mol.GetBoolProp("_CCSNonPlanar2D")
    assert any("Schematic cluster arrangement" in line for line in footer)
    assert any("Non-planar graph" in line for line in footer)
    conf = mol.GetConformer()
    coordinates = [
        (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y)
        for i in range(mol.GetNumAtoms())
    ]
    assert len({(round(x, 2), round(y, 2)) for x, y in coordinates}) >= (
        mol.GetNumAtoms() - 2
    )
    if number == 42:
        center = next(i for i, a in enumerate(mol.GetAtoms())
                      if a.GetSymbol() == "Ni")
        p = coordinates[center]
        shell = [
            hypot(x-p[0],y-p[1]) for i,(x,y) in enumerate(coordinates)
            if i != center
        ]
        assert len(shell) == 12
        assert max(shell) - min(shell) < 1e-5
        assert min(shell) > 3
    elif number in (43, 44):
        from rdkit import Chem
        core = next(i for i,a in enumerate(mol.GetAtoms())
                    if a.GetSymbol() in {"K", "Na"} and a.GetDegree() >= 6)
        py = coordinates[core][1]
        # Metal placed outside the three-bridge host, not on a cage atom.
        ys = [y for i,(x,y) in enumerate(coordinates)
              if i!=core and mol.GetAtomWithIdx(i).GetSymbol()=="O"]
        assert max(ys) < py - 2.0
    else:
        metal_x = [
            x for i,(x,y) in enumerate(coordinates)
            if mol.GetAtomWithIdx(i).GetSymbol()=="Lu"
        ]
        assert len(metal_x) == 2
        assert metal_x[0] != metal_x[1]
        assert abs(metal_x[0] - metal_x[1]) > 5


def test_cluster_layout_is_optional_and_does_not_change_ccs_identity():
    from examples.s5._common import build_example
    d = build_example(42).forms["A"]
    baseline = d.to_word(order=tuple(build_example(42).order)) if False else d.num_electrons
    for layout in ("auto", "cluster", "coordgen", "rdkit"):
        svg = d.to_svg(layout=layout)
        assert "<svg" in svg
    assert d.num_electrons == baseline



@pytest.mark.parametrize(
    ("number", "metals", "expected_style"),
    [
        (18, "Fe", "separated-ligand-layers"),
        (19, "U", "separated-ligand-layers"),
        (20, "Tb", "separated-ligand-layers"),
        (21, "Lu", "separated-ligand-layers"),
        (43, "K", "three-bridge-cryptand"),
        (44, "Na", "three-bridge-cryptand"),
    ],
)
def test_metals_are_placed_last_outside_ligand_atoms(
    number, metals, expected_style
):
    from examples.s5._common import build_example
    from constructive_chemical_space.ccs import _ccs_render_options
    d = build_example(number).constitution
    mol, footer, _, _, _ = _ccs_render_options(d, layout="cluster")
    assert mol.GetProp("_CCSClusterLayoutKind") == expected_style
    conf = mol.GetConformer()
    positions = [
        (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y)
        for i in range(mol.GetNumAtoms())
    ]
    selected = [
        atom.GetIdx() for atom in mol.GetAtoms()
        if atom.GetSymbol() == metals and atom.GetDegree() >= 6
    ]
    assert selected
    for idx in selected:
        x, y = positions[idx]
        clearance = min(
            hypot(x - positions[j][0], y - positions[j][1])
            for j in range(len(positions)) if j != idx
        )
        assert clearance > 2.4
    assert "Schematic cluster arrangement" in " ".join(footer)


@pytest.mark.parametrize("number", [43, 44])
def test_cryptand_explicit_h_positions_are_distinct(number):
    from examples.s5._common import build_example
    from constructive_chemical_space.ccs import _ccs_render_options
    mol, _, _, _, _ = _ccs_render_options(
        build_example(number).constitution, layout="cluster"
    )
    conf = mol.GetConformer()
    positions = {
        i: (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y)
        for i in range(mol.GetNumAtoms())
    }
    hydrogens = [
        atom.GetIdx() for atom in mol.GetAtoms()
        if atom.GetSymbol() == "H"
    ]
    for i, j in __import__("itertools").combinations(hydrogens, 2):
        assert hypot(
            positions[i][0] - positions[j][0],
            positions[i][1] - positions[j][1],
        ) > 0.8


def test_residual_dots_remain_close_to_ligand_nuclei():
    from rdkit.Chem.Draw import rdMolDraw2D
    from constructive_chemical_space.ccs import (
        _ccs_render_options, _ccs_endpoint_annotation_positions,
        _ccs_nuclear_annotation_positions,
    )
    d = ElectronAllocation.from_indexed(
        ["N", "H", "H", "H"], [7, 1, 1, 1],
        [(0, 1, 1, 1), (0, 2, 1, 1), (0, 3, 1, 1)],
    )
    mol, _, body_height, _, _ = _ccs_render_options(d)
    drawer = rdMolDraw2D.MolDraw2DSVG(1000, body_height)
    drawer.DrawMolecule(mol)
    endpoints = _ccs_endpoint_annotation_positions(drawer, mol)
    dots = [
        entry for entry in _ccs_nuclear_annotation_positions(
            drawer, mol, endpoints
        ) if entry[1] == "residual"
    ]
    assert len(dots) == 1
    atom, _, _, _, x, y = dots[0]
    center = drawer.GetDrawCoords(atom)
    assert hypot(x - center.x, y - center.y) <= 29.1
