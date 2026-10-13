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


def test_2d_planar_K4_preserves_chemical_ring_layout():
    # K4 is planar but many conventional molecule sketchers draw it crossing.
    from constructive_chemical_space.ccs import _ccs_render_options
    m = Constitution.from_indexed(
        ["C"] * 4,
        [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
        24,
    )
    mol, _, _, _, _ = _ccs_render_options(m, layout="auto")
    assert mol.GetBoolProp("_CCSPlanarGraph")
    assert mol.GetProp("_CCSLayoutMethod") in {"rdkit", "coordgen"}
    assert mol.GetIntProp("_CCS2DCrossings") >= 0
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


@pytest.mark.parametrize("number", [18, 19, 20, 21, 42, 43, 44])
def test_dense_s5_layout_is_independent_of_element_category(number):
    from examples.s5._common import build_example
    from constructive_chemical_space.ccs import _ccs_render_options
    from networkx import Graph, check_planarity

    example = build_example(number)
    mol, footer, _, _, _ = _ccs_render_options(
        example.constitution, layout="auto", order=example.order
    )
    assert mol.GetNumAtoms() == len(example.order)
    assert mol.GetNumBonds() == len(example.constitution.connections)
    assert mol.GetProp("_CCSLayoutMethod") in {
        "rdkit", "coordgen", "coordgen-light-hubs", "spring"
    }
    assert not mol.HasProp("_CCSClusterLayoutKind")
    graph = Graph()
    graph.add_nodes_from(range(mol.GetNumAtoms()))
    graph.add_edges_from(
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()
    )
    planar = check_planarity(graph)[0]
    assert mol.GetBoolProp("_CCSPlanarGraph") == planar
    assert mol.GetIntProp("_CCS2DCrossings") >= 0
    if planar:
        assert not any("Non-planar graph" in line for line in footer)
    else:
        assert any("Non-planar graph" in line for line in footer)


def test_old_cluster_specific_layout_flag_is_not_supported():
    from examples.s5._common import build_example
    with pytest.raises(ValueError, match="layout"):
        build_example(18).constitution.to_svg(layout="cluster")


@pytest.mark.parametrize("number", [43, 44])
def test_generic_spacing_keeps_explicit_h_nuclei_distinct(number):
    from itertools import combinations
    from statistics import median
    from examples.s5._common import build_example
    from constructive_chemical_space.ccs import _ccs_render_options

    example = build_example(number)
    mol, _, _, _, _ = _ccs_render_options(
        example.constitution, layout="auto", order=example.order
    )
    conf = mol.GetConformer()
    positions = [
        (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y)
        for i in range(mol.GetNumAtoms())
    ]
    sizes = [
        hypot(positions[b.GetBeginAtomIdx()][0] -
              positions[b.GetEndAtomIdx()][0],
              positions[b.GetBeginAtomIdx()][1] -
              positions[b.GetEndAtomIdx()][1])
        for b in mol.GetBonds()
    ]
    typical = median(sizes)
    hydrogens = [
        a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "H"
    ]
    assert all(
        hypot(positions[a][0]-positions[b][0],
              positions[a][1]-positions[b][1]) > 0.10 * typical
        for a, b in combinations(hydrogens, 2)
    )


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


def test_generic_layout_planarity_does_not_depend_on_element_identity():
    from itertools import combinations
    from constructive_chemical_space.ccs import _ccs_render_options
    edges = list(combinations(range(5), 2))
    for identities in (["C"] * 5, ["N"] * 5):
        obj = Constitution.from_indexed(identities, edges, 30)
        mol, _, _, _, _ = _ccs_render_options(obj)
        assert mol.GetNumBonds() == 10
        assert not mol.GetBoolProp("_CCSPlanarGraph")
        assert mol.GetIntProp("_CCS2DCrossings") > 0


def test_s5_ferrocene_cp_rings_are_pentagonal():
    """Both Cp rings remain legible despite allowed Fe–C line intersections."""
    from examples.s5._common import build_example
    from constructive_chemical_space.ccs import _ccs_render_options
    import networkx as nx
    ex = build_example(18)
    mol, _, _, _, _ = _ccs_render_options(
        ex.constitution, order=ex.order, layout="auto"
    )
    graph = nx.Graph()
    graph.add_nodes_from(range(mol.GetNumAtoms()))
    graph.add_edges_from(
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()
    )
    cycles = [r for r in nx.cycle_basis(
        graph.subgraph(v for v, d in graph.degree() if d <= 5)
    ) if len(r) == 5]
    assert len(cycles) == 2
    con = mol.GetConformer()
    for ring in cycles:
        p = [(con.GetAtomPosition(i).x, con.GetAtomPosition(i).y)
             for i in ring]
        lengths = [
            hypot(p[k][0]-p[(k+1)%5][0],p[k][1]-p[(k+1)%5][1])
            for k in range(5)
        ]
        mean = sum(lengths)/5
        variation = (sum((d-mean)**2 for d in lengths)/5)**0.5/mean
        assert variation < 0.18
    assert mol.GetProp("_CCSLayoutMethod") in {"coordgen", "coordgen-light-hubs", "rdkit"}


@pytest.mark.parametrize("number", [47, 48])
def test_s5_perfluorocubane_square_faces_do_not_collapse(number):
    """Perfluorocubane must retain recognizable nondegenerate 4-rings."""
    from statistics import median
    from examples.s5._common import build_example
    from constructive_chemical_space.ccs import _ccs_render_options
    import networkx as nx
    ex = build_example(number)
    mol, _, _, _, _ = _ccs_render_options(
        ex.constitution, order=ex.order, layout="auto"
    )
    con = mol.GetConformer()
    pts = [(con.GetAtomPosition(i).x, con.GetAtomPosition(i).y)
           for i in range(mol.GetNumAtoms())]
    g = nx.Graph()
    g.add_nodes_from(range(mol.GetNumAtoms()))
    g.add_edges_from(
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()
    )
    scale = median(
        hypot(pts[a][0]-pts[b][0],pts[a][1]-pts[b][1])
        for a, b in g.edges()
    )
    rings = [r for r in nx.cycle_basis(
        g.subgraph(v for v, d in g.degree() if d <= 5)
    ) if len(r) == 4]
    assert len(rings) >= 4
    for ring in rings:
        area = abs(sum(
            pts[ring[k]][0] * pts[ring[(k+1)%4]][1]
            - pts[ring[(k+1)%4]][0] * pts[ring[k]][1]
            for k in range(4)
        ))/2.0
        assert area/scale**2 > 0.15


def test_short_ring_quality_is_priority_over_no_crossings():
    from examples.s5._common import build_example
    from constructive_chemical_space.ccs import _ccs_render_options
    for number in (18,47):
        ex=build_example(number)
        mol,_,_,_,_=_ccs_render_options(
            ex.constitution, order=ex.order, layout="auto"
        )
        assert mol.GetProp("_CCSLayoutMethod") in {"rdkit","coordgen"}
        assert mol.GetDoubleProp("_CCSRingDistortion") < 90.0

@pytest.mark.parametrize(
    ("number", "expected_method"),
    [(18, "rdkit"), (47, "coordgen"), (48, "coordgen")],
)
def test_s5_auto_selects_generic_ring_hub_or_cage_projection(number, expected_method):
    """Compare actual S5 candidates, independently of allocation bond glyphs."""
    from examples.s5._common import build_example
    from constructive_chemical_space.ccs import _ccs_render_options

    example = build_example(number)
    reference_points = None
    for source in (example.constitution, *example.forms.values()):
        mol, _, _, _, _ = _ccs_render_options(
            source, order=example.order, layout="auto"
        )
        assert mol.GetProp("_CCSLayoutMethod") == expected_method
        scores = {
            method: float(value)
            for method, value in (
                entry.split("=")
                for entry in mol.GetProp("_CCSLayoutCandidateScores").split(";")
            )
        }
        competitor = "coordgen" if expected_method == "rdkit" else "rdkit"
        assert scores[expected_method] < scores[competitor]
        assert mol.GetDoubleProp("_CCSLayoutScore") == pytest.approx(
            scores[expected_method], rel=1e-8
        )

        # The nuclear skeleton is common to M and all D over M.
        conf = mol.GetConformer()
        points = [
            coordinate
            for i in range(mol.GetNumAtoms())
            for coordinate in (
                conf.GetAtomPosition(i).x,
                conf.GetAtomPosition(i).y,
            )
        ]
        if reference_points is None:
            reference_points = points
        else:
            assert points == pytest.approx(reference_points, abs=1e-6)

        if number == 18 and source is not example.constitution:
            # A temporary all-SINGLE coordinate model cannot mutate the
            # final rendered allocation's genuine line/arrow styles.
            styles = {bond.GetBondType() for bond in mol.GetBonds()}
            assert {Chem.BondType.SINGLE, Chem.BondType.DATIVE,
                    Chem.BondType.ZERO} <= styles

    for method in ("coordgen", "rdkit"):
        forced, _, _, _, _ = _ccs_render_options(
            example.constitution, order=example.order, layout=method
        )
        assert forced.GetProp("_CCSLayoutMethod") == method
        assert forced.GetDoubleProp("_CCSLayoutScore") == pytest.approx(
            scores[method], rel=1e-8
        )
