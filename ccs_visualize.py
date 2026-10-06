"""
CCS visualization helpers.

This module is deliberately downstream of the CCS core.  It never changes
CCS identity or replay semantics.  RDKit bond glyphs, atom charges, temporary
display indices, and layout are visualization aids only.

Default depiction rules
------------------------
Constitution:
    - draw only the connection graph;
    - do not infer bond orders or local charges;
    - show N_e in a footer.

Electron-allocation form:
    - q_v = Z(X(v)) - n_v is shown explicitly as the atom formal charge;
    - (1,1), (2,2), (3,3) are drawn as single/double/triple;
    - (2,0)/(0,2) are drawn as dative donor -> acceptor;
    - every other allocation is kept as a graph connection and labeled
      literally, e.g. "(0,0)", "1e (1,0)", "(3,2)", "(5,5)", "(0,6)";
    - n_empty is shown in the footer when nonzero.

Word / PracticalWord:
    - replay first, then draw the resulting electron-allocation form;
    - construction addresses are reused as temporary display indices when
      show_indices=True.

No aromaticity, oxidation state, spin, bond order, or valence is inferred by
CCS.  Conventional RDKit glyphs are used only for the small set of endpoint
allocations listed above.
"""

from __future__ import annotations

from html import escape as _xml_escape
from pathlib import Path
import re
from typing import Any, Iterable, Sequence, TYPE_CHECKING

if TYPE_CHECKING:
    from ccs_s1_s4 import (
        AllocatedConnection,
        Constitution,
        ElectronAllocation,
        Nucleus,
        PracticalWord,
        Word,
    )
else:
    AllocatedConnection = Constitution = ElectronAllocation = Nucleus = PracticalWord = Word = Any


class CCSVisualizationError(RuntimeError):
    """Raised when a CCS object cannot be rendered."""


def _load_rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError as exc:
        raise CCSVisualizationError(
            "RDKit is required for CCS visualization but is not installed."
        ) from exc
    return Chem, rdDepictor, rdMolDraw2D


def _ccs_kind(obj: Any) -> str:
    """Classify CCS objects by their public protocol, not module identity.

    This intentionally accepts objects from both ``ccs_s1_s4`` and the
    compatibility module ``ccs_s1_s4_rdkit``.
    """
    if hasattr(obj, "replay_state") and hasattr(obj, "tokens"):
        return "word"
    if hasattr(obj, "expand") and hasattr(obj, "items"):
        return "practical-word"
    if (
        hasattr(obj, "assigned_electrons")
        and hasattr(obj, "allocated_connections")
        and hasattr(obj, "num_unassigned_electron")
    ):
        return "allocation"
    if (
        hasattr(obj, "connections")
        and hasattr(obj, "num_electrons")
        and hasattr(obj, "total_charge")
        and hasattr(obj, "nuclei")
    ):
        return "constitution"
    raise TypeError(
        "expected a CCS Constitution, ElectronAllocation, Word, or PracticalWord"
    )


def _validate_order(
    nuclei: Iterable[Nucleus],
    order: Sequence[Nucleus] | None,
) -> tuple[Nucleus, ...]:
    nuclei = tuple(nuclei)
    nuclei_set = set(nuclei)
    if order is not None:
        order = tuple(order)
        if len(order) != len(nuclei) or set(order) != nuclei_set:
            raise ValueError("order must enumerate every nucleus exactly once")
        return order

    # Direct Constitution / ElectronAllocation objects are intrinsically
    # index-independent.  This order is only a display enumeration.  Ties
    # between automorphism-equivalent nuclei are intentionally not promoted
    # to chemical identity.
    def key(v: Nucleus):
        return (
            v.identity.atomic_number,
            v.identity.isotope or 0,
            v.identity.label,
            id(v),
        )

    return tuple(sorted(nuclei, key=key))


def _resolve_source(
    obj: Any,
    *,
    view: str,
    order: Sequence[Any] | None,
):
    if view not in {"auto", "constitution", "allocation"}:
        raise ValueError("view must be 'auto', 'constitution', or 'allocation'")

    kind = _ccs_kind(obj)

    if kind == "word":
        state = obj.replay_state()
        form = state.form
        source_order = state.address_table
        if order is not None:
            source_order = _validate_order(form.nuclei, order)
        if view == "constitution":
            return form.to_constitution(), source_order, "constitution"
        return form, source_order, "allocation"

    if kind == "practical-word":
        state = obj.expand().replay_state()
        form = state.form
        source_order = state.address_table
        if order is not None:
            source_order = _validate_order(form.nuclei, order)
        if view == "constitution":
            return form.to_constitution(), source_order, "constitution"
        return form, source_order, "allocation"

    if kind == "allocation":
        source_order = _validate_order(obj.nuclei, order)
        if view == "constitution":
            return obj.to_constitution(), source_order, "constitution"
        return obj, source_order, "allocation"

    if kind == "constitution":
        source_order = _validate_order(obj.nuclei, order)
        if view == "allocation":
            raise TypeError(
                "a Constitution has no electron-allocation data; "
                "use view='constitution' or provide an ElectronAllocation"
            )
        return obj, source_order, "constitution"

    raise AssertionError("unreachable CCS object kind")


def _allocation_style(Chem: Any, edge: AllocatedConnection):
    """Return (bond_type, donor, acceptor, conventional).

    Only chemically familiar display patterns are given conventional glyphs.
    Endpoint allocations remain the authoritative CCS data.
    """
    ai, aj = edge.allocation1, edge.allocation2

    if (ai, aj) == (1, 1):
        return Chem.BondType.SINGLE, None, None, True
    if (ai, aj) == (2, 2):
        return Chem.BondType.DOUBLE, None, None, True
    if (ai, aj) == (3, 3):
        return Chem.BondType.TRIPLE, None, None, True
    if (ai, aj) == (2, 0):
        return Chem.BondType.DATIVE, edge.endpoint1, edge.endpoint2, True
    if (ai, aj) == (0, 2):
        return Chem.BondType.DATIVE, edge.endpoint2, edge.endpoint1, True

    return Chem.BondType.ZERO, None, None, False


def _literal_allocation_label(edge: AllocatedConnection) -> str:
    ai, aj = edge.allocation1, edge.allocation2
    total = ai + aj
    if total in {0, 1}:
        return f"{total}e ({ai},{aj})"
    return f"({ai},{aj})"


def _set_atom_note(
    atom,
    nucleus: Nucleus,
    display_index: int,
    *,
    form: ElectronAllocation | None,
    show_indices: bool,
    diagnostics: bool,
) -> None:
    parts: list[str] = []

    if show_indices:
        parts.append(str(display_index))

    if diagnostics and form is not None:
        n = form.assigned_electrons(nucleus)
        r = form.residual_electrons(nucleus)
        q = form.formal_charge(nucleus)
        parts.append(f"n={n} r={r} q={q:+d}")

    if parts:
        atom.SetProp("atomNote", "  ".join(parts))


def _build_constitution_mol(
    constitution: Constitution,
    order: Sequence[Nucleus],
    *,
    show_hydrogens: bool,
    show_indices: bool,
):
    Chem, _, _ = _load_rdkit()

    visible = tuple(
        v for v in order
        if show_hydrogens or v.identity.symbol != "H"
    )
    visible_set = set(visible)

    rw = Chem.RWMol()
    index: dict[Nucleus, int] = {}
    display_index = {v: i for i, v in enumerate(order)}

    for nucleus in visible:
        atom = Chem.Atom(nucleus.identity.symbol)
        atom.SetNoImplicit(True)
        if nucleus.identity.isotope is not None:
            atom.SetIsotope(nucleus.identity.isotope)
        _set_atom_note(
            atom,
            nucleus,
            display_index[nucleus],
            form=None,
            show_indices=show_indices,
            diagnostics=False,
        )
        atom.SetProp("_CCSNuclearIdentity", nucleus.identity.label)
        atom.SetIntProp("_CCSDisplayIndex", display_index[nucleus])
        index[nucleus] = rw.AddAtom(atom)

    for edge in constitution.connections:
        if edge.endpoint1 not in visible_set or edge.endpoint2 not in visible_set:
            continue
        # Solid line = CCS direct connection only; NOT a bond-order inference.
        rw.AddBond(
            index[edge.endpoint1],
            index[edge.endpoint2],
            Chem.BondType.SINGLE,
        )
        bond = rw.GetBondBetweenAtoms(index[edge.endpoint1], index[edge.endpoint2])
        if bond is not None:
            bond.SetProp("_CCSDepictionRole", "constitution-connection")

    mol = rw.GetMol()
    mol.SetIntProp("_CCSTotalElectrons", constitution.num_electrons)
    mol.SetIntProp("_CCSTotalCharge", constitution.total_charge)
    mol.SetProp(
        "_CCSProjectionNote",
        "Solid lines in constitution view denote direct connections only; "
        "they are not CCS bond orders.",
    )
    return mol


def _build_allocation_mol(
    form: ElectronAllocation,
    order: Sequence[Nucleus],
    *,
    show_hydrogens: bool,
    show_indices: bool,
    diagnostics: bool,
    label_all_allocations: bool,
):
    Chem, _, _ = _load_rdkit()

    visible = tuple(
        v for v in order
        if show_hydrogens or v.identity.symbol != "H"
    )
    visible_set = set(visible)

    rw = Chem.RWMol()
    index: dict[Nucleus, int] = {}
    display_index = {v: i for i, v in enumerate(order)}

    for nucleus in visible:
        atom = Chem.Atom(nucleus.identity.symbol)
        atom.SetNoImplicit(True)
        if nucleus.identity.isotope is not None:
            atom.SetIsotope(nucleus.identity.isotope)

        # Explicitly show the CCS-derived local charge q_v = Z - n.
        atom.SetFormalCharge(form.formal_charge(nucleus))

        _set_atom_note(
            atom,
            nucleus,
            display_index[nucleus],
            form=form,
            show_indices=show_indices,
            diagnostics=diagnostics,
        )

        atom.SetIntProp("_CCSAssignedElectrons", form.assigned_electrons(nucleus))
        atom.SetIntProp("_CCSResidualElectrons", form.residual_electrons(nucleus))
        atom.SetIntProp("_CCSAuxCharge", form.formal_charge(nucleus))
        atom.SetIntProp("_CCSDisplayIndex", display_index[nucleus])
        atom.SetProp("_CCSNuclearIdentity", nucleus.identity.label)

        index[nucleus] = rw.AddAtom(atom)

    for edge in form.allocated_connections:
        if edge.endpoint1 not in visible_set or edge.endpoint2 not in visible_set:
            continue

        bond_type, donor, acceptor, conventional = _allocation_style(Chem, edge)

        if bond_type == Chem.BondType.DATIVE:
            assert donor is not None and acceptor is not None
            begin_nucleus, end_nucleus = donor, acceptor
        else:
            begin_nucleus, end_nucleus = edge.endpoint1, edge.endpoint2

        begin = index[begin_nucleus]
        end = index[end_nucleus]
        rw.AddBond(begin, end, bond_type)
        bond = rw.GetBondBetweenAtoms(begin, end)
        if bond is None:
            continue

        # Store allocation in the actual drawing direction.
        bond.SetIntProp(
            "_CCSBeginAllocation",
            edge.allocation_for(begin_nucleus),
        )
        bond.SetIntProp(
            "_CCSEndAllocation",
            edge.allocation_for(end_nucleus),
        )
        bond.SetProp(
            "_CCSAllocation",
            f"({edge.allocation1},{edge.allocation2})",
        )

        if label_all_allocations or not conventional:
            bond.SetProp("bondNote", _literal_allocation_label(edge))

        if conventional:
            bond.SetProp("_CCSDepictionRole", "conventional-glyph")
        else:
            bond.SetProp("_CCSDepictionRole", "literal-allocation")

    mol = rw.GetMol()
    mol.SetIntProp("_CCSTotalElectrons", form.num_electrons)
    mol.SetIntProp("_CCSUnassignedElectrons", form.num_unassigned_electron)
    mol.SetIntProp("_CCSTotalCharge", form.to_constitution().total_charge)
    mol.SetProp(
        "_CCSProjectionNote",
        "Bond glyphs and formal charges are depiction aids. "
        "Endpoint allocations remain the CCS data.",
    )
    return mol


def to_depiction_mol(
    obj: Constitution | ElectronAllocation | Word | PracticalWord,
    *,
    view: str = "auto",
    order: Sequence[Nucleus] | None = None,
    show_hydrogens: bool = True,
    show_indices: bool = False,
    diagnostics: bool = False,
    label_all_allocations: bool = False,
):
    """Build an unsanitized RDKit Mol for CCS drawing.

    For Word and PracticalWord inputs, replay is performed first and the
    construction address order is preserved for temporary display indices.
    """
    source, source_order, resolved_view = _resolve_source(
        obj,
        view=view,
        order=order,
    )

    if resolved_view == "constitution":
        return _build_constitution_mol(
            source,
            source_order,
            show_hydrogens=show_hydrogens,
            show_indices=show_indices,
        )

    return _build_allocation_mol(
        source,
        source_order,
        show_hydrogens=show_hydrogens,
        show_indices=show_indices,
        diagnostics=diagnostics,
        label_all_allocations=label_all_allocations,
    )


def _metadata_lines(
    source: Any,
    *,
    source_kind: str,
    show_total_charge: bool,
) -> list[str]:
    if source_kind == "constitution":
        lines = [f"N_e = {source.num_electrons}"]
        if show_total_charge:
            lines.append(f"Q = {source.total_charge:+d}")
        return lines

    lines = [f"N_e = {source.num_electrons}"]
    if show_total_charge:
        lines.append(f"Q = {source.to_constitution().total_charge:+d}")
    if source.num_unassigned_electron:
        lines.append(f"n_empty = {source.num_unassigned_electron}")
    return lines


def _compute_2d(mol, bond_length_scale: float) -> None:
    if bond_length_scale <= 0:
        raise ValueError("bond_length_scale must be positive")

    _, rdDepictor, _ = _load_rdkit()
    try:
        mol.UpdatePropertyCache(strict=False)
    except Exception:
        pass

    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)

    if bond_length_scale == 1:
        return

    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(
            i,
            (p.x * bond_length_scale, p.y * bond_length_scale, 0.0),
        )


def _drawing_options(drawer, *, annotation_scale: float):
    opts = drawer.drawOptions()
    opts.annotationFontScale = annotation_scale
    opts.bondLineWidth = 2.2
    opts.additionalAtomLabelPadding = 0.12
    opts.padding = 0.05
    return opts


def _svg_with_footer(
    mol,
    *,
    width: int,
    height: int,
    metadata: Sequence[str],
    annotation_scale: float,
    footer_font_size: int,
) -> str:
    _, _, rdMolDraw2D = _load_rdkit()

    footer_height = max(44, 10 + len(metadata) * (footer_font_size + 4))
    mol_height = max(120, height - footer_height)

    drawer = rdMolDraw2D.MolDraw2DSVG(width, mol_height)
    _drawing_options(drawer, annotation_scale=annotation_scale)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    # Expand the RDKit SVG canvas downward and add a right-aligned footer.
    svg = re.sub(
        rf"height=['\"]{mol_height}px['\"]",
        f"height='{height}px'",
        svg,
        count=1,
    )
    svg = re.sub(
        rf"viewBox=['\"]0 0 {width} {mol_height}['\"]",
        f"viewBox='0 0 {width} {height}'",
        svg,
        count=1,
    )

    y = mol_height + footer_font_size + 4
    footer = [
        f"<rect x='0' y='{mol_height}' width='{width}' height='{footer_height}' "
        "fill='white'/>"
    ]
    for line in metadata:
        footer.append(
            f"<text x='{width - 16}' y='{y}' text-anchor='end' "
            f"font-family='sans-serif' font-size='{footer_font_size}px'>"
            f"{_xml_escape(line)}</text>"
        )
        y += footer_font_size + 4

    return svg.replace("</svg>", "\n" + "\n".join(footer) + "\n</svg>")


def _png_with_footer(
    mol,
    *,
    width: int,
    height: int,
    metadata: Sequence[str],
    annotation_scale: float,
    footer_font_size: int,
) -> bytes:
    _, _, rdMolDraw2D = _load_rdkit()

    footer_height = max(44, 10 + len(metadata) * (footer_font_size + 4))
    mol_height = max(120, height - footer_height)

    drawer = rdMolDraw2D.MolDraw2DCairo(width, mol_height)
    _drawing_options(drawer, annotation_scale=annotation_scale)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    raw = drawer.GetDrawingText()

    try:
        from io import BytesIO
        from PIL import Image, ImageDraw, ImageFont
    except ImportError:
        # Pillow is optional.  If absent, return the RDKit PNG without footer.
        return raw

    base = Image.open(BytesIO(raw)).convert("RGBA")
    canvas = Image.new("RGBA", (width, height), "white")
    canvas.paste(base, (0, 0))

    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype("DejaVuSans.ttf", footer_font_size)
    except Exception:
        font = ImageFont.load_default()

    y = mol_height + 4
    for line in metadata:
        bbox = draw.textbbox((0, 0), line, font=font)
        text_width = bbox[2] - bbox[0]
        draw.text(
            (width - 16 - text_width, y),
            line,
            font=font,
            fill="black",
        )
        y += footer_font_size + 4

    out = BytesIO()
    canvas.convert("RGB").save(out, format="PNG")
    return out.getvalue()


def render(
    obj: Constitution | ElectronAllocation | Word | PracticalWord,
    *,
    view: str = "auto",
    fmt: str = "svg",
    order: Sequence[Nucleus] | None = None,
    show_hydrogens: bool = True,
    show_indices: bool = False,
    diagnostics: bool = False,
    label_all_allocations: bool = False,
    show_total_charge: bool = False,
    width: int = 1200,
    height: int = 800,
    bond_length_scale: float = 1.25,
    annotation_scale: float = 0.8,
    footer_font_size: int = 20,
):
    """Render a CCS object as SVG text or PNG bytes."""
    if width <= 0 or height <= 0:
        raise ValueError("width and height must be positive")
    if fmt.lower() not in {"svg", "png"}:
        raise ValueError("fmt must be 'svg' or 'png'")

    source, source_order, resolved_view = _resolve_source(
        obj,
        view=view,
        order=order,
    )

    if resolved_view == "constitution":
        mol = _build_constitution_mol(
            source,
            source_order,
            show_hydrogens=show_hydrogens,
            show_indices=show_indices,
        )
    else:
        mol = _build_allocation_mol(
            source,
            source_order,
            show_hydrogens=show_hydrogens,
            show_indices=show_indices,
            diagnostics=diagnostics,
            label_all_allocations=label_all_allocations,
        )

    _compute_2d(mol, bond_length_scale)
    metadata = _metadata_lines(
        source,
        source_kind=resolved_view,
        show_total_charge=show_total_charge,
    )

    if fmt.lower() == "svg":
        return _svg_with_footer(
            mol,
            width=width,
            height=height,
            metadata=metadata,
            annotation_scale=annotation_scale,
            footer_font_size=footer_font_size,
        )

    return _png_with_footer(
        mol,
        width=width,
        height=height,
        metadata=metadata,
        annotation_scale=annotation_scale,
        footer_font_size=footer_font_size,
    )


def draw(
    obj: Constitution | ElectronAllocation | Word | PracticalWord,
    **kwargs,
):
    """Alias for :func:`render`."""
    return render(obj, **kwargs)


def draw_constitution(
    constitution: Constitution,
    **kwargs,
):
    """Render a constitution as a graph-only CCS view."""
    kwargs["view"] = "constitution"
    return render(constitution, **kwargs)


def draw_allocation(
    form: ElectronAllocation,
    **kwargs,
):
    """Render an electron-allocation form."""
    kwargs["view"] = "allocation"
    return render(form, **kwargs)


def draw_word(
    word: Word | PracticalWord,
    **kwargs,
):
    """Replay and render a construction word."""
    kwargs.setdefault("view", "allocation")
    return render(word, **kwargs)


def save_drawing(
    obj: Constitution | ElectronAllocation | Word | PracticalWord,
    path: str | Path,
    **kwargs,
) -> Path:
    """Render to a file; extension selects SVG or PNG unless fmt is supplied."""
    path = Path(path)
    fmt = kwargs.pop("fmt", path.suffix.lower().lstrip("."))
    if fmt not in {"svg", "png"}:
        raise ValueError("output path must end in .svg or .png, or pass fmt='svg'/'png'")

    data = render(obj, fmt=fmt, **kwargs)
    path.parent.mkdir(parents=True, exist_ok=True)

    if fmt == "svg":
        path.write_text(data, encoding="utf-8")
    else:
        path.write_bytes(data)

    return path


__all__ = [
    "CCSVisualizationError",
    "to_depiction_mol",
    "render",
    "draw",
    "draw_constitution",
    "draw_allocation",
    "draw_word",
    "save_drawing",
]