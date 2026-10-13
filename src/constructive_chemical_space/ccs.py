"""
Constructive Chemical Space (CCS)
Reference implementation of Supporting Information Sections S1-S4 for
Constructive Chemical Space (CCS).

Implemented layers
------------------
S1. Intrinsic, index-independent CCS constitutions and electron-allocation forms.
S2. Core construction tokens Nuc / Ele / Con and total left-to-right replay.
S3. Exact construction witnesses and the all-unassigned reference lift.
S4. Terminal End completion, explicit-H core expansion, finite macros,
    and optional search-policy filtering.

Important design rule
---------------------
Nuclear indices are NOT stored in Constitution or ElectronAllocation.
Python Nucleus objects are the abstract members of V. Construction-time integer
addresses live only in ExecutionState / construction words and are discarded
when replay returns an intrinsic ElectronAllocation.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import ceil, hypot, cos, sin, pi
from pathlib import Path
import re
from types import MappingProxyType
import warnings
from typing import Any, Callable, Iterable, Iterator, Mapping, Sequence


# ---------------------------------------------------------------------------
# Nuclear identities
# ---------------------------------------------------------------------------

_ELEMENT_SYMBOLS = (
    None,
    "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
)
_ATOMIC_NUMBER = {
    symbol: z for z, symbol in enumerate(_ELEMENT_SYMBOLS) if symbol is not None
}

_GROUP = {}
for g, symbols in {
    1: ("H", "Li", "Na", "K", "Rb", "Cs", "Fr"),
    2: ("Be", "Mg", "Ca", "Sr", "Ba", "Ra"),
    13: ("B", "Al", "Ga", "In", "Tl", "Nh"),
    14: ("C", "Si", "Ge", "Sn", "Pb", "Fl"),
    15: ("N", "P", "As", "Sb", "Bi", "Mc"),
    16: ("O", "S", "Se", "Te", "Po", "Lv"),
    17: ("F", "Cl", "Br", "I", "At", "Ts"),
    18: ("He", "Ne", "Ar", "Kr", "Xe", "Rn", "Og"),
}.items():
    for symbol in symbols:
        _GROUP[symbol] = g


@dataclass(frozen=True, slots=True)
class NuclearIdentity:
    """A nuclear-identity label X.

    Isotope-unspecified and isotope-resolved identities are distinct values.
    The isotope is part of X, not a separate electronic field.
    """

    symbol: str
    isotope: int | None = None

    def __post_init__(self) -> None:
        if self.symbol not in _ATOMIC_NUMBER:
            raise ValueError(f"unknown element symbol: {self.symbol!r}")
        if self.isotope is not None:
            if not _is_nonnegative_int(self.isotope) or self.isotope == 0:
                raise ValueError("isotope mass number must be a positive integer")

    @property
    def atomic_number(self) -> int:
        return _ATOMIC_NUMBER[self.symbol]

    @property
    def group(self) -> int | None:
        return _GROUP.get(self.symbol)

    @property
    def label(self) -> str:
        return self.symbol if self.isotope is None else f"{self.isotope}{self.symbol}"

    def __str__(self) -> str:
        return self.label

    @classmethod
    def coerce(cls, value: str | "NuclearIdentity") -> "NuclearIdentity":
        if isinstance(value, NuclearIdentity):
            return value
        if not isinstance(value, str):
            raise TypeError("nuclear identity must be a string or NuclearIdentity")

        text = value.strip()
        # Accept C, 13C, ^13C, [13C], C-13.
        if text.startswith("[") and text.endswith("]"):
            text = text[1:-1].strip()
        text = text.replace("^{", "").replace("}", "")
        if text.startswith("^"):
            text = text[1:]

        m = re.fullmatch(r"(\d+)?([A-Z][a-z]?)", text)
        if m:
            isotope = int(m.group(1)) if m.group(1) else None
            return cls(m.group(2), isotope)

        m = re.fullmatch(r"([A-Z][a-z]?)-(\d+)", text)
        if m:
            return cls(m.group(1), int(m.group(2)))

        raise ValueError(f"cannot parse nuclear identity: {value!r}")


def _is_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _is_nonnegative_int(value: object) -> bool:
    return _is_int(value) and value >= 0



# ---------------------------------------------------------------------------
# Optional RDKit projection
# ---------------------------------------------------------------------------

class RDKitConversionError(RuntimeError):
    """Raised when a CCS electron-allocation form cannot be projected to RDKit."""


def _load_rdkit():
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise RDKitConversionError(
            "RDKit is required for this operation but is not installed."
        ) from exc
    return Chem


def _rdkit_bond_type(Chem: Any, edge: "AllocatedConnection"):
    """Infer an RDKit bond representation from a CCS endpoint allocation.

    Returns
    -------
    (bond_type, donor_endpoint, acceptor_endpoint, inference_label)

    The mapping is deliberately conservative:
      (0,0)       -> UNSPECIFIED
      (1,1)       -> SINGLE
      (2,2)       -> DOUBLE
      (3,3)       -> TRIPLE
      (4,4)       -> QUADRUPLE
      (5,5)       -> QUINTUPLE
      (6,6)       -> HEXTUPLE
      (2,0)/(0,2) -> DATIVE, directed donor -> acceptor
      otherwise   -> UNSPECIFIED

    CCS endpoint allocations are not bond orders, so unfamiliar patterns are
    not forced into a conventional multiplicity.
    """
    ai, aj = edge.allocation1, edge.allocation2

    if ai == aj:
        if ai == 0:
            return Chem.BondType.UNSPECIFIED, None, None, "zero-allocation"
        bond_types = {
            1: Chem.BondType.SINGLE,
            2: Chem.BondType.DOUBLE,
            3: Chem.BondType.TRIPLE,
            4: Chem.BondType.QUADRUPLE,
            5: Chem.BondType.QUINTUPLE,
            6: Chem.BondType.HEXTUPLE,
        }
        if ai in bond_types:
            return bond_types[ai], None, None, f"symmetric-{ai}"
        return Chem.BondType.UNSPECIFIED, None, None, "unsupported-symmetric"

    if (ai, aj) == (2, 0):
        return Chem.BondType.DATIVE, edge.endpoint1, edge.endpoint2, "dative"
    if (ai, aj) == (0, 2):
        return Chem.BondType.DATIVE, edge.endpoint2, edge.endpoint1, "dative"

    return Chem.BondType.UNSPECIFIED, None, None, "asymmetric-unspecified"


def _validate_nucleus_order(
    form: "ElectronAllocation",
    order: Sequence["Nucleus"] | None,
) -> tuple["Nucleus", ...]:
    if order is None:
        return tuple(form.nuclei)
    order = tuple(order)
    if len(order) != len(form.nuclei) or set(order) != set(form.nuclei):
        raise ValueError("order must enumerate every nucleus exactly once")
    return order


def _build_rdkit_mol_once(
    form: "ElectronAllocation",
    *,
    order: Sequence["Nucleus"] | None,
    dative_charge_shift: bool,
    include_ccs_properties: bool,
):
    Chem = _load_rdkit()
    order = _validate_nucleus_order(form, order)

    rw = Chem.RWMol()
    index: dict[Nucleus, int] = {}

    # Add nuclei as explicit RDKit atoms.  Implicit H is disabled because
    # hydrogen is explicit in CCS.
    for nucleus in order:
        atom = Chem.Atom(nucleus.identity.symbol)
        atom.SetNoImplicit(True)

        if nucleus.identity.isotope is not None:
            atom.SetIsotope(nucleus.identity.isotope)

        # RDKit requires a toolkit formal-charge field.  We use CCS q_v only
        # as an initial projection heuristic; RDKit sanitization (including the
        # optional dative retry) may yield different formal charges.  CCS q_v
        # is an assignment imbalance, never a Lewis-charge assertion.
        q_v = form.assignment_imbalance(nucleus)
        atom.SetFormalCharge(q_v)

        if include_ccs_properties:
            atom.SetIntProp("_CCSAssignedElectrons", form.assigned_electrons(nucleus))
            atom.SetIntProp("_CCSResidualElectrons", form.residual_electrons(nucleus))
            atom.SetIntProp("_CCSAssignmentImbalance", q_v)
            atom.SetIntProp("_CCSAuxCharge", q_v)  # backward-compatible RDKit property
            atom.SetProp("_CCSNuclearIdentity", nucleus.identity.label)

        index[nucleus] = rw.AddAtom(atom)

    dative_pairs: list[tuple[Nucleus, Nucleus]] = []

    for edge in form.allocated_connections:
        bond_type, donor, acceptor, inference = _rdkit_bond_type(Chem, edge)

        if bond_type == Chem.BondType.DATIVE:
            assert donor is not None and acceptor is not None
            begin = index[donor]
            end = index[acceptor]
            dative_pairs.append((donor, acceptor))
        else:
            begin = index[edge.endpoint1]
            end = index[edge.endpoint2]

        rw.AddBond(begin, end, bond_type)
        bond = rw.GetBondBetweenAtoms(begin, end)

        if include_ccs_properties and bond is not None:
            # Store the allocations in the actual RDKit bond direction.
            begin_nucleus = donor if bond_type == Chem.BondType.DATIVE else edge.endpoint1
            end_nucleus = acceptor if bond_type == Chem.BondType.DATIVE else edge.endpoint2
            bond.SetIntProp(
                "_CCSBeginAllocation",
                edge.allocation_for(begin_nucleus),
            )
            bond.SetIntProp(
                "_CCSEndAllocation",
                edge.allocation_for(end_nucleus),
            )
            bond.SetProp("_CCSBondInference", inference)

    mol = rw.GetMol()

    # A donor allocation (2,0) can require the conventional donor+/acceptor-
    # bookkeeping for RDKit valence sanitization (e.g. NH3->BH3).  This is
    # tried only as a fallback after the direct q_v projection fails.
    if dative_charge_shift:
        for donor, acceptor in dative_pairs:
            donor_atom = mol.GetAtomWithIdx(index[donor])
            acceptor_atom = mol.GetAtomWithIdx(index[acceptor])
            donor_atom.SetFormalCharge(donor_atom.GetFormalCharge() + 1)
            acceptor_atom.SetFormalCharge(acceptor_atom.GetFormalCharge() - 1)

        if include_ccs_properties and dative_pairs:
            mol.SetProp("_CCSRDKitChargeAdjustment", "dative_retry")

    if include_ccs_properties:
        mol.SetIntProp("_CCSTotalElectrons", form.num_electrons)
        mol.SetIntProp("_CCSUnassignedElectrons", form.num_unassigned_electrons)
        mol.SetProp(
            "_CCSProjectionNote",
            "RDKit bond types/formal charges are a downstream best-effort projection; "
            "CCS endpoint allocations remain stored as properties.",
        )

    return mol


def _electron_allocation_to_rdkit_mol(
    form: "ElectronAllocation",
    *,
    sanitize: bool = True,
    allow_unassigned: bool = False,
    order: Sequence["Nucleus"] | None = None,
    include_ccs_properties: bool = True,
    strict: bool = False,
):
    """Best-effort projection of an ElectronAllocation to an RDKit Mol.

    If the form contains unassigned electrons, the default behavior is to
    return None because ordinary RDKit molecular graphs/SMILES have no
    index-independent place to store n_empty.  Set allow_unassigned=True to
    obtain a structural Mol carrying n_empty only as a molecule property.

    The initial RDKit formal-charge projection uses the CCS assignment
    imbalance q_v = Z - n as a heuristic, not a CCS definition of Lewis charge.
    If sanitization fails and a (2,0)/(0,2) dative allocation is present,
    a second attempt shifts donor/acceptor charges by +1/-1, respectively.
    _CCSAssignmentImbalance retains the unchanged CCS quantity.
    """
    if form.num_unassigned_electrons and not allow_unassigned:
        message = (
            "the CCS form contains unassigned electrons (n_empty > 0), which "
            "cannot be represented faithfully in an ordinary RDKit molecular graph"
        )
        if strict:
            raise RDKitConversionError(message)
        return None

    Chem = _load_rdkit()

    base = _build_rdkit_mol_once(
        form,
        order=order,
        dative_charge_shift=False,
        include_ccs_properties=include_ccs_properties,
    )

    if not sanitize:
        return base

    try:
        Chem.SanitizeMol(base)
        if include_ccs_properties:
            base.SetBoolProp("_CCSSanitized", True)
        return base
    except Exception as first_error:
        has_dative = any(
            (edge.allocation1, edge.allocation2) in {(2, 0), (0, 2)}
            for edge in form.allocated_connections
        )

        if has_dative:
            retry = _build_rdkit_mol_once(
                form,
                order=order,
                dative_charge_shift=True,
                include_ccs_properties=include_ccs_properties,
            )
            try:
                Chem.SanitizeMol(retry)
                if include_ccs_properties:
                    retry.SetBoolProp("_CCSSanitized", True)
                return retry
            except Exception as second_error:
                if strict:
                    raise RDKitConversionError(
                        "RDKit sanitization failed for both the direct and "
                        "dative-charge-adjusted projections"
                    ) from second_error

        if strict:
            raise RDKitConversionError(
                "RDKit sanitization failed for the CCS electron-allocation form"
            ) from first_error

        # Best-effort Mol: preserve the graph and CCS properties even when RDKit
        # cannot assign a chemically sanitized valence model.
        if include_ccs_properties:
            base.SetBoolProp("_CCSSanitized", False)
            base.SetProp("_CCSSanitizeError", str(first_error))
        return base


def _electron_allocation_to_smiles(
    form: "ElectronAllocation",
    *,
    canonical: bool = True,
    isomeric: bool = True,
    remove_hs: bool = False,
    sanitize: bool = True,
    allow_unassigned: bool = False,
    strict: bool = False,
) -> str | None:
    """Convert a CCS allocation form to a best-effort RDKit SMILES.

    This is not a lossless CCS serialization.  Endpoint-allocation properties,
    residual counts, and n_empty are not encoded in ordinary SMILES.
    """
    Chem = _load_rdkit()
    mol = _electron_allocation_to_rdkit_mol(
        form,
        sanitize=sanitize,
        allow_unassigned=allow_unassigned,
        include_ccs_properties=True,
        strict=strict,
    )
    if mol is None:
        return None

    if remove_hs:
        try:
            mol = Chem.RemoveHs(mol, sanitize=sanitize)
        except Exception as exc:
            if strict:
                raise RDKitConversionError("RDKit could not remove explicit H atoms") from exc

    try:
        return Chem.MolToSmiles(
            mol,
            canonical=canonical,
            isomericSmiles=isomeric,
        )
    except Exception as exc:
        if strict:
            raise RDKitConversionError("RDKit could not produce SMILES") from exc
        return None


# ---------------------------------------------------------------------------
# S1: intrinsic nuclear collections, constitutions, and allocation forms
# ---------------------------------------------------------------------------

@dataclass(frozen=True, eq=False, slots=True)
class Nucleus:
    """One abstract member v of the finite carrier V.

    Equality is object identity, not nuclear identity.  Two H nuclei are two
    distinct carrier elements even though both have X(v) == H.
    """

    identity: NuclearIdentity

    def __init__(self, identity: str | NuclearIdentity):
        object.__setattr__(self, "identity", NuclearIdentity.coerce(identity))

    # eq=False intentionally leaves object identity semantics.

    def __repr__(self) -> str:
        return f"Nucleus({self.identity.label})"


@dataclass(frozen=True, eq=False, slots=True)
class NuclearCollection:
    """C_X = (V, X), represented intrinsically as an unordered finite carrier."""

    nuclei: frozenset[Nucleus]

    def __init__(self, nuclei: Iterable[Nucleus] = ()):
        fs = frozenset(nuclei)
        if not all(isinstance(v, Nucleus) for v in fs):
            raise TypeError("NuclearCollection requires Nucleus objects")
        object.__setattr__(self, "nuclei", fs)

    def __len__(self) -> int:
        return len(self.nuclei)

    def __iter__(self) -> Iterator[Nucleus]:
        return iter(self.nuclei)

    def identity_multiset(self) -> Counter[NuclearIdentity]:
        return Counter(v.identity for v in self.nuclei)

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, NuclearCollection)
            and self.identity_multiset() == other.identity_multiset()
        )

    __hash__ = None


@dataclass(frozen=True, eq=False, slots=True)
class Connection:
    """One unordered pair {v,w} of distinct nuclei."""

    endpoint1: Nucleus
    endpoint2: Nucleus

    def __post_init__(self) -> None:
        if not isinstance(self.endpoint1, Nucleus) or not isinstance(self.endpoint2, Nucleus):
            raise TypeError("connection endpoints must be Nucleus objects")
        if self.endpoint1 is self.endpoint2:
            raise ValueError("self-connections are not allowed")

    @property
    def endpoints(self) -> frozenset[Nucleus]:
        return frozenset((self.endpoint1, self.endpoint2))

    def contains(self, nucleus: Nucleus) -> bool:
        return nucleus is self.endpoint1 or nucleus is self.endpoint2

    def other(self, nucleus: Nucleus) -> Nucleus:
        if nucleus is self.endpoint1:
            return self.endpoint2
        if nucleus is self.endpoint2:
            return self.endpoint1
        raise KeyError("nucleus is not an endpoint of this connection")

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Connection) and self.endpoints == other.endpoints

    def __hash__(self) -> int:
        return hash(self.endpoints)

    def __repr__(self) -> str:
        return f"Connection({self.endpoint1!r}, {self.endpoint2!r})"


@dataclass(frozen=True, eq=False, slots=True)
class AllocatedConnection:
    """A connection together with its two endpoint electron allocations."""

    endpoint1: Nucleus
    endpoint2: Nucleus
    allocation1: int
    allocation2: int

    def __post_init__(self) -> None:
        if self.endpoint1 is self.endpoint2:
            raise ValueError("self-connections are not allowed")
        if not _is_nonnegative_int(self.allocation1):
            raise ValueError("endpoint allocation must be a nonnegative integer")
        if not _is_nonnegative_int(self.allocation2):
            raise ValueError("endpoint allocation must be a nonnegative integer")

    @property
    def connection(self) -> Connection:
        return Connection(self.endpoint1, self.endpoint2)

    def allocation_for(self, nucleus: Nucleus) -> int:
        if nucleus is self.endpoint1:
            return self.allocation1
        if nucleus is self.endpoint2:
            return self.allocation2
        raise KeyError("nucleus is not an endpoint of this connection")

    def other(self, nucleus: Nucleus) -> Nucleus:
        return self.connection.other(nucleus)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, AllocatedConnection):
            return False
        return (
            (
                self.endpoint1 is other.endpoint1
                and self.endpoint2 is other.endpoint2
                and self.allocation1 == other.allocation1
                and self.allocation2 == other.allocation2
            )
            or (
                self.endpoint1 is other.endpoint2
                and self.endpoint2 is other.endpoint1
                and self.allocation1 == other.allocation2
                and self.allocation2 == other.allocation1
            )
        )

    def __hash__(self) -> int:
        return hash(
            frozenset(
                (
                    (self.endpoint1, self.allocation1),
                    (self.endpoint2, self.allocation2),
                )
            )
        )


@dataclass(frozen=True, eq=False, init=False, slots=True)
class Constitution:
    """M = (C_X, S, N_e), intrinsically index-independent.

    Immutable after construction, so the S1.3 invariants cannot be bypassed
    through ordinary attribute reassignment.
    """

    collection: NuclearCollection
    connections: frozenset[Connection]
    num_electrons: int

    def __init__(
        self,
        collection: NuclearCollection | Iterable[Nucleus],
        connections: Iterable[Connection] = (),
        num_electrons: int = 0,
    ) -> None:
        normalized = (
            collection if isinstance(collection, NuclearCollection)
            else NuclearCollection(collection)
        )
        if not _is_nonnegative_int(num_electrons):
            raise ValueError("num_electrons must be a nonnegative integer")

        conns = frozenset(connections)
        for edge in conns:
            if not isinstance(edge, Connection):
                raise TypeError("connections must contain Connection objects")
            if edge.endpoint1 not in normalized.nuclei or edge.endpoint2 not in normalized.nuclei:
                raise ValueError("every connection endpoint must belong to the nuclear collection")

        object.__setattr__(self, "collection", normalized)
        object.__setattr__(self, "connections", conns)
        object.__setattr__(self, "num_electrons", num_electrons)

    @property
    def nuclei(self) -> frozenset[Nucleus]:
        return self.collection.nuclei

    @property
    def num_total_electron(self) -> int:
        return self.num_electrons

    @property
    def total_charge(self) -> int:
        return sum(v.identity.atomic_number for v in self.nuclei) - self.num_electrons

    def degree(self, nucleus: Nucleus) -> int:
        self._require_nucleus(nucleus)
        return sum(edge.contains(nucleus) for edge in self.connections)

    def neighbors(self, nucleus: Nucleus) -> frozenset[Nucleus]:
        self._require_nucleus(nucleus)
        return frozenset(edge.other(nucleus) for edge in self.connections if edge.contains(nucleus))

    def _require_nucleus(self, nucleus: Nucleus) -> None:
        if nucleus not in self.nuclei:
            raise KeyError("nucleus does not belong to this constitution")

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Constitution) and _constitution_isomorphic(self, other)

    __hash__ = None

    @classmethod
    def from_indexed(
        cls,
        identities: Sequence[str | NuclearIdentity],
        connections: Iterable[tuple[int, int]],
        num_electrons: int,
    ) -> "Constitution":
        """Convenience input using temporary indices; the resulting object stores no indices."""
        nuclei = [Nucleus(x) for x in identities]
        edges = []
        for i, j in connections:
            if not _is_nonnegative_int(i) or not _is_nonnegative_int(j):
                raise ValueError("temporary input indices must be nonnegative integers")
            try:
                edges.append(Connection(nuclei[i], nuclei[j]))
            except IndexError as exc:
                raise ValueError("temporary input index is out of range") from exc
        return cls(nuclei, edges, num_electrons)

    # backwards-friendly convenience name
    from_data = from_indexed

    def all_unassigned_form(self) -> "ElectronAllocation":
        """σ_∅(M): n_v=0, all endpoint allocations 0, n_∅=N_e."""
        assigned = {v: 0 for v in self.nuclei}
        allocated = [
            AllocatedConnection(edge.endpoint1, edge.endpoint2, 0, 0)
            for edge in self.connections
        ]
        return ElectronAllocation(
            self.collection,
            assigned,
            allocated,
            self.num_electrons,
        )

    zero_allocation = all_unassigned_form

    def show(self, **kwargs) -> None:
        """Display an optional RDKit-rendered CCS view in IPython/Jupyter."""
        _show_ccs(self, **kwargs)

    def to_svg(self, **kwargs) -> str:
        """Return the rendered CCS view as an SVG string (requires RDKit)."""
        return _render_ccs_svg(self, **kwargs)

    def to_png(self, **kwargs) -> bytes:
        """Return PNG bytes with the CCS electron-count footer."""
        return _render_ccs_png(self, **kwargs)

    def save(self, path: str | Path, **kwargs) -> None:
        """Export a CCS view to .svg or .png, selected by filename suffix."""
        _save_ccs_visualization(self, path, **kwargs)

    def __repr__(self) -> str:
        counts = Counter(v.identity.label for v in self.nuclei)
        return (
            f"Constitution(nuclei={dict(sorted(counts.items()))}, "
            f"connections={len(self.connections)}, N_e={self.num_electrons})"
        )


@dataclass(frozen=True, eq=False, init=False, slots=True)
class ElectronAllocation:
    """D = (C_X, S, n, a, n_empty), intrinsically index-independent.

    Immutable after construction.  The assigned map is also read-only;
    constructor arguments are copied so external mutation cannot change D.
    """

    collection: NuclearCollection
    _assigned: Mapping[Nucleus, int]
    allocated_connections: frozenset[AllocatedConnection]
    num_unassigned_electrons: int

    def __init__(
        self,
        collection: NuclearCollection | Iterable[Nucleus],
        assigned_electrons: Mapping[Nucleus, int],
        allocated_connections: Iterable[AllocatedConnection] = (),
        num_unassigned_electrons: int = 0,
    ) -> None:
        normalized = (
            collection if isinstance(collection, NuclearCollection)
            else NuclearCollection(collection)
        )
        nuclei = normalized.nuclei
        assigned = dict(assigned_electrons)
        if set(assigned) != set(nuclei):
            raise ValueError("assigned_electrons must contain exactly one entry for every nucleus")
        if any(not _is_nonnegative_int(value) for value in assigned.values()):
            raise ValueError("assigned electron counts must be nonnegative integers")
        if not _is_nonnegative_int(num_unassigned_electrons):
            raise ValueError("num_unassigned_electrons must be a nonnegative integer")

        allocs = tuple(allocated_connections)
        seen_pairs: set[Connection] = set()
        totals = {v: 0 for v in nuclei}
        for edge in allocs:
            if not isinstance(edge, AllocatedConnection):
                raise TypeError("allocated_connections must contain AllocatedConnection objects")
            if edge.endpoint1 not in nuclei or edge.endpoint2 not in nuclei:
                raise ValueError("allocation endpoint does not belong to the nuclear collection")
            if edge.connection in seen_pairs:
                raise ValueError("parallel/duplicate connections are not allowed")
            seen_pairs.add(edge.connection)
            totals[edge.endpoint1] += edge.allocation1
            totals[edge.endpoint2] += edge.allocation2
        if any(totals[v] > assigned[v] for v in nuclei):
            raise ValueError(
                "electron-allocation form is not well-formed: "
                "endpoint allocations exceed assigned electrons"
            )

        object.__setattr__(self, "collection", normalized)
        object.__setattr__(self, "_assigned", MappingProxyType(assigned))
        object.__setattr__(self, "allocated_connections", frozenset(allocs))
        object.__setattr__(self, "num_unassigned_electrons", num_unassigned_electrons)

    @classmethod
    def empty(cls) -> "ElectronAllocation":
        return cls((), {}, (), 0)

    @classmethod
    def from_indexed(
        cls,
        identities: Sequence[str | NuclearIdentity],
        assigned_electrons: Sequence[int],
        allocated_connections: Iterable[tuple[int, int, int, int]] = (),
        num_unassigned_electrons: int = 0,
    ) -> "ElectronAllocation":
        if len(identities) != len(assigned_electrons):
            raise ValueError("identities and assigned_electrons must have the same length")
        nuclei = [Nucleus(x) for x in identities]
        assigned = {v: n for v, n in zip(nuclei, assigned_electrons)}
        allocs = []
        for i, j, ai, aj in allocated_connections:
            if not _is_nonnegative_int(i) or not _is_nonnegative_int(j):
                raise ValueError("temporary input indices must be nonnegative integers")
            try:
                allocs.append(AllocatedConnection(nuclei[i], nuclei[j], ai, aj))
            except IndexError as exc:
                raise ValueError("temporary input index is out of range") from exc
        return cls(nuclei, assigned, allocs, num_unassigned_electrons)

    from_data = from_indexed

    @property
    def nuclei(self) -> frozenset[Nucleus]:
        return self.collection.nuclei

    @property
    def connections(self) -> frozenset[Connection]:
        return frozenset(edge.connection for edge in self.allocated_connections)

    def assigned_electrons(self, nucleus: Nucleus) -> int:
        try:
            return self._assigned[nucleus]
        except KeyError as exc:
            raise KeyError("nucleus does not belong to this allocation form") from exc

    def assigned_mapping(self) -> dict[Nucleus, int]:
        return dict(self._assigned)

    def allocation_on(self, nucleus1: Nucleus, nucleus2: Nucleus) -> AllocatedConnection | None:
        for edge in self.allocated_connections:
            if edge.connection == Connection(nucleus1, nucleus2):
                return edge
        return None

    def allocated_sum(self, nucleus: Nucleus) -> int:
        if nucleus not in self.nuclei:
            raise KeyError("nucleus does not belong to this allocation form")
        return sum(
            edge.allocation_for(nucleus)
            for edge in self.allocated_connections
            if edge.connection.contains(nucleus)
        )

    def residual_electrons(self, nucleus: Nucleus) -> int:
        return self.assigned_electrons(nucleus) - self.allocated_sum(nucleus)

    def assignment_imbalance(self, nucleus: Nucleus) -> int:
        """S1.4 auxiliary q_v = Z(X(v)) - n_v (not Lewis formal charge)."""
        return nucleus.identity.atomic_number - self.assigned_electrons(nucleus)

    def formal_charge(self, nucleus: Nucleus) -> int:
        """Deprecated compatibility alias for assignment_imbalance().

        This does NOT compute an RDKit/Lewis formal charge or oxidation state.
        """
        warnings.warn(
            "ElectronAllocation.formal_charge() is deprecated; use "
            "assignment_imbalance() for the CCS q_v quantity. This is not "
            "necessarily a conventional Lewis formal charge.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.assignment_imbalance(nucleus)

    @property
    def num_electrons(self) -> int:
        return sum(self._assigned.values()) + self.num_unassigned_electrons

    @property
    def num_total_electron(self) -> int:
        return self.num_electrons

    def to_constitution(self) -> Constitution:
        """Q_el(D)."""
        return Constitution(self.collection, self.connections, self.num_electrons)

    def to_rdkit_mol(
        self,
        *,
        sanitize: bool = True,
        allow_unassigned: bool = False,
        order: Sequence[Nucleus] | None = None,
        include_ccs_properties: bool = True,
        strict: bool = False,
    ):
        """Best-effort RDKit projection of this electron-allocation form.

        Hydrogen is kept explicit and implicit H generation is disabled.
        Common symmetric allocations are mapped to ordinary bond
        multiplicities; (2,0)/(0,2) is mapped to a directed dative bond.
        Other allocation patterns retain the connection as an UNSPECIFIED
        RDKit bond and keep the original allocation as bond properties.

        Returns None when n_empty > 0 unless allow_unassigned=True.
        """
        return _electron_allocation_to_rdkit_mol(
            self,
            sanitize=sanitize,
            allow_unassigned=allow_unassigned,
            order=order,
            include_ccs_properties=include_ccs_properties,
            strict=strict,
        )

    def to_smiles(
        self,
        *,
        canonical: bool = True,
        isomeric: bool = True,
        remove_hs: bool = False,
        sanitize: bool = True,
        allow_unassigned: bool = False,
        strict: bool = False,
    ) -> str | None:
        """Return a best-effort RDKit SMILES, or None when projection fails.

        The default keeps CCS-explicit H atoms.  Set remove_hs=True for a
        conventional compact SMILES when RDKit can remove those hydrogens.
        """
        return _electron_allocation_to_smiles(
            self,
            canonical=canonical,
            isomeric=isomeric,
            remove_hs=remove_hs,
            sanitize=sanitize,
            allow_unassigned=allow_unassigned,
            strict=strict,
        )

    def __eq__(self, other: object) -> bool:
        return isinstance(other, ElectronAllocation) and _allocation_isomorphic(self, other)

    __hash__ = None

    def to_word(self, order: Sequence[Nucleus] | None = None) -> "Word":
        """Construct the exact S5 W(D) for the chosen temporary enumeration.

        The nuclei are installed in the supplied order, followed by n_empty
        Ele tokens and Con tokens sorted lexicographically by (i, j), i < j.
        Endpoint allocations are reoriented to match the increasing addresses.
        Without order, a noncanonical enumeration is chosen: intrinsic forms
        deliberately have no preferred nuclear numbering (S1.2, S3.2).
        """
        order = _validate_nucleus_order(self, order)
        address = {v: i for i, v in enumerate(order)}
        tokens: list[CoreTokenType] = [
            Nuc(v.identity, self.assigned_electrons(v)) for v in order
        ]
        tokens.extend(Ele() for _ in range(self.num_unassigned_electrons))

        connections: list[tuple[int, int, int, int]] = []
        for edge in self.allocated_connections:
            i, j = address[edge.endpoint1], address[edge.endpoint2]
            if i > j:
                i, j = j, i
            connections.append((
                i, j,
                edge.allocation_for(order[i]),
                edge.allocation_for(order[j]),
            ))
        tokens.extend(Con(i, j, ai, aj) for i, j, ai, aj in sorted(connections))
        return Word(tokens)

    def show(self, **kwargs) -> None:
        """Display an optional RDKit-rendered CCS view in IPython/Jupyter."""
        _show_ccs(self, **kwargs)

    def to_svg(self, **kwargs) -> str:
        """Return the rendered CCS view as an SVG string (requires RDKit)."""
        return _render_ccs_svg(self, **kwargs)

    def to_png(self, **kwargs) -> bytes:
        """Return PNG bytes with the CCS electron-count footer."""
        return _render_ccs_png(self, **kwargs)

    def save(self, path: str | Path, **kwargs) -> None:
        """Export a CCS view to .svg or .png, selected by filename suffix."""
        _save_ccs_visualization(self, path, **kwargs)

    def __repr__(self) -> str:
        counts = Counter(v.identity.label for v in self.nuclei)
        return (
            f"ElectronAllocation(nuclei={dict(sorted(counts.items()))}, "
            f"connections={len(self.allocated_connections)}, "
            f"N_e={self.num_electrons}, n_empty={self.num_unassigned_electrons})"
        )


# ---------------------------------------------------------------------------
# Index-independent equality helpers
# ---------------------------------------------------------------------------

def _adjacency(connections: Iterable[Connection], nuclei: Iterable[Nucleus]) -> dict[Nucleus, set[Nucleus]]:
    adj = {v: set() for v in nuclei}
    for edge in connections:
        adj[edge.endpoint1].add(edge.endpoint2)
        adj[edge.endpoint2].add(edge.endpoint1)
    return adj


def _constitution_isomorphic(a: Constitution, b: Constitution) -> bool:
    if a.num_electrons != b.num_electrons:
        return False
    if len(a.nuclei) != len(b.nuclei) or len(a.connections) != len(b.connections):
        return False
    if a.collection.identity_multiset() != b.collection.identity_multiset():
        return False

    adj_a = _adjacency(a.connections, a.nuclei)
    adj_b = _adjacency(b.connections, b.nuclei)

    nodes_a = tuple(a.nuclei)
    nodes_b = tuple(b.nuclei)

    candidates = {
        v: [
            w for w in nodes_b
            if w.identity == v.identity and len(adj_b[w]) == len(adj_a[v])
        ]
        for v in nodes_a
    }
    if any(not c for c in candidates.values()):
        return False

    mapping: dict[Nucleus, Nucleus] = {}
    used: set[Nucleus] = set()

    def choose_next() -> Nucleus:
        remaining = [v for v in nodes_a if v not in mapping]
        return min(remaining, key=lambda v: sum(w not in used for w in candidates[v]))

    def dfs() -> bool:
        if len(mapping) == len(nodes_a):
            return True
        v = choose_next()
        for w in candidates[v]:
            if w in used:
                continue
            ok = True
            for v2, w2 in mapping.items():
                if ((v2 in adj_a[v]) != (w2 in adj_b[w])):
                    ok = False
                    break
            if not ok:
                continue
            mapping[v] = w
            used.add(w)
            if dfs():
                return True
            used.remove(w)
            del mapping[v]
        return False

    return dfs()


def _allocation_isomorphic(a: ElectronAllocation, b: ElectronAllocation) -> bool:
    if a.num_unassigned_electrons != b.num_unassigned_electrons:
        return False
    if len(a.nuclei) != len(b.nuclei):
        return False
    if len(a.allocated_connections) != len(b.allocated_connections):
        return False

    con_a = a.to_constitution()
    con_b = b.to_constitution()
    if con_a.num_electrons != con_b.num_electrons:
        return False

    adj_a = _adjacency(a.connections, a.nuclei)
    adj_b = _adjacency(b.connections, b.nuclei)

    nodes_a = tuple(a.nuclei)
    nodes_b = tuple(b.nuclei)

    candidates = {
        v: [
            w for w in nodes_b
            if (
                w.identity == v.identity
                and b.assigned_electrons(w) == a.assigned_electrons(v)
                and len(adj_b[w]) == len(adj_a[v])
            )
        ]
        for v in nodes_a
    }
    if any(not c for c in candidates.values()):
        return False

    def edge_map(form: ElectronAllocation) -> dict[frozenset[Nucleus], AllocatedConnection]:
        return {
            frozenset((edge.endpoint1, edge.endpoint2)): edge
            for edge in form.allocated_connections
        }

    ea = edge_map(a)
    eb = edge_map(b)
    mapping: dict[Nucleus, Nucleus] = {}
    used: set[Nucleus] = set()

    def choose_next() -> Nucleus:
        remaining = [v for v in nodes_a if v not in mapping]
        return min(remaining, key=lambda v: sum(w not in used for w in candidates[v]))

    def dfs() -> bool:
        if len(mapping) == len(nodes_a):
            return True
        v = choose_next()
        for w in candidates[v]:
            if w in used:
                continue
            ok = True
            for v2, w2 in mapping.items():
                key_a = frozenset((v, v2))
                key_b = frozenset((w, w2))
                has_a = key_a in ea
                has_b = key_b in eb
                if has_a != has_b:
                    ok = False
                    break
                if has_a:
                    edge_a = ea[key_a]
                    edge_b = eb[key_b]
                    if edge_a.allocation_for(v) != edge_b.allocation_for(w):
                        ok = False
                        break
                    if edge_a.allocation_for(v2) != edge_b.allocation_for(w2):
                        ok = False
                        break
            if not ok:
                continue
            mapping[v] = w
            used.add(w)
            if dfs():
                return True
            used.remove(w)
            del mapping[v]
        return False

    return dfs()


# ---------------------------------------------------------------------------
# S2: core construction alphabet and replay
# ---------------------------------------------------------------------------

class CoreToken:
    pass


@dataclass(frozen=True, init=False, slots=True)
class Nuc(CoreToken):
    """Nuc(X; n, u).  If n is omitted, n = Z(X)."""

    identity: NuclearIdentity
    n: int
    u: int

    def __init__(
        self,
        identity: str | NuclearIdentity,
        n: int | None = None,
        u: int = 0,
    ) -> None:
        identity = NuclearIdentity.coerce(identity)
        if n is None:
            n = identity.atomic_number
        if not _is_nonnegative_int(n):
            raise ValueError("assigned-electron count n must be a nonnegative integer")
        if not _is_nonnegative_int(u):
            raise ValueError("reservation u must be a nonnegative integer")
        object.__setattr__(self, "identity", identity)
        object.__setattr__(self, "n", n)
        object.__setattr__(self, "u", u)

    @property
    def assigned_electrons(self) -> int:
        return self.n

    @property
    def reservation(self) -> int:
        return self.u

    def __repr__(self) -> str:
        if self.n == self.identity.atomic_number and self.u == 0:
            return f"Nuc({self.identity.label})"
        if self.u == 0:
            return f"Nuc({self.identity.label}; {self.n})"
        return f"Nuc({self.identity.label}; {self.n}, {self.u})"


@dataclass(frozen=True, slots=True)
class Ele(CoreToken):
    def __repr__(self) -> str:
        return "Ele"


@dataclass(frozen=True, slots=True)
class Con(CoreToken):
    i: int
    j: int
    ai: int
    aj: int

    def __post_init__(self) -> None:
        if not _is_int(self.i) or not _is_int(self.j):
            raise TypeError("Con addresses i and j must be integers")
        if not _is_nonnegative_int(self.ai) or not _is_nonnegative_int(self.aj):
            raise ValueError("Con endpoint allocations must be nonnegative integers")

    @property
    def address1(self) -> int:
        return self.i

    @property
    def address2(self) -> int:
        return self.j

    @property
    def allocation1(self) -> int:
        return self.ai

    @property
    def allocation2(self) -> int:
        return self.aj

    def __repr__(self) -> str:
        return f"Con({self.i}, {self.j}; {self.ai}, {self.aj})"


CoreTokenType = Nuc | Ele | Con


@dataclass(frozen=True, slots=True)
class ExecutionState:
    """Intrinsic form plus temporary construction metadata."""

    form: ElectronAllocation
    address_table: tuple[Nucleus, ...]
    reservations: tuple[int, ...]

    @classmethod
    def empty(cls) -> "ExecutionState":
        return cls(ElectronAllocation.empty(), (), ())

    @property
    def num_nuclei(self) -> int:
        return len(self.address_table)

    def resolve_address(self, address: int) -> Nucleus | None:
        if not _is_int(address):
            raise TypeError("address must be an integer")
        m = self.num_nuclei
        if 0 <= address < m:
            return self.address_table[address]
        if -m <= address < 0:
            return self.address_table[m + address]
        return None

    def reservation_at(self, address: int) -> int:
        nucleus = self.resolve_address(address)
        if nucleus is None:
            raise KeyError("address does not resolve")
        # address_table is injective
        for k, v in enumerate(self.address_table):
            if v is nucleus:
                return self.reservations[k]
        raise RuntimeError("inconsistent address table")

    def apply(self, token: CoreTokenType) -> "ExecutionState":
        if isinstance(token, Nuc):
            nucleus = Nucleus(token.identity)
            nuclei = set(self.form.nuclei)
            nuclei.add(nucleus)
            assigned = self.form.assigned_mapping()
            assigned[nucleus] = token.n
            new_form = ElectronAllocation(
                nuclei,
                assigned,
                self.form.allocated_connections,
                self.form.num_unassigned_electrons,
            )
            return ExecutionState(
                new_form,
                self.address_table + (nucleus,),
                self.reservations + (token.u,),
            )

        if isinstance(token, Ele):
            new_form = ElectronAllocation(
                self.form.collection,
                self.form.assigned_mapping(),
                self.form.allocated_connections,
                self.form.num_unassigned_electrons + 1,
            )
            return ExecutionState(new_form, self.address_table, self.reservations)

        if isinstance(token, Con):
            v = self.resolve_address(token.i)
            w = self.resolve_address(token.j)

            # Inapplicable Con: the entire state, including addresses, is unchanged.
            if v is None or w is None or v is w:
                return self

            edge_key = Connection(v, w)
            if edge_key in self.form.connections:
                return self

            if token.ai > self.form.residual_electrons(v):
                return self
            if token.aj > self.form.residual_electrons(w):
                return self

            new_allocs = set(self.form.allocated_connections)
            new_allocs.add(AllocatedConnection(v, w, token.ai, token.aj))
            new_form = ElectronAllocation(
                self.form.collection,
                self.form.assigned_mapping(),
                new_allocs,
                self.form.num_unassigned_electrons,
            )
            return ExecutionState(new_form, self.address_table, self.reservations)

        raise TypeError(f"not a core token: {token!r}")


class Word:
    """A finite word over the core alphabet A."""

    __slots__ = ("tokens",)

    def __init__(self, tokens: Iterable[CoreTokenType] = ()) -> None:
        tokens = tuple(tokens)
        if not all(isinstance(t, (Nuc, Ele, Con)) for t in tokens):
            raise TypeError("Word may contain only Nuc, Ele, and Con tokens")
        self.tokens = tokens

    def __iter__(self) -> Iterator[CoreTokenType]:
        return iter(self.tokens)

    def __len__(self) -> int:
        return len(self.tokens)

    def append(self, token: CoreTokenType) -> "Word":
        if not isinstance(token, (Nuc, Ele, Con)):
            raise TypeError("Word may contain only Nuc, Ele, and Con tokens")
        return Word(self.tokens + (token,))

    def replay_state(self) -> ExecutionState:
        state = ExecutionState.empty()
        for token in self.tokens:
            state = state.apply(token)
        return state

    def replay(self) -> ElectronAllocation:
        # Temporary address table and u metadata are forgotten here.
        return self.replay_state().form

    def to_constitution(self) -> Constitution:
        return self.replay().to_constitution()

    def to_rdkit_mol(self, **kwargs):
        return self.replay().to_rdkit_mol(**kwargs)

    def to_smiles(self, **kwargs) -> str | None:
        return self.replay().to_smiles(**kwargs)

    def show(self, **kwargs) -> None:
        """Display an optional RDKit-rendered CCS view in IPython/Jupyter."""
        _show_ccs(self, **kwargs)

    def to_svg(self, **kwargs) -> str:
        """Return the rendered CCS view as an SVG string (requires RDKit)."""
        return _render_ccs_svg(self, **kwargs)

    def to_png(self, **kwargs) -> bytes:
        """Return PNG bytes with the CCS electron-count footer."""
        return _render_ccs_png(self, **kwargs)

    def save(self, path: str | Path, **kwargs) -> None:
        """Export a CCS view to .svg or .png, selected by filename suffix."""
        _save_ccs_visualization(self, path, **kwargs)

    def __repr__(self) -> str:
        return "Word([" + ", ".join(repr(t) for t in self.tokens) + "])"


# ---------------------------------------------------------------------------
# S4: End completion
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class End:
    def __repr__(self) -> str:
        return "End"


@dataclass(frozen=True, slots=True)
class CompletionInfo:
    address: int
    nucleus: Nucleus
    in_completion_domain: bool
    eligible: bool
    b: int | None
    B: int
    T: int
    E: int | None
    L: int | None
    h: int


def _in_completion_domain(identity: NuclearIdentity) -> bool:
    # S4.1: underlying elements in groups 13-17.  H itself is excluded.
    return identity.symbol != "H" and identity.group in {13, 14, 15, 16, 17}


def _reference_valence(identity: NuclearIdentity) -> int:
    if not _in_completion_domain(identity):
        raise ValueError("identity is outside the S4 completion domain")
    assert identity.group is not None
    return identity.group - 10


def _core_count(identity: NuclearIdentity) -> int:
    return identity.atomic_number - _reference_valence(identity)


def _position_cost(edge: AllocatedConnection) -> int:
    # rho_s = max(1, floor((a_i,s + a_j,s)/2))
    return max(1, (edge.allocation1 + edge.allocation2) // 2)


def completion_plan(state: ExecutionState) -> tuple[CompletionInfo, ...]:
    result = []
    form = state.form

    for address, nucleus in enumerate(state.address_table):
        incident = [
            edge
            for edge in form.allocated_connections
            if edge.connection.contains(nucleus)
        ]
        B = sum(edge.allocation_for(nucleus) for edge in incident)
        T = sum(_position_cost(edge) for edge in incident)

        if not _in_completion_domain(nucleus.identity):
            result.append(
                CompletionInfo(address, nucleus, False, False, None, B, T, None, None, 0)
            )
            continue

        b = form.assigned_electrons(nucleus) - _core_count(nucleus.identity)
        u = state.reservations[address]
        E = b - u - B

        eligible = (0 <= b <= 8) and (0 <= u <= b) and (E >= 0)
        if not eligible:
            result.append(
                CompletionInfo(address, nucleus, True, False, b, B, T, E, None, 0)
            )
            continue

        L = max(4 - u - T, ceil(E / 2))
        h = min(E, 2 * L - E)
        result.append(
            CompletionInfo(address, nucleus, True, True, b, B, T, E, L, h)
        )

    return tuple(result)


# ---------------------------------------------------------------------------
# S4.5: finite macros
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class Arg:
    """Reference to one macro invocation argument by zero-based position."""
    index: int

    def __post_init__(self) -> None:
        if not _is_nonnegative_int(self.index):
            raise ValueError("macro argument index must be a nonnegative integer")


MacroEndpoint = int | Arg


@dataclass(frozen=True, slots=True)
class MacroCon:
    """Connection token inside a macro body.

    Endpoints may be Arg(k) or negative relative addresses.
    Nonnegative literal addresses are intentionally forbidden.
    """

    i: MacroEndpoint
    j: MacroEndpoint
    ai: int
    aj: int

    def __post_init__(self) -> None:
        for endpoint in (self.i, self.j):
            if isinstance(endpoint, Arg):
                continue
            if not _is_int(endpoint):
                raise TypeError("macro endpoints must be Arg(k) or integer relative addresses")
            if endpoint >= 0:
                raise ValueError(
                    "nonnegative address literals are not permitted inside a macro body; "
                    "use Arg(k) for external nuclei"
                )
        if not _is_nonnegative_int(self.ai) or not _is_nonnegative_int(self.aj):
            raise ValueError("macro endpoint allocations must be nonnegative integers")


MacroBodyToken = Nuc | Ele | MacroCon | Con


@dataclass(frozen=True, slots=True)
class MacroCall:
    macro: "Macro"
    arguments: tuple[int, ...]

    def __repr__(self) -> str:
        args = ", ".join(map(str, self.arguments))
        return f"Macro({self.macro.name}" + (f"; {args})" if args else ")")


class Macro:
    """A finite abbreviation of ordinary core construction tokens."""

    __slots__ = ("name", "num_args", "body")

    def __init__(
        self,
        name: str,
        body: Iterable[MacroBodyToken],
        num_args: int = 0,
    ) -> None:
        if not isinstance(name, str) or not name:
            raise ValueError("macro name must be a nonempty string")
        if not _is_nonnegative_int(num_args):
            raise ValueError("num_args must be a nonnegative integer")

        normalized: list[Nuc | Ele | MacroCon] = []
        for token in body:
            if isinstance(token, (Nuc, Ele)):
                normalized.append(token)
            elif isinstance(token, MacroCon):
                normalized.append(token)
            elif isinstance(token, Con):
                # Ordinary Con is allowed as shorthand only when both addresses
                # are negative relative addresses.
                if token.i >= 0 or token.j >= 0:
                    raise ValueError(
                        "nonnegative address literals are not permitted inside a macro body"
                    )
                normalized.append(MacroCon(token.i, token.j, token.ai, token.aj))
            else:
                raise TypeError("macro body may contain only Nuc, Ele, and macro-local Con")

        for token in normalized:
            if isinstance(token, MacroCon):
                for endpoint in (token.i, token.j):
                    if isinstance(endpoint, Arg) and endpoint.index >= num_args:
                        raise ValueError(
                            f"macro body references Arg({endpoint.index}) "
                            f"but num_args={num_args}"
                        )

        self.name = name
        self.num_args = num_args
        self.body = tuple(normalized)

    def __call__(self, *arguments: int) -> MacroCall:
        if len(arguments) != self.num_args:
            raise ValueError(
                f"macro {self.name!r} expects {self.num_args} arguments, "
                f"got {len(arguments)}"
            )
        for value in arguments:
            if not _is_nonnegative_int(value):
                raise ValueError("macro invocation arguments must be nonnegative addresses")
        return MacroCall(self, tuple(arguments))

    def expand(self, arguments: Sequence[int], pre_invocation_nuclei: int) -> tuple[CoreTokenType, ...]:
        if len(arguments) != self.num_args:
            raise ValueError("wrong number of macro arguments")
        if not _is_nonnegative_int(pre_invocation_nuclei):
            raise ValueError("pre_invocation_nuclei must be nonnegative")
        if any(not _is_nonnegative_int(k) for k in arguments):
            raise ValueError("macro arguments must be nonnegative integer addresses")

        # An out-of-range external address makes the invocation contribute no tokens.
        if any(k >= pre_invocation_nuclei for k in arguments):
            return ()

        def resolve(endpoint: MacroEndpoint) -> int:
            if isinstance(endpoint, Arg):
                return arguments[endpoint.index]
            return endpoint  # negative relative address

        result: list[CoreTokenType] = []
        for token in self.body:
            if isinstance(token, (Nuc, Ele)):
                result.append(token)
            else:
                result.append(
                    Con(
                        resolve(token.i),
                        resolve(token.j),
                        token.ai,
                        token.aj,
                    )
                )
        return tuple(result)

    def __repr__(self) -> str:
        return f"MacroDefinition({self.name!r}, num_args={self.num_args})"


# ---------------------------------------------------------------------------
# Practical input: core tokens + macros + End
# ---------------------------------------------------------------------------

PracticalItem = CoreTokenType | MacroCall | End


class PracticalWord:
    """Finite practical construction input.

    - Core tokens replay with the S2 semantics.
    - Macro calls expand to finite ordinary core tokens.
    - The first End performs completion exactly once and terminates execution.
    - Every item after the first End is ignored.
    - With no End, the input has ordinary core replay after macro expansion.
    """

    __slots__ = ("items",)

    def __init__(self, items: Iterable[PracticalItem] = ()) -> None:
        items = tuple(items)
        if not all(isinstance(x, (Nuc, Ele, Con, MacroCall, End)) for x in items):
            raise TypeError("unsupported practical construction item")
        self.items = items

    def __iter__(self) -> Iterator[PracticalItem]:
        return iter(self.items)

    def append(self, item: PracticalItem) -> "PracticalWord":
        return PracticalWord(self.items + (item,))

    @property
    def has_end(self) -> bool:
        return any(isinstance(item, End) for item in self.items)

    def _effective_prefix(self) -> tuple[tuple[CoreTokenType, ...], ExecutionState, bool]:
        state = ExecutionState.empty()
        tokens: list[CoreTokenType] = []

        for item in self.items:
            if isinstance(item, End):
                return tuple(tokens), state, True

            if isinstance(item, MacroCall):
                expanded = item.macro.expand(item.arguments, state.num_nuclei)
                for token in expanded:
                    tokens.append(token)
                    state = state.apply(token)
                continue

            # ordinary core token
            tokens.append(item)
            state = state.apply(item)

        return tuple(tokens), state, False

    def pre_end_state(self) -> ExecutionState:
        """State immediately before the first End, or after all items if no End occurs."""
        return self._effective_prefix()[1]

    def completion_plan(self) -> tuple[CompletionInfo, ...]:
        """Hypothetical/actual End plan for the current effective prefix."""
        return completion_plan(self.pre_end_state())

    def completion_counts(self) -> tuple[int, ...]:
        return tuple(info.h for info in self.completion_plan())

    def expand(self) -> Word:
        prefix, state, has_end = self._effective_prefix()
        if not has_end:
            return Word(prefix)

        tokens: list[CoreTokenType] = list(prefix)
        for info in completion_plan(state):
            for _ in range(info.h):
                tokens.append(Nuc("H"))
                tokens.append(Con(info.address, -1, 1, 1))
        return Word(tokens)

    def replay(self) -> ElectronAllocation:
        return self.expand().replay()

    def to_constitution(self) -> Constitution:
        return self.replay().to_constitution()

    def to_rdkit_mol(self, **kwargs):
        return self.replay().to_rdkit_mol(**kwargs)

    def to_smiles(self, **kwargs) -> str | None:
        return self.replay().to_smiles(**kwargs)

    def complete(self) -> "PracticalWord":
        """Append End only when the input does not already contain one."""
        if self.has_end:
            return self
        return PracticalWord(self.items + (End(),))

    def show(self, **kwargs) -> None:
        """Display an optional RDKit-rendered CCS view in IPython/Jupyter."""
        _show_ccs(self, **kwargs)

    def to_svg(self, **kwargs) -> str:
        """Return the rendered CCS view as an SVG string (requires RDKit)."""
        return _render_ccs_svg(self, **kwargs)

    def to_png(self, **kwargs) -> bytes:
        """Return PNG bytes with the CCS electron-count footer."""
        return _render_ccs_png(self, **kwargs)

    def save(self, path: str | Path, **kwargs) -> None:
        """Export a CCS view to .svg or .png, selected by filename suffix."""
        _save_ccs_visualization(self, path, **kwargs)

    def __repr__(self) -> str:
        return "PracticalWord([" + ", ".join(repr(x) for x in self.items) + "])"


# ---------------------------------------------------------------------------
# S4.6: optional search profile (proposal policy only)
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class SearchProfile:
    """Optional proposal filter; it never changes CCS replay semantics.

    All fields are optional.  A token rejected by this policy remains a valid
    CCS construction token whenever it belongs to the underlying alphabet.
    """

    allowed_identities: frozenset[NuclearIdentity] | None = None
    allowed_assigned_counts: frozenset[int] | None = None
    allowed_reservations: frozenset[int] | None = None
    allowed_endpoint_allocations: frozenset[tuple[int, int]] | None = None
    allowed_macro_names: frozenset[str] | None = None
    max_nuclei: int | None = None
    allow_ele: bool = True
    extra_filter: Callable[[PracticalItem, ExecutionState], bool] | None = None

    @classmethod
    def create(
        cls,
        *,
        allowed_identities: Iterable[str | NuclearIdentity] | None = None,
        allowed_assigned_counts: Iterable[int] | None = None,
        allowed_reservations: Iterable[int] | None = None,
        allowed_endpoint_allocations: Iterable[tuple[int, int]] | None = None,
        allowed_macro_names: Iterable[str] | None = None,
        max_nuclei: int | None = None,
        allow_ele: bool = True,
        extra_filter: Callable[[PracticalItem, ExecutionState], bool] | None = None,
    ) -> "SearchProfile":
        identities = (
            None
            if allowed_identities is None
            else frozenset(NuclearIdentity.coerce(x) for x in allowed_identities)
        )
        counts = None if allowed_assigned_counts is None else frozenset(allowed_assigned_counts)
        reservations = None if allowed_reservations is None else frozenset(allowed_reservations)
        allocations = (
            None
            if allowed_endpoint_allocations is None
            else frozenset(tuple(pair) for pair in allowed_endpoint_allocations)
        )
        macros = None if allowed_macro_names is None else frozenset(allowed_macro_names)

        if max_nuclei is not None and not _is_nonnegative_int(max_nuclei):
            raise ValueError("max_nuclei must be a nonnegative integer or None")
        return cls(
            identities,
            counts,
            reservations,
            allocations,
            macros,
            max_nuclei,
            allow_ele,
            extra_filter,
        )

    def allows(self, proposal: PracticalItem, state: ExecutionState) -> bool:
        if isinstance(proposal, Nuc):
            if self.max_nuclei is not None and state.num_nuclei >= self.max_nuclei:
                return False
            if self.allowed_identities is not None and proposal.identity not in self.allowed_identities:
                return False
            if self.allowed_assigned_counts is not None and proposal.n not in self.allowed_assigned_counts:
                return False
            if self.allowed_reservations is not None and proposal.u not in self.allowed_reservations:
                return False

        elif isinstance(proposal, Con):
            if (
                self.allowed_endpoint_allocations is not None
                and (proposal.ai, proposal.aj) not in self.allowed_endpoint_allocations
            ):
                return False

        elif isinstance(proposal, Ele):
            if not self.allow_ele:
                return False

        elif isinstance(proposal, MacroCall):
            if self.allowed_macro_names is not None and proposal.macro.name not in self.allowed_macro_names:
                return False
            if self.max_nuclei is not None:
                expanded = proposal.macro.expand(proposal.arguments, state.num_nuclei)
                added = sum(isinstance(token, Nuc) for token in expanded)
                if state.num_nuclei + added > self.max_nuclei:
                    return False

        # End is never restricted by the built-in fields.

        if self.extra_filter is not None and not self.extra_filter(proposal, state):
            return False
        return True



# ---------------------------------------------------------------------------
# Optional visualization and export, in this same implementation module.
# No depiction field enters S1-S4 molecular identity or construction replay.
# ---------------------------------------------------------------------------

class CCSVisualizationError(RuntimeError):
    """Optional CCS depiction could not be rendered."""


def _visual_source(obj, *, view: str = "auto", order=None):
    """Return an intrinsic source and a purely temporary display enumeration."""
    if view not in {"auto", "constitution", "allocation"}:
        raise ValueError("view must be 'auto', 'constitution', or 'allocation'")

    if isinstance(obj, (Word, PracticalWord)):
        state = (
            obj.replay_state() if isinstance(obj, Word)
            else obj.expand().replay_state()
        )
        source = state.form
        default_order = state.address_table
    elif isinstance(obj, ElectronAllocation):
        source = obj
        default_order = tuple(obj.nuclei)
    elif isinstance(obj, Constitution):
        source = obj
        default_order = tuple(obj.nuclei)
    else:
        raise TypeError("visualization supports Constitution, ElectronAllocation, "
                        "Word, and PracticalWord")

    if order is None:
        display_order = tuple(default_order)
    else:
        display_order = tuple(order)
        if (len(display_order) != len(source.nuclei)
                or set(display_order) != set(source.nuclei)):
            raise ValueError("order must enumerate each nucleus of this object exactly once")

    if isinstance(source, Constitution):
        if view == "allocation":
            raise TypeError("Constitution has no electron-allocation data")
        return source, display_order, "constitution"
    if view == "constitution":
        return source.to_constitution(), display_order, "constitution"
    return source, display_order, "allocation"


def _ccs_visual_style(Chem, edge: AllocatedConnection):
    """Purely illustrative RDKit stroke; never defines a CCS bond order."""
    a, b = edge.allocation1, edge.allocation2
    if (a, b) == (1, 1):
        return Chem.BondType.SINGLE, None, None, True
    if (a, b) == (2, 2):
        return Chem.BondType.DOUBLE, None, None, True
    if (a, b) == (3, 3):
        return Chem.BondType.TRIPLE, None, None, True
    if (a, b) == (2, 0):
        return Chem.BondType.DATIVE, edge.endpoint1, edge.endpoint2, True
    if (a, b) == (0, 2):
        return Chem.BondType.DATIVE, edge.endpoint2, edge.endpoint1, True
    return Chem.BondType.ZERO, None, None, False


def _build_ccs_visual_mol(
    source,
    order: Sequence[Nucleus],
    *,
    source_kind: str,
    show_hydrogens: bool,
    show_indices: bool,
    diagnostics: bool,
    label_all_allocations: bool,
    allocation_label_mode: str = "bond",
):
    """Build an unsanitized drawing model without mutating CCS source.

    The internal 'bond' default retains legacy bondNote expectations; the
    public drawing interface defaults to endpoint-local numeric labels.
    """
    if allocation_label_mode not in {"endpoint", "bond", "none"}:
        raise ValueError("allocation_label_mode must be endpoint, bond or none")
    Chem = _load_rdkit()
    visible = tuple(
        v for v in order if show_hydrogens or v.identity.symbol != "H"
    )
    visible_set = set(visible)
    display_index = {v: i for i, v in enumerate(order)}
    rw = Chem.RWMol()
    indices = {}

    for nucleus in visible:
        atom = Chem.Atom(nucleus.identity.symbol)
        atom.SetNoImplicit(True)
        atom.SetIntProp("_CCSTemporaryIndex", display_index[nucleus])
        if nucleus.identity.isotope is not None:
            atom.SetIsotope(nucleus.identity.isotope)
        note = [f"#{display_index[nucleus]}"] if show_indices else []
        if source_kind == "allocation":
            # IMPORTANT: CCS q_v is the S1.4 auxiliary assignment imbalance,
            # not RDKit/Lewis formal charge. Do not call SetFormalCharge(q_v).
            q_v = source.assignment_imbalance(nucleus)
            n_v = source.assigned_electrons(nucleus)
            r_v = source.residual_electrons(nucleus)
            atom.SetIntProp("_CCSAssignmentImbalance", q_v)
            atom.SetIntProp("_CCSAssignedElectrons", n_v)
            atom.SetIntProp("_CCSResidualElectrons", r_v)
            # This valence-electron dot count is an S1 allocation diagnostic,
            # not a claim about spin, lone-pair pairing or orbital occupancy.
            # It intentionally excludes H and elements outside groups 13–17.
            if _in_completion_domain(nucleus.identity):
                valence_remainder = r_v - _core_count(nucleus.identity)
                atom.SetIntProp("_CCSResidualValenceElectrons", valence_remainder)
            if diagnostics:
                note.append(f"n={n_v} r={r_v} q={q_v:+d}")
        if note:
            atom.SetProp("atomNote", "  ".join(note))
        indices[nucleus] = rw.AddAtom(atom)

    if source_kind == "constitution":
        for edge in sorted(
            source.connections,
            key=lambda e: tuple(sorted((
                display_index[e.endpoint1], display_index[e.endpoint2]
            ))),
        ):
            if edge.endpoint1 in visible_set and edge.endpoint2 in visible_set:
                # In this view, a plain stroke means direct CCS connection,
                # not conventional bond order or electron sharing.
                rw.AddBond(indices[edge.endpoint1], indices[edge.endpoint2],
                           Chem.BondType.SINGLE)
    else:
        for edge in sorted(
            source.allocated_connections,
            key=lambda e: tuple(sorted((
                display_index[e.endpoint1], display_index[e.endpoint2]
            ))),
        ):
            if edge.endpoint1 not in visible_set or edge.endpoint2 not in visible_set:
                continue
            bond_type, donor, acceptor, conventional = _ccs_visual_style(Chem, edge)
            begin = donor if donor is not None else edge.endpoint1
            end = acceptor if acceptor is not None else edge.endpoint2
            rw.AddBond(indices[begin], indices[end], bond_type)
            bond = rw.GetBondBetweenAtoms(indices[begin], indices[end])
            if bond is not None:
                bond.SetIntProp("_CCSBeginAllocation", edge.allocation_for(begin))
                bond.SetIntProp("_CCSEndAllocation", edge.allocation_for(end))
                if label_all_allocations or not conventional:
                    if allocation_label_mode == "bond":
                        i, j = sorted((display_index[edge.endpoint1],
                                       display_index[edge.endpoint2]))
                        bond.SetProp(
                            "bondNote",
                            f"{i}:{edge.allocation_for(order[i])}, "
                            f"{j}:{edge.allocation_for(order[j])}",
                        )
                    elif allocation_label_mode == "endpoint":
                        # Display the two counts near their corresponding
                        # atoms, following RDKit's actual drawn orientation.
                        bond.SetBoolProp("_CCSDisplayEndpointAllocations", True)
                    if edge.allocation1 + edge.allocation2 == 1:
                        bond.SetProp("_CCSOneElectron", "true")

    mol = rw.GetMol()
    mol.SetProp("_CCSView", source_kind)
    mol.SetBoolProp("_CCSLossyHydrogenOmission",
                    not show_hydrogens and any(v.identity.symbol == "H" for v in order))
    return mol



def _ccs_endpoint_annotation_positions(drawer, mol):
    """Place each allocation next to its *own* bond endpoint in drawing space.

    A deterministic local collision search keeps numerical labels away from
    other atoms and previously placed numbers without changing the underlying
    CCS allocations, drawing coordinates, or nucleus-index correspondence.
    """
    annotations = []
    atom_centers = [
        (drawer.GetDrawCoords(i).x, drawer.GetDrawCoords(i).y)
        for i in range(mol.GetNumAtoms())
    ]
    min_x, max_x = 13.0, float(drawer.Width()) - 13.0
    min_y, max_y = 14.0, float(drawer.Height()) - 14.0
    for bond in mol.GetBonds():
        if not (bond.HasProp("_CCSDisplayEndpointAllocations")
                and bond.GetBoolProp("_CCSDisplayEndpointAllocations")):
            continue
        begin, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        p, q = atom_centers[begin], atom_centers[end]
        dx, dy = q[0] - p[0], q[1] - p[1]
        length = hypot(dx, dy)
        if length < 1.0e-8:
            nx, ny = 0.0, 1.0
        else:
            nx, ny = -dy / length, dx / length

        for atom_idx, own, other, electrons in (
            (begin, p, q, bond.GetIntProp("_CCSBeginAllocation")),
            (end, q, p, bond.GetIntProp("_CCSEndAllocation")),
        ):
            best = None
            text_radius = 6.0 + 3.5 * len(str(electrons))
            # Offset from the nucleus along its own connection, then shift
            # perpendicular to avoid nearby nuclei, other labels and edges.
            for fraction in (0.23, 0.31, 0.16):
                for side in (1, -1):
                    for displacement in (13.0, 21.0, 29.0, 38.0):
                        x = own[0] + fraction * (other[0] - own[0]) + side * displacement * nx
                        y = own[1] + fraction * (other[1] - own[1]) + side * displacement * ny
                        penalty = 0.18 * displacement + 28.0 * abs(fraction - 0.23)
                        if not min_x < x < max_x or not min_y < y < max_y:
                            penalty += 2.0e4
                        for atom_number, (ax, ay) in enumerate(atom_centers):
                            distance = hypot(x - ax, y - ay)
                            threshold = (18.0 if atom_number == atom_idx else 21.0) + text_radius * 0.55
                            if distance < threshold:
                                penalty += 170.0 * (threshold - distance + 1.0)
                        for _, _, previous_electrons, px, py in annotations:
                            target = 13.0 + text_radius + 3.5 * len(str(previous_electrons))
                            separation = hypot(x - px, y - py)
                            if separation < target:
                                penalty += 220.0 * (target - separation + 1.0)
                        candidate = (penalty, x, y)
                        if best is None or candidate < best:
                            best = candidate
            atom = mol.GetAtomWithIdx(atom_idx)
            annotations.append((
                bond.GetIdx(),
                atom.GetIntProp("_CCSTemporaryIndex"),
                electrons,
                best[1],
                best[2],
            ))
    return annotations


def _ccs_segment_intersections(mol) -> int:
    """Count straight-connection conflicts in the temporary 2D depiction.

    Connected edges may meet at their shared nucleus; unrelated segments
    cannot cross, coincide, or pass through an unrelated nucleus.
    """
    if mol.GetNumAtoms() < 3 or mol.GetNumBonds() == 0:
        return 0
    conf = mol.GetConformer()
    pts = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    segments = [
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()
    ]

    def orient(p, q, r):
        return (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x)

    def hit(p, q, r):
        length = hypot(q.x - p.x, q.y - p.y)
        if length < 1e-10:
            return hypot(r.x - p.x, r.y - p.y) < 1e-7
        if abs(orient(p, q, r)) > 1e-7 * length:
            return False
        dot = ((r.x - p.x) * (q.x - p.x)
               + (r.y - p.y) * (q.y - p.y)) / (length * length)
        return 1e-6 < dot < 1 - 1e-6

    conflicts = 0
    for k, (a, b) in enumerate(segments):
        p, q = pts[a], pts[b]
        for c, d in segments[k + 1:]:
            if len({a, b, c, d}) < 4:
                continue
            r, t = pts[c], pts[d]
            o1, o2 = orient(p, q, r), orient(p, q, t)
            o3, o4 = orient(r, t, p), orient(r, t, q)
            crossing = (o1 * o2 < -1e-12 and o3 * o4 < -1e-12)
            touching = (hit(p, q, r) or hit(p, q, t)
                        or hit(r, t, p) or hit(r, t, q))
            if crossing or touching:
                conflicts += 1
        # An edge passing through an unconnected atom is also ambiguous.
        for i, other in enumerate(pts):
            if i not in (a, b) and hit(p, q, other):
                conflicts += 1
    return conflicts



def _ccs_layout_2d(mol, layout: str) -> None:
    """Choose an intelligible 2D nuclear graph without chemical exceptions.

    Score polygon shape, fully incident ring hubs, relative bond lengths
    and separation. All candidate coordinates use a temporary untyped
    connection graph; electron allocations cannot alter the geometry.
    """
    from math import isfinite
    from statistics import median
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    import networkx as nx

    graph = nx.Graph()
    graph.add_nodes_from(range(mol.GetNumAtoms()))
    graph.add_edges_from(
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()
    )
    planar = nx.check_planarity(graph)[0]
    mol.SetBoolProp("_CCSPlanarGraph", bool(planar))
    # Ignore hub-to-edge triangles using only graph degree, not element type.
    ring_graph = graph.subgraph(
        v for v, degree in graph.degree() if degree <= 5
    )
    short_rings = [
        ring for ring in nx.cycle_basis(ring_graph)
        if 4 <= len(ring) <= 8
    ]
    # Any vertex joined to an entire short cycle belongs near its projected
    # center. This rule uses graph incidence alone, not element or bond type.
    ring_hubs = tuple(
        tuple(
            hub for hub in graph
            if hub not in ring and all(graph.has_edge(hub, v) for v in ring)
        )
        for ring in short_rings
    )

    def capture(candidate):
        if not candidate.GetNumAtoms():
            return ()
        if not candidate.GetNumConformers():
            return None
        con = candidate.GetConformer()
        xy = tuple(
            (float(con.GetAtomPosition(i).x),
             float(con.GetAtomPosition(i).y))
            for i in range(mol.GetNumAtoms())
        )
        return xy if all(isfinite(x) and isfinite(y)
                         for x, y in xy) else None

    def install(xy):
        con = Chem.Conformer(mol.GetNumAtoms())
        for i, (x, y) in enumerate(xy):
            con.SetAtomPosition(i, (float(x), float(y), 0.0))
        mol.RemoveAllConformers()
        mol.AddConformer(con, assignId=True)

    def score(xy, intersections):
        if not graph.number_of_edges():
            return 0.0, 0.0, 0.0
        lengths = [
            hypot(xy[i][0] - xy[j][0], xy[i][1] - xy[j][1])
            for i, j in graph.edges()
        ]
        typical = median(lengths)
        if typical <= 1e-10 or not isfinite(typical):
            return float("inf"), float("inf"), float("inf")
        ring_cost = 0.0
        hub_cost = 0.0
        for ring, hubs in zip(short_rings, ring_hubs):
            n = len(ring)
            positions = [xy[i] for i in ring]
            edges = [
                hypot(positions[k][0] - positions[(k+1) % n][0],
                      positions[k][1] - positions[(k+1) % n][1])
                for k in range(n)
            ]
            avg = sum(edges) / n
            edge_variance = sum((v - avg)**2 for v in edges) / n
            center_x = sum(p[0] for p in positions) / n
            center_y = sum(p[1] for p in positions) / n
            radii = [
                hypot(x - center_x, y - center_y) for x, y in positions
            ]
            radius_avg = sum(radii) / n
            radius_variance = sum(
                (v - radius_avg)**2 for v in radii
            ) / n
            area = abs(sum(
                positions[k][0] * positions[(k+1) % n][1]
                - positions[(k+1) % n][0] * positions[k][1]
                for k in range(n)
            )) / 2.0
            ring_cost += (
                80.0 * edge_variance / max(avg*avg, 1e-12)
                + 45.0 * radius_variance / max(radius_avg**2, 1e-12)
                + 80.0 * max(0.0, 0.40 - area/typical**2)**2
            )
            # Perfectly regular rings can still be misleading if a vertex
            # incident to their whole perimeter is depicted far outside.
            # Allow a small off-center tolerance in units of the ring radius.
            for hub in hubs:
                offset = hypot(
                    xy[hub][0] - center_x, xy[hub][1] - center_y
                ) / max(radius_avg, 1e-8 * typical)
                hub_cost += 12.0 * max(0.0, offset - 0.55)**2
        # Keep metal/long-range incidence strokes: extreme lengths are
        # capped so they do not outweigh an intact Cp ring or cage face.
        bond_cost = sum(
            4.0 * max(0.0, 0.55 - d/typical)**2
            + 0.10 * min(2.5, max(0.0, d/typical - 3.0))**2
            for d in lengths
        )
        overlap_cost = 0.0
        for i in range(len(xy)):
            for j in range(i+1, len(xy)):
                if graph.has_edge(i, j):
                    continue
                distance = hypot(xy[i][0]-xy[j][0],
                                 xy[i][1]-xy[j][1])
                overlap_cost += 25.0 * max(
                    0.0, 0.72-distance/typical
                )**2
        return (ring_cost + hub_cost + bond_cost + overlap_cost
                + 0.08*intersections), ring_cost, hub_cost

    def spread_overlapping_leaves(points):
        """Separate coincident terminal atoms without moving cycle vertices.

        Based solely on graph degree and 2D geometry; no H/metal-specific
        exceptions. Other candidate layouts and CCS allocations are unchanged.
        """
        if not graph.number_of_edges():
            return points
        from math import atan2, sin, cos, pi
        p = [list(pair) for pair in points]
        bond_lengths = [
            hypot(p[i][0]-p[j][0],p[i][1]-p[j][1])
            for i,j in graph.edges()
        ]
        scale = median(bond_lengths)
        if scale <= 1e-8:
            return points
        cx = sum(v[0] for v in p)/len(p)
        cy = sum(v[1] for v in p)/len(p)
        for parent in sorted(graph):
            leaves = sorted(
                v for v in graph.neighbors(parent) if graph.degree(v) == 1
            )
            if len(leaves) < 2:
                continue
            if all(hypot(
                p[leaves[i]][0]-p[leaves[j]][0],
                p[leaves[i]][1]-p[leaves[j]][1]
            ) >= 0.45*scale
                   for i in range(len(leaves))
                   for j in range(i+1,len(leaves))):
                continue
            x0,y0 = p[parent]
            dx,dy = 0.0,0.0
            for other in graph.neighbors(parent):
                if other in leaves:
                    continue
                ux,uy = p[other][0]-x0,p[other][1]-y0
                d = hypot(ux,uy)
                if d>1e-8:
                    dx-=ux/d
                    dy-=uy/d
            if hypot(dx,dy)<1e-5:
                dx,dy=x0-cx,y0-cy
            if hypot(dx,dy)<1e-5:
                dx,dy=0.0,1.0
            baseline=atan2(dy,dx)
            for k,leaf in enumerate(leaves):
                angle=baseline+(k-(len(leaves)-1)/2.0)*1.15
                best=None
                for delta in (0.0,-0.25,0.25,-0.55,0.55,1.,-1.):
                    for length in (1.,1.25,1.50):
                        th=angle+delta
                        x=x0+length*scale*cos(th)
                        y=y0+length*scale*sin(th)
                        separation=min(
                            hypot(x-p[v][0],y-p[v][1])
                            for v in graph if v not in (parent,leaf)
                        )
                        value=(
                            150*max(0.0,1.0*scale-separation)**2
                            + 3*abs(delta)+2*abs(length-1.0)
                        )
                        if best is None or value<best[0]:
                            best=(value,x,y)
                p[leaf]=[best[1],best[2]]
        return tuple((x,y) for x,y in p)

    # A highly incident graph hub may deform otherwise regular rings when
    # *all* edges are treated as conventional single bonds by CoordGen.
    # In one extra purely topological candidate, weaken high-degree edges
    # during COORDINATE GENERATION ONLY, then restore the original graph for
    # every visible connection and electron-allocation glyph.
    max_degree=max((degree for _,degree in graph.degree()),default=0)
    candidates = []
    methods = (
        ("coordgen", "coordgen-light-hubs", "rdkit")
        if layout == "auto" and max_degree >= 8
        else (("coordgen","rdkit") if layout=="auto" else (layout,))
    )
    for method in methods:
        try:
            # Derive coordinates from the connection graph alone. A copy
            # prevents DATIVE/ZERO and other electron-allocation bond glyphs
            # from changing the 2D geometry of the same CCS constitution.
            working = Chem.Mol(mol)
            working.RemoveAllConformers()
            for edge in working.GetBonds():
                i, j = edge.GetBeginAtomIdx(), edge.GetEndAtomIdx()
                weakened = (
                    method == "coordgen-light-hubs"
                    and max(graph.degree(i), graph.degree(j)) >= 8
                )
                edge.SetBondType(
                    Chem.BondType.ZERO if weakened else Chem.BondType.SINGLE
                )
            working.UpdatePropertyCache(strict=False)
            if method in {"coordgen", "coordgen-light-hubs"}:
                from rdkit.Chem import rdCoordGen
                rdCoordGen.AddCoords(working)
            else:
                rdDepictor.Compute2DCoords(working)
            xy = capture(working)
            if xy is not None:
                # Preserve skeleton cycles; separate only overlapping leaves.
                xy = spread_overlapping_leaves(xy)
                install(xy)
                crossing = _ccs_segment_intersections(mol)
                energy, rings, hub_penalty = score(xy, crossing)
                candidates.append((
                    (energy, method not in {"coordgen", "coordgen-light-hubs"}),
                    method, xy, crossing, rings, hub_penalty
                ))
        except (ImportError, RuntimeError, ValueError) as exc:
            if layout != "auto":
                raise CCSVisualizationError(
                    f"2D depiction with {layout!r} failed"
                ) from exc

    if not candidates:
        # No 3D molecular geometry, only a generic graph fallback.
        pos = nx.spring_layout(graph, seed=23, iterations=90)
        xy = tuple(tuple(map(float, pos[i]))
                   for i in range(mol.GetNumAtoms()))
        install(xy)
        crossing = _ccs_segment_intersections(mol)
        energy, rings, hub_penalty = score(xy, crossing)
        candidates.append((
            (energy, True), "spring", xy, crossing, rings, hub_penalty
        ))

    ranking, method, xy, crossing, rings, hub_penalty = min(
        candidates, key=lambda candidate: candidate[0]
    )
    install(xy)
    mol.SetProp("_CCSLayoutMethod", method)
    mol.SetDoubleProp("_CCSLayoutScore", float(ranking[0]))
    mol.SetProp(
        "_CCSLayoutCandidateScores",
        ";".join(f"{candidate[1]}={candidate[0][0]:.10g}" for candidate in candidates),
    )
    mol.SetDoubleProp("_CCSRingDistortion", float(rings))
    mol.SetDoubleProp("_CCSRingHubPenalty", float(hub_penalty))
    mol.SetIntProp("_CCS2DCrossings", crossing)
    mol.SetBoolProp("_CCSNonPlanar2D", not planar)


def _ccs_nuclear_annotation_positions(drawer, mol, endpoint_annotations):
    """Layout auxiliary q signs and residual valence dots near each nucleus.

    Uses the same temporary RDKit atom addresses as endpoint labels.
    No atom charge or radical is assigned to RDKit's chemistry fields.
    """
    centers = [
        (drawer.GetDrawCoords(i).x, drawer.GetDrawCoords(i).y)
        for i in range(mol.GetNumAtoms())
    ]
    segments = [
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        for bond in mol.GetBonds()
    ]
    occupied = [
        (x, y, 11.0 + len(str(value)) * 2.0)
        for _, _, value, x, y in endpoint_annotations
    ]
    labels = []
    width, height = float(drawer.Width()), float(drawer.Height())

    def segment_distance(x, y, a, b):
        ax, ay = centers[a]
        bx, by = centers[b]
        dx, dy = bx - ax, by - ay
        divisor = dx * dx + dy * dy
        if divisor < 1e-10:
            return hypot(x - ax, y - ay)
        t = min(1.0, max(0.0, ((x - ax) * dx + (y - ay) * dy) / divisor))
        return hypot(x - (ax + t * dx), y - (ay + t * dy))

    for atom in mol.GetAtoms():
        if not atom.HasProp("_CCSAssignedElectrons"):
            continue
        idx = atom.GetIdx()
        q_value = atom.GetIntProp("_CCSAssignmentImbalance")
        remaining = (atom.GetIntProp("_CCSResidualValenceElectrons")
                     if atom.HasProp("_CCSResidualValenceElectrons") else None)
        items = []
        if q_value:
            sign = "+" if q_value > 0 else "\u2212"
            glyph = (str(abs(q_value)) if abs(q_value) > 1 else "") + sign
            items.append(("charge", glyph, q_value))
        if remaining is not None and remaining > 0:
            if remaining <= 8:
                count = remaining
                while count:
                    n = min(2, count)
                    items.append(("residual", "\u2022" * n, n))
                    count -= n
            else:
                items.append(("residual_overflow", f"r={remaining}", remaining))
        if not items:
            continue
        cx, cy = centers[idx]
        for kind, glyph, count in items:
            candidates = []
            # Residual electrons belong visually next to their nucleus,
            # unlike a charge label which can be slightly farther away.
            # Candidate positions avoid the bond fan and endpoint labels.
            radii = (
                (11.0, 15.0, 21.0, 29.0) if kind == "residual"
                else (18.0, 26.0, 36.0, 48.0)
            )
            for radius in radii:
                for angle_index in range(16):
                    theta = pi * angle_index / 8.0
                    x = cx + radius * cos(theta)
                    y = cy + radius * sin(theta)
                    score = 0.16 * radius + 0.25 * angle_index
                    if not (12.0 < x < width - 12.0
                            and 12.0 < y < height - 12.0):
                        score += 1e6
                    for j, (ax, ay) in enumerate(centers):
                        threshold = (
                            (15.0 if j != idx else 9.0)
                            if kind == "residual"
                            else (20.0 if j != idx else 16.0)
                        )
                        gap = hypot(x - ax, y - ay)
                        if gap < threshold:
                            score += 80.0 * (threshold - gap + 1.0)
                    for a, b in segments:
                        distance = segment_distance(x, y, a, b)
                        bond_clearance = 7.0 if kind == "residual" else 11.0
                        if distance < bond_clearance:
                            score += 70.0 * (bond_clearance + 1.0 - distance)
                    radius_text = 8.0 + len(glyph) * 3.5
                    for ox, oy, ore in occupied:
                        gap = hypot(x - ox, y - oy)
                        if gap < radius_text + ore:
                            score += 120.0 * (radius_text + ore - gap + 1.0)
                    candidates.append((score, x, y))
            _, x, y = min(candidates)
            labels.append((idx, kind, glyph, count, x, y))
            occupied.append((x, y, 8.0 + len(glyph) * 3.5))
    return labels


def _ccs_footer(source, source_kind: str, *,
                show_total_charge: bool, show_hydrogens: bool) -> list[str]:
    lines = [f"N_e = {source.num_electrons}"]
    if show_total_charge:
        constitution = (source if source_kind == "constitution"
                        else source.to_constitution())
        lines.append(f"Q = {constitution.total_charge:+d}")
    if source_kind == "allocation" and source.num_unassigned_electrons:
        lines.append(f"n_empty = {source.num_unassigned_electrons}")
    if not show_hydrogens and any(v.identity.symbol == "H" for v in source.nuclei):
        lines.append("H omitted - lossy view (not full CCS graph)")
    return lines


def _ccs_render_options(
    obj, *, view="auto", order=None, show_hydrogens=True,
    show_indices=False, diagnostics=False, label_all_allocations=False,
    allocation_label_mode="endpoint", layout="auto",
    show_total_charge=False, width=1000, height=700,
    bond_length_scale=1.35, annotation_scale=0.8, footer_font_size=18,
):
    if not (_is_int(width) and _is_int(height) and width >= 120 and height >= 150):
        raise ValueError("width and height must be integers (at least 120 and 150)")
    if not _is_int(footer_font_size) or footer_font_size <= 0:
        raise ValueError("footer_font_size must be a positive integer")
    if bond_length_scale <= 0 or annotation_scale <= 0:
        raise ValueError("bond_length_scale and annotation_scale must be positive")
    if not isinstance(show_hydrogens, bool):
        raise TypeError("show_hydrogens must be a boolean; default is True")
    if allocation_label_mode not in {"endpoint", "bond", "none"}:
        raise ValueError("allocation_label_mode must be endpoint, bond or none")
    if layout not in {"auto", "coordgen", "rdkit"}:
        raise ValueError("layout must be 'auto', 'coordgen' or 'rdkit'")

    source, display_order, kind = _visual_source(obj, view=view, order=order)
    mol = _build_ccs_visual_mol(
        source, display_order, source_kind=kind,
        show_hydrogens=show_hydrogens, show_indices=show_indices,
        diagnostics=diagnostics, label_all_allocations=label_all_allocations,
        allocation_label_mode=allocation_label_mode,
    )
    # Enforce the planarity constraint independently of RDKit bond styles.
    mol.UpdatePropertyCache(strict=False)
    _ccs_layout_2d(mol, layout)
    if bond_length_scale != 1 and mol.GetNumAtoms():
        conf = mol.GetConformer()
        for index in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(index)
            conf.SetAtomPosition(
                index, (p.x * bond_length_scale, p.y * bond_length_scale, 0.0)
            )
    footer = _ccs_footer(
        source, kind, show_total_charge=show_total_charge,
        show_hydrogens=show_hydrogens,
    )
    if mol.GetBoolProp("_CCSNonPlanar2D"):
        footer.append("Non-planar graph: crossings unavoidable in 2D")
    footer_height = max(40, 8 + len(footer) * (footer_font_size + 3))
    drawing_height = height - footer_height
    if drawing_height < 80:
        raise ValueError("height is too small for the requested CCS footer")
    return mol, footer, drawing_height, annotation_scale, footer_font_size


def _render_ccs_svg(
    obj, *, view="auto", order=None, show_hydrogens=True,
    show_indices=False, diagnostics=False, label_all_allocations=False,
    allocation_label_mode="endpoint", layout="auto",
    show_total_charge=False, width=1000, height=700,
    bond_length_scale=1.35, annotation_scale=0.8, footer_font_size=18,
) -> str:
    """Export the full-graph CCS view as SVG; does not require IPython."""
    from rdkit.Chem.Draw import rdMolDraw2D
    import html

    mol, footer, body_height, annotation_scale, font_size = _ccs_render_options(
        obj, view=view, order=order, show_hydrogens=show_hydrogens,
        show_indices=show_indices, diagnostics=diagnostics,
        label_all_allocations=label_all_allocations,
        allocation_label_mode=allocation_label_mode, layout=layout,
        show_total_charge=show_total_charge,
        width=width, height=height, bond_length_scale=bond_length_scale,
        annotation_scale=annotation_scale, footer_font_size=footer_font_size,
    )
    drawer = rdMolDraw2D.MolDraw2DSVG(width, body_height)
    opts = drawer.drawOptions()
    # Carbon is normally hidden by chemical drawing conventions. When a CCS
    # atom bears a derived charge sign or residual valence dots, show its
    # element symbol so that these labels have an unambiguous visible owner.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C" or not atom.HasProp("_CCSAssignedElectrons"):
            continue
        q = atom.GetIntProp("_CCSAssignmentImbalance")
        remainder = (atom.GetIntProp("_CCSResidualValenceElectrons")
                     if atom.HasProp("_CCSResidualValenceElectrons") else 0)
        if q or remainder > 0:
            opts.atomLabels[atom.GetIdx()] = "C"
    opts.annotationFontScale = annotation_scale
    opts.bondLineWidth = 2.2
    opts.additionalAtomLabelPadding = 0.12
    opts.padding = 0.05
    drawer.DrawMolecule(mol)
    endpoint_annotations = _ccs_endpoint_annotation_positions(drawer, mol)
    nuclear_annotations = _ccs_nuclear_annotation_positions(
        drawer, mol, endpoint_annotations
    )
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    # Individual electron counts remain next to the endpoint that contributed
    # them, and exact temporary endpoint indices are retained as SVG metadata.
    endpoint_glyphs = ["<g class='ccs-endpoint-allocations'>"]
    for bond_idx, atom_idx, electrons, x, y in endpoint_annotations:
        endpoint_glyphs.append(
            f"<text class='ccs-endpoint-allocation' data-bond='{bond_idx}' "
            f"data-atom='{atom_idx}' data-electrons='{electrons}' "
            f"x='{x:.2f}' y='{y + 4.5:.2f}' text-anchor='middle' "
            "font-family='sans-serif' font-size='13px' font-weight='bold' "
            "fill='#405bbb' stroke='white' stroke-width='3.5' "
            "stroke-linejoin='round' paint-order='stroke fill'>"
            f"{electrons}</text>"
        )
    endpoint_glyphs.append("</g>")
    nuclear_glyphs = ["<g class='ccs-nuclear-annotations'>"]
    for atom_index, kind, glyph, count, x, y in nuclear_annotations:
        if kind == "residual":
            # Vector circles avoid all UTF-8 / bullet-font serialization and
            # SVG-client glyph-rendering differences. Metadata remains exact.
            nuclear_glyphs.append(
                f"<g class='ccs-nuclear-annotation' "
                f"data-atom='{atom_index}' data-kind='residual' "
                f"data-count='{count}'>"
            )
            for index in range(count):
                dot_x = x + (index - (count - 1) / 2.0) * 7.4
                nuclear_glyphs.append(
                    f"<circle cx='{dot_x:.2f}' cy='{y:.2f}' "
                    "r='2.7' fill='#222222'/>"
                )
            nuclear_glyphs.append("</g>")
        else:
            nuclear_glyphs.append(
                f"<text class='ccs-nuclear-annotation' "
                f"data-atom='{atom_index}' data-kind='{kind}' "
                f"data-count='{count}' x='{x:.2f}' y='{y + 4.3:.2f}' "
                "text-anchor='middle' font-family='sans-serif' "
                "font-size='16px' fill='#222222' stroke='white' "
                "stroke-width='3px' stroke-linejoin='round' "
                "paint-order='stroke fill'>"
                f"{html.escape(glyph)}</text>"
            )
    nuclear_glyphs.append("</g>")

    # Extend the RDKit canvas with a separate full-system metadata footer.
    svg, count_height = re.subn(
        rf"height=['\"]{body_height}px['\"]",
        f"height='{height}px'", svg, count=1,
    )
    svg, count_viewbox = re.subn(
        rf"viewBox=['\"]0 0 {width} {body_height}['\"]",
        f"viewBox='0 0 {width} {height}'", svg, count=1,
    )
    if count_height != 1 or count_viewbox != 1:
        raise CCSVisualizationError("Could not extend RDKit SVG canvas")

    y = body_height + font_size + 2
    fragments = [
        f"<rect x='0' y='{body_height}' width='{width}' "
        f"height='{height - body_height}' fill='white'/>"
    ]
    for line in footer:
        fragments.append(
            f"<text x='{width - 14}' y='{y}' text-anchor='end' "
            f"font-family='sans-serif' font-size='{font_size}px'>"
            f"{html.escape(line)}</text>"
        )
        y += font_size + 3
    # RDKit may convert annotation text to SVG paths. Preserve the exact
    # S1.4 auxiliary values as machine-readable XML metadata as well.
    source, ordered, kind = _visual_source(obj, view=view, order=order)
    if kind == "allocation":
        imbalance_entries = [
            f"{i}:{v.identity.label}:q={source.assignment_imbalance(v):+d}"
            for i, v in enumerate(ordered)
        ]
        fragments.insert(
            0,
            "<desc>CCS assignment imbalance (not Lewis formal charge): "
            + html.escape("; ".join(imbalance_entries))
            + "</desc>",
        )
    return svg.replace(
        "</svg>",
        "\n" + "\n".join(endpoint_glyphs + nuclear_glyphs + fragments)
        + "\n</svg>",
    )


def _render_ccs_png(
    obj, *, view="auto", order=None, show_hydrogens=True,
    show_indices=False, diagnostics=False, label_all_allocations=False,
    allocation_label_mode="endpoint", layout="auto",
    show_total_charge=False, width=1000, height=700,
    bond_length_scale=1.35, annotation_scale=0.8, footer_font_size=18,
) -> bytes:
    """Export PNG with the same CCS metadata semantics (RDKit + Pillow)."""
    from io import BytesIO
    from rdkit.Chem.Draw import rdMolDraw2D
    try:
        from PIL import Image, ImageDraw, ImageFont
    except ImportError as exc:
        raise CCSVisualizationError(
            "Pillow is required for PNG export; install the visualization extra"
        ) from exc

    mol, footer, body_height, annotation_scale, font_size = _ccs_render_options(
        obj, view=view, order=order, show_hydrogens=show_hydrogens,
        show_indices=show_indices, diagnostics=diagnostics,
        label_all_allocations=label_all_allocations,
        allocation_label_mode=allocation_label_mode, layout=layout,
        show_total_charge=show_total_charge,
        width=width, height=height, bond_length_scale=bond_length_scale,
        annotation_scale=annotation_scale, footer_font_size=footer_font_size,
    )
    drawer = rdMolDraw2D.MolDraw2DCairo(width, body_height)
    opts = drawer.drawOptions()
    # Carbon is normally hidden by chemical drawing conventions. When a CCS
    # atom bears a derived charge sign or residual valence dots, show its
    # element symbol so that these labels have an unambiguous visible owner.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C" or not atom.HasProp("_CCSAssignedElectrons"):
            continue
        q = atom.GetIntProp("_CCSAssignmentImbalance")
        remainder = (atom.GetIntProp("_CCSResidualValenceElectrons")
                     if atom.HasProp("_CCSResidualValenceElectrons") else 0)
        if q or remainder > 0:
            opts.atomLabels[atom.GetIdx()] = "C"
    opts.annotationFontScale = annotation_scale
    opts.bondLineWidth = 2.2
    opts.additionalAtomLabelPadding = 0.12
    opts.padding = 0.05
    drawer.DrawMolecule(mol)
    endpoint_annotations = _ccs_endpoint_annotation_positions(drawer, mol)
    nuclear_annotations = _ccs_nuclear_annotation_positions(
        drawer, mol, endpoint_annotations
    )
    drawer.FinishDrawing()
    canvas = Image.new("RGB", (width, height), "white")
    with Image.open(BytesIO(drawer.GetDrawingText())) as body:
        canvas.paste(body.convert("RGB"), (0, 0))
    text_draw = ImageDraw.Draw(canvas)
    # PNG uses exactly the same atom-to-endpoint association as SVG.
    try:
        endpoint_font = ImageFont.truetype("DejaVuSans-Bold.ttf", 13)
    except OSError:
        endpoint_font = ImageFont.load_default()
    for _, _, electrons, x, y in endpoint_annotations:
        number = str(electrons)
        bounds = text_draw.textbbox((0, 0), number, font=endpoint_font)
        tw, th = bounds[2] - bounds[0], bounds[3] - bounds[1]
        left, top = x - tw / 2, y - th / 2
        text_draw.rectangle(
            (left - 2, top - 2, left + tw + 2, top + th + 2),
            fill="white",
        )
        text_draw.text(
            (left - bounds[0], top - bounds[1]), number,
            fill="#405bbb", font=endpoint_font,
        )
    try:
        nuclear_font = ImageFont.truetype("DejaVuSans.ttf", 16)
    except OSError:
        nuclear_font = ImageFont.load_default()
    for atom_index, kind, glyph, count, x, y in nuclear_annotations:
        if kind == "residual":
            for index in range(count):
                dx = x + (index - (count - 1) / 2.0) * 7.4
                text_draw.ellipse(
                    (dx - 2.7, y - 2.7, dx + 2.7, y + 2.7),
                    fill="#222222",
                )
        else:
            box = text_draw.textbbox((0, 0), glyph, font=nuclear_font)
            tw, th = box[2] - box[0], box[3] - box[1]
            left, top = x - tw / 2, y - th / 2
            text_draw.rectangle(
                (left - 2, top - 2, left + tw + 2, top + th + 2),
                fill="white",
            )
            text_draw.text(
                (left - box[0], top - box[1]), glyph, fill="#222222",
                font=nuclear_font,
            )
    try:
        font = ImageFont.truetype("DejaVuSans.ttf", font_size)
    except OSError:
        font = ImageFont.load_default()
    y = body_height + 2
    for line in footer:
        bbox = text_draw.textbbox((0, 0), line, font=font)
        text_width = bbox[2] - bbox[0]
        text_draw.text((width - 14 - text_width, y), line, fill="black", font=font)
        y += font_size + 3
    buffer = BytesIO()
    canvas.save(buffer, format="PNG")
    return buffer.getvalue()


def _show_ccs(obj, **kwargs) -> None:
    """Display a CCS SVG in a notebook; rendering itself does not need IPython."""
    try:
        from IPython.display import SVG, display
    except ImportError as exc:
        raise CCSVisualizationError(
            "IPython is required for show(); use to_svg()/save() without it"
        ) from exc
    display(SVG(data=_render_ccs_svg(obj, **kwargs)))


def _save_ccs_visualization(obj, path: str | Path, **kwargs) -> None:
    destination = Path(path)
    suffix = destination.suffix.lower()
    if suffix == ".svg":
        destination.write_text(_render_ccs_svg(obj, **kwargs), encoding="utf-8")
    elif suffix == ".png":
        destination.write_bytes(_render_ccs_png(obj, **kwargs))
    else:
        raise ValueError("CCS visualization export requires .svg or .png filename")


__all__ = [
    # S1
    "NuclearIdentity",
    "Nucleus",
    "NuclearCollection",
    "Connection",
    "AllocatedConnection",
    "Constitution",
    "ElectronAllocation",
    "RDKitConversionError",
    "CCSVisualizationError",
    # S2
    "Nuc",
    "Ele",
    "Con",
    "ExecutionState",
    "Word",
    # S4
    "End",
    "CompletionInfo",
    "completion_plan",
    "Arg",
    "MacroCon",
    "Macro",
    "MacroCall",
    "PracticalWord",
    "SearchProfile",
]