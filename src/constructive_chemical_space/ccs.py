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
from math import ceil
import re
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

        # q_v = Z(X(v)) - n_v is the natural first projection to an RDKit
        # atom formal charge.  This is only a downstream representation;
        # CCS itself does not identify q_v with a conventional Lewis charge
        # for every allocation.
        atom.SetFormalCharge(form.formal_charge(nucleus))

        if include_ccs_properties:
            atom.SetIntProp("_CCSAssignedElectrons", form.assigned_electrons(nucleus))
            atom.SetIntProp("_CCSResidualElectrons", form.residual_electrons(nucleus))
            atom.SetIntProp("_CCSAuxCharge", form.formal_charge(nucleus))
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

    Sanitization is tried first using q_v = Z - n as RDKit formal charges.
    If that fails and a (2,0)/(0,2) dative allocation is present, a second
    attempt shifts donor/acceptor charges by +1/-1, respectively.
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


class Constitution:
    """M = (C_X, S, N_e), intrinsically index-independent."""

    __slots__ = ("collection", "connections", "num_electrons")

    def __init__(
        self,
        collection: NuclearCollection | Iterable[Nucleus],
        connections: Iterable[Connection] = (),
        num_electrons: int = 0,
    ) -> None:
        self.collection = (
            collection if isinstance(collection, NuclearCollection)
            else NuclearCollection(collection)
        )
        if not _is_nonnegative_int(num_electrons):
            raise ValueError("num_electrons must be a nonnegative integer")
        self.num_electrons = num_electrons

        conns = frozenset(connections)
        for edge in conns:
            if not isinstance(edge, Connection):
                raise TypeError("connections must contain Connection objects")
            if edge.endpoint1 not in self.collection.nuclei or edge.endpoint2 not in self.collection.nuclei:
                raise ValueError("every connection endpoint must belong to the nuclear collection")
        self.connections = conns

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

    def __repr__(self) -> str:
        counts = Counter(v.identity.label for v in self.nuclei)
        return (
            f"Constitution(nuclei={dict(sorted(counts.items()))}, "
            f"connections={len(self.connections)}, N_e={self.num_electrons})"
        )


class ElectronAllocation:
    """D = (C_X, S, n, a, n_empty), intrinsically index-independent."""

    __slots__ = ("collection", "_assigned", "allocated_connections", "num_unassigned_electrons")

    def __init__(
        self,
        collection: NuclearCollection | Iterable[Nucleus],
        assigned_electrons: Mapping[Nucleus, int],
        allocated_connections: Iterable[AllocatedConnection] = (),
        num_unassigned_electrons: int = 0,
    ) -> None:
        self.collection = (
            collection if isinstance(collection, NuclearCollection)
            else NuclearCollection(collection)
        )
        nuclei = self.collection.nuclei

        assigned = dict(assigned_electrons)
        if set(assigned) != set(nuclei):
            raise ValueError("assigned_electrons must contain exactly one entry for every nucleus")
        if any(not _is_nonnegative_int(value) for value in assigned.values()):
            raise ValueError("assigned electron counts must be nonnegative integers")

        if not _is_nonnegative_int(num_unassigned_electrons):
            raise ValueError("num_unassigned_electrons must be a nonnegative integer")

        allocs = tuple(allocated_connections)
        seen_pairs: set[Connection] = set()
        for edge in allocs:
            if not isinstance(edge, AllocatedConnection):
                raise TypeError("allocated_connections must contain AllocatedConnection objects")
            if edge.endpoint1 not in nuclei or edge.endpoint2 not in nuclei:
                raise ValueError("allocation endpoint does not belong to the nuclear collection")
            if edge.connection in seen_pairs:
                raise ValueError("parallel/duplicate connections are not allowed")
            seen_pairs.add(edge.connection)

        self._assigned = assigned
        self.allocated_connections = frozenset(allocs)
        self.num_unassigned_electrons = num_unassigned_electrons

        for nucleus in nuclei:
            if self.allocated_sum(nucleus) > self.assigned_electrons(nucleus):
                raise ValueError(
                    "electron-allocation form is not well-formed: "
                    "endpoint allocations exceed assigned electrons"
                )

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

    def formal_charge(self, nucleus: Nucleus) -> int:
        """Auxiliary q_v = Z(X(v)) - n_v; not an independently stored field."""
        return nucleus.identity.atomic_number - self.assigned_electrons(nucleus)

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
        """Construct one exact witness word as in S3.2.

        The supplied order is a temporary enumeration only.  If omitted, an
        arbitrary finite enumeration of V is chosen.
        """
        if order is None:
            order = tuple(self.nuclei)
        else:
            order = tuple(order)
            if len(order) != len(self.nuclei) or set(order) != set(self.nuclei):
                raise ValueError("order must enumerate every nucleus exactly once")

        address = {v: i for i, v in enumerate(order)}
        tokens: list[CoreToken] = [
            Nuc(v.identity, self.assigned_electrons(v)) for v in order
        ]
        tokens.extend(Ele() for _ in range(self.num_unassigned_electrons))
        for edge in self.allocated_connections:
            tokens.append(
                Con(
                    address[edge.endpoint1],
                    address[edge.endpoint2],                    edge.allocation1,
                    edge.allocation2,
                )
            )
        return Word(tokens)

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