"""
Constructive Chemical Space (CCS)
Reference data structures and replay semantics for SI Sections S1-S4.

This module follows the specification in:
    supporting_information_ks_r001_260930.pdf

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
from typing import Callable, Iterable, Iterator, Mapping, Sequence


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
            if not _is_int(i) or not _is_int(j):
                raise TypeError("temporary input indices must be integers")
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
        """Display this constitution in a Jupyter/IPython notebook."""
        kwargs.setdefault("view", "constitution")
        _show_ccs(self, **kwargs)

    def __repr__(self) -> str:
        counts = Counter(v.identity.label for v in self.nuclei)
        return (
            f"Constitution(nuclei={dict(sorted(counts.items()))}, "
            f"connections={len(self.connections)}, N_e={self.num_electrons})"
        )


class ElectronAllocation:
    """D = (C_X, S, n, a, n_empty), intrinsically index-independent."""

    __slots__ = ("collection", "_assigned", "allocated_connections", "num_unassigned_electron")

    def __init__(
        self,
        collection: NuclearCollection | Iterable[Nucleus],
        assigned_electrons: Mapping[Nucleus, int],
        allocated_connections: Iterable[AllocatedConnection] = (),
        num_unassigned_electron: int = 0,
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

        if not _is_nonnegative_int(num_unassigned_electron):
            raise ValueError("num_unassigned_electron must be a nonnegative integer")

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
        self.num_unassigned_electron = num_unassigned_electron

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
        num_unassigned_electron: int = 0,
    ) -> "ElectronAllocation":
        if len(identities) != len(assigned_electrons):
            raise ValueError("identities and assigned_electrons must have the same length")
        nuclei = [Nucleus(x) for x in identities]
        assigned = {v: n for v, n in zip(nuclei, assigned_electrons)}
        allocs = []
        for i, j, ai, aj in allocated_connections:
            try:
                allocs.append(AllocatedConnection(nuclei[i], nuclei[j], ai, aj))
            except IndexError as exc:
                raise ValueError("temporary input index is out of range") from exc
        return cls(nuclei, assigned, allocs, num_unassigned_electron)

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
        return sum(self._assigned.values()) + self.num_unassigned_electron

    @property
    def num_total_electron(self) -> int:
        return self.num_electrons

    def to_constitution(self) -> Constitution:
        """Q_el(D)."""
        return Constitution(self.collection, self.connections, self.num_electrons)

    def show(self, **kwargs) -> None:
        """Display this electron-allocation form in a Jupyter/IPython notebook."""
        kwargs.setdefault("view", "allocation")
        _show_ccs(self, **kwargs)

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
        tokens.extend(Ele() for _ in range(self.num_unassigned_electron))
        for edge in self.allocated_connections:
            tokens.append(
                Con(
                    address[edge.endpoint1],
                    address[edge.endpoint2],
                    edge.allocation1,
                    edge.allocation2,
                )
            )
        return Word(tokens)

    def __repr__(self) -> str:
        counts = Counter(v.identity.label for v in self.nuclei)
        return (
            f"ElectronAllocation(nuclei={dict(sorted(counts.items()))}, "
            f"connections={len(self.allocated_connections)}, "
            f"N_e={self.num_electrons}, n_empty={self.num_unassigned_electron})"
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
    if a.num_unassigned_electron != b.num_unassigned_electron:
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
                self.form.num_unassigned_electron,
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
                self.form.num_unassigned_electron + 1,
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
                self.form.num_unassigned_electron,
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

    def show(self, **kwargs) -> None:
        """Replay and display this core construction word."""
        _show_ccs(self, **kwargs)

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

    def show(self, **kwargs) -> None:
        """Expand/replay and display this practical construction word."""
        _show_ccs(self, **kwargs)

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



# ---------------------------------------------------------------------------
# Optional notebook visualization (RDKit)
# ---------------------------------------------------------------------------

class CCSVisualizationError(RuntimeError):
    """Raised when an optional RDKit/IPython depiction cannot be produced."""


def _load_visualization_dependencies():
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError as exc:
        raise CCSVisualizationError(
            "RDKit is required for .show() visualization."
        ) from exc
    try:
        from IPython.display import SVG, display
    except ImportError as exc:
        raise CCSVisualizationError(
            "IPython is required for .show() notebook display."
        ) from exc
    return Chem, rdDepictor, rdMolDraw2D, SVG, display


def _display_order(nuclei, order=None):
    nuclei = tuple(nuclei)
    if order is not None:
        order = tuple(order)
        if len(order) != len(nuclei) or set(order) != set(nuclei):
            raise ValueError("order must enumerate every nucleus exactly once")
        return order

    # Direct intrinsic objects carry no indices.  This is only a display order.
    return tuple(sorted(
        nuclei,
        key=lambda v: (
            v.identity.atomic_number,
            v.identity.isotope or 0,
            v.identity.label,
            id(v),
        ),
    ))


def _show_source(obj, *, view="auto", order=None):
    if view not in {"auto", "constitution", "allocation"}:
        raise ValueError("view must be 'auto', 'constitution', or 'allocation'")

    if isinstance(obj, Word):
        state = obj.replay_state()
        form = state.form
        display_order = tuple(state.address_table) if order is None else _display_order(form.nuclei, order)
        if view == "constitution":
            return form.to_constitution(), display_order, "constitution"
        return form, display_order, "allocation"

    if isinstance(obj, PracticalWord):
        state = obj.expand().replay_state()
        form = state.form
        display_order = tuple(state.address_table) if order is None else _display_order(form.nuclei, order)
        if view == "constitution":
            return form.to_constitution(), display_order, "constitution"
        return form, display_order, "allocation"

    if isinstance(obj, ElectronAllocation):
        display_order = _display_order(obj.nuclei, order)
        if view == "constitution":
            return obj.to_constitution(), display_order, "constitution"
        return obj, display_order, "allocation"

    if isinstance(obj, Constitution):
        if view == "allocation":
            raise TypeError("Constitution has no electron-allocation data")
        return obj, _display_order(obj.nuclei, order), "constitution"

    raise TypeError(
        ".show() supports Constitution, ElectronAllocation, Word, and PracticalWord"
    )


def _auto_show_hydrogens(order, show_hydrogens):
    if show_hydrogens is not None:
        return bool(show_hydrogens)
    # Conventional compact default for molecules containing heavy atoms;
    # all-hydrogen species such as H2+ and H3+ remain visible.
    return not any(v.identity.symbol != "H" for v in order)


def _show_allocation_style(Chem, edge):
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
    # Unusual CCS allocations are kept literal instead of being interpreted
    # as a conventional bond order.
    return Chem.BondType.ZERO, None, None, False


def _show_allocation_label(edge):
    ai, aj = edge.allocation1, edge.allocation2
    if ai + aj == 0:
        return "(0,0)"
    if ai + aj == 1:
        return f"1e ({ai},{aj})"
    return f"({ai},{aj})"


def _build_show_mol(
    source,
    order,
    *,
    source_kind,
    show_hydrogens,
    show_indices,
    diagnostics,
    label_all_allocations,
):
    Chem, _, _, _, _ = _load_visualization_dependencies()

    visible = tuple(
        v for v in order
        if show_hydrogens or v.identity.symbol != "H"
    )
    visible_set = set(visible)
    display_index = {v: i for i, v in enumerate(order)}

    rw = Chem.RWMol()
    index = {}

    for nucleus in visible:
        atom = Chem.Atom(nucleus.identity.symbol)
        atom.SetNoImplicit(True)
        if nucleus.identity.isotope is not None:
            atom.SetIsotope(nucleus.identity.isotope)

        note = []
        if show_indices:
            note.append(str(display_index[nucleus]))

        if source_kind == "allocation":
            # q_v = Z(X(v)) - n_v: derived display quantity, not CCS identity.
            atom.SetFormalCharge(source.formal_charge(nucleus))
            if diagnostics:
                note.append(
                    f"n={source.assigned_electrons(nucleus)} "
                    f"r={source.residual_electrons(nucleus)} "
                    f"q={source.formal_charge(nucleus):+d}"
                )

        if note:
            atom.SetProp("atomNote", "  ".join(note))

        index[nucleus] = rw.AddAtom(atom)

    if source_kind == "constitution":
        for edge in source.connections:
            if edge.endpoint1 in visible_set and edge.endpoint2 in visible_set:
                # A solid line means only "direct CCS connection".
                rw.AddBond(
                    index[edge.endpoint1],
                    index[edge.endpoint2],
                    Chem.BondType.SINGLE,
                )
        return rw.GetMol()

    for edge in source.allocated_connections:
        if edge.endpoint1 not in visible_set or edge.endpoint2 not in visible_set:
            continue

        bond_type, donor, acceptor, conventional = _show_allocation_style(Chem, edge)
        if bond_type == Chem.BondType.DATIVE:
            begin_nucleus, end_nucleus = donor, acceptor
        else:
            begin_nucleus, end_nucleus = edge.endpoint1, edge.endpoint2

        rw.AddBond(index[begin_nucleus], index[end_nucleus], bond_type)
        bond = rw.GetBondBetweenAtoms(index[begin_nucleus], index[end_nucleus])

        if bond is not None and (label_all_allocations or not conventional):
            bond.SetProp("bondNote", _show_allocation_label(edge))

    return rw.GetMol()


def _show_metadata(source, source_kind, show_total_charge):
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


def _render_show_svg(
    obj,
    *,
    view="auto",
    order=None,
    show_hydrogens=None,
    show_indices=False,
    diagnostics=False,
    label_all_allocations=False,
    show_total_charge=False,
    width=1000,
    height=700,
    bond_length_scale=1.35,
    annotation_scale=0.8,
    footer_font_size=18,
):
    """Return the SVG used by .show(); primarily useful for tests."""
    if width <= 0 or height <= 0 or bond_length_scale <= 0:
        raise ValueError("width, height, and bond_length_scale must be positive")

    Chem, rdDepictor, rdMolDraw2D, _, _ = _load_visualization_dependencies()
    source, display_order, source_kind = _show_source(
        obj, view=view, order=order
    )
    show_hydrogens = _auto_show_hydrogens(display_order, show_hydrogens)

    mol = _build_show_mol(
        source,
        display_order,
        source_kind=source_kind,
        show_hydrogens=show_hydrogens,
        show_indices=show_indices,
        diagnostics=diagnostics,
        label_all_allocations=label_all_allocations,
    )

    try:
        mol.UpdatePropertyCache(strict=False)
    except Exception:
        pass

    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)

    if bond_length_scale != 1:
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(
                i,
                (p.x * bond_length_scale, p.y * bond_length_scale, 0.0),
            )

    metadata = _show_metadata(source, source_kind, show_total_charge)
    footer_height = max(40, 8 + len(metadata) * (footer_font_size + 3))
    mol_height = max(120, height - footer_height)

    drawer = rdMolDraw2D.MolDraw2DSVG(width, mol_height)
    opts = drawer.drawOptions()
    opts.annotationFontScale = annotation_scale
    opts.bondLineWidth = 2.2
    opts.additionalAtomLabelPadding = 0.12
    opts.padding = 0.05
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    # Extend the canvas and add N_e / optional Q / n_empty at bottom right.
    import html
    import re as _re
    svg = _re.sub(
        rf"height=['\"]{mol_height}px['\"]",
        f"height='{height}px'",
        svg,
        count=1,
    )
    svg = _re.sub(
        rf"viewBox=['\"]0 0 {width} {mol_height}['\"]",
        f"viewBox='0 0 {width} {height}'",
        svg,
        count=1,
    )

    y = mol_height + footer_font_size + 2
    footer = [
        f"<rect x='0' y='{mol_height}' width='{width}' height='{footer_height}' fill='white'/>"
    ]
    for line in metadata:
        footer.append(
            f"<text x='{width - 14}' y='{y}' text-anchor='end' "
            f"font-family='sans-serif' font-size='{footer_font_size}px'>"
            f"{html.escape(line)}</text>"
        )
        y += footer_font_size + 3

    return svg.replace("</svg>", "\n" + "\n".join(footer) + "\n</svg>")


def _show_ccs(obj, **kwargs) -> None:
    """Notebook display backend used by CCS object .show() methods."""
    _, _, _, SVG, display = _load_visualization_dependencies()
    display(SVG(data=_render_show_svg(obj, **kwargs)))


__all__ = [
    # S1
    "NuclearIdentity",
    "Nucleus",
    "NuclearCollection",
    "Connection",
    "AllocatedConnection",
    "Constitution",
    "ElectronAllocation",
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