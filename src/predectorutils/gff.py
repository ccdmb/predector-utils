#!/usr/bin/env python3

from typing import Optional, Union
from typing import Set, List, Dict
from typing import Sequence, Mapping
from typing import Iterable, Iterator
from typing import TypeVar
from typing import cast
from typing import TextIO

from enum import Enum
from collections import deque

from predectorutils.higher import fmap
from predectorutils.parsers import (
    ParseError,
    LineParseError,
    FieldParseError,
    raise_it,
    parse_field,
    parse_int,
    parse_float,
    parse_or_none,
    is_one_of,
    parse_str,
    parse_bool_options,
)

T = TypeVar('T')

GFF3_KEY_TO_ATTR: Dict[str, str] = {
    "ID": "id",
    "Name": "name",
    "Alias": "alias",
    "Parent": "parent",
    "Target": "target",
    "Gap": "gap",
    "Derives_from": "derives_from",
    "Note": "note",
    "Dbxref": "dbxref",
    "Ontology_term": "ontology_term",
    "Is_circular": "is_circular",
}

GFF3_ATTR_TO_KEY: Dict[str, str] = {v: k for k, v in GFF3_KEY_TO_ATTR.items()}

GFF3_WRITE_ORDER: List[str] = [
    "ID",
    "Name",
    "Alias",
    "Parent",
    "Target",
    "Gap",
    "Derives_from",
    "Note",
    "Dbxref",
    "Ontology_term",
    "Is_circular",
]

rec_seqid = raise_it(parse_field(parse_str, "seqid"))
rec_source = raise_it(parse_field(parse_str, "source"))
rec_type = raise_it(parse_field(parse_str, "type"))
rec_start = raise_it(parse_field(parse_int, "start"))
rec_end = raise_it(parse_field(parse_int, "end"))
rec_score = raise_it(parse_field(parse_or_none(parse_float, "."), "score"))
rec_strand = raise_it(parse_field(is_one_of(["-", "+", ".", "?"]), "strand"))
rec_phase = raise_it(parse_field(is_one_of(["0", "1", "2", "."]), "phase"))


def parse_attr_list(string: str) -> List[str]:
    return list(f.strip() for f in string.strip(", ").split(","))


attr_is_circular = raise_it(parse_field(
    parse_bool_options(["true", "TRUE", "True"],
                       ["false", "FALSE", "False"]),
    "is_circular",
    "attributes"
))

attr_target_id = raise_it(parse_field(parse_str, "target.id" "attributes"))
attr_target_start = raise_it(parse_field(
    parse_int,
    "target.start",
    "attributes"
))
attr_target_end = raise_it(parse_field(parse_int, "target.end", "attributes"))
attr_target_strand = raise_it(parse_field(
    is_one_of(["+", "-"]),
    "target.strand",
    "attributes"
))
attr_gap_code = raise_it(parse_field(
    is_one_of(["M", "I", "D", "F", "R"]),
    "gap",
    "attributes"
))
attr_gap_len = raise_it(parse_field(parse_int, "gap", "attributes"))


def flatten_list_of_lists(li: Iterable[Iterable[T]]) -> Iterator[T]:
    for i in li:
        for j in i:
            yield j
    return


def attr_escape(string: str) -> str:
    return (
        string
        .replace("%", "%25")
        .replace(";", "%3B")
        .replace(",", "%2C")
        .replace("=", "%3D")
        .replace("\t", "%09")
        .replace("\n", "%0A")
        .replace("\r", "%0D")
    )


def attr_unescape(string: str) -> str:
    return (
        string
        .replace("%3B", ";")
        .replace("%2C", ",")
        .replace("%3D", "=")
        .replace("%09", "\t")
        .replace("%0A", "\n")
        .replace("%0D", "\r")
        .replace("%25", "%")
    )


class Strand(Enum):
    PLUS = 0
    MINUS = 1
    UNSTRANDED = 2
    UNKNOWN = 3

    def __str__(self):
        into_str_map: List[str] = ["+", "-", ".", "?"]
        return into_str_map[self.value]

    def __repr__(self):
        return f"Strand.{self.name}"

    @classmethod
    def parse(cls, string: str) -> "Strand":
        from_str_map: Dict[str, Strand] = {
            "+": cls.PLUS,
            "-": cls.MINUS,
            ".": cls.UNSTRANDED,
            "?": cls.UNKNOWN,
        }

        try:
            return from_str_map[string]
        except KeyError:
            valid = list(from_str_map.keys())
            raise ValueError(f"Invalid option. Must be one of {valid}")


class Phase(Enum):
    FIRST = 0
    SECOND = 1
    THIRD = 2
    NOT_CDS = 3

    def __str__(self):
        into_str_map: List[str] = ["0", "1", "2", "."]
        return into_str_map[self.value]

    def __repr__(self):
        return f"Phase.{self.name}"

    @classmethod
    def parse(cls, string: str) -> "Phase":
        from_str_map: Dict[str, Phase] = {
            "0": cls.FIRST,
            "1": cls.SECOND,
            "2": cls.THIRD,
            ".": cls.NOT_CDS,
        }

        try:
            return from_str_map[string]
        except KeyError:
            valid = list(from_str_map.keys())
            raise ValueError(f"Invalid option. Must be one of {valid}")


class TargetStrand(Enum):
    PLUS = 0
    MINUS = 1

    def __str__(self) -> str:
        into_str_map = ["+", "-"]
        return into_str_map[self.value]

    def __repr__(self) -> str:
        return "TargetStrand.{self.name}"

    @classmethod
    def parse(cls, string: str) -> "TargetStrand":
        from_str_map = {
            "+": cls.PLUS,
            "-": cls.MINUS,
        }

        try:
            return from_str_map[string]
        except KeyError:
            valid = list(from_str_map.keys())
            raise ValueError(f"Invalid option. Must be one of {valid}")


class Target(object):

    """
    Indicates the target of a nucleotide-to-nucleotide or
    protein-to-nucleotide alignment.
    The format of the value is "target_id start end [strand]",
    where strand is optional and may be "+" or "-".
    If the target_id contains spaces, they must be escaped as hex escape %20.
    """

    def __init__(
        self,
        target_id: str,
        start: int,
        end: int,
        strand: Optional[TargetStrand] = None,
    ) -> None:
        self.target_id = target_id
        self.start = start
        self.end = end
        self.strand = strand
        return

    def __repr__(self) -> str:
        if self.strand is None:
            return f"Target('{self.target_id}', {self.start}, {self.end})"
        else:
            strand = repr(self.strand)
            return (
                f"Target('{self.target_id}', {self.start}, "
                f"{self.end}, {strand})"
            )

    def __str__(self) -> str:
        target_id = attr_escape(self.target_id).replace(" ", "%20")

        # Recode back to 1 based (inclusive)
        start = self.start + 1

        if self.strand is None:
            return "{} {} {}".format(target_id, start, self.end)
        else:
            return "{} {} {} {}".format(
                target_id,
                start,
                self.end,
                self.strand
            )

    @classmethod
    def parse(cls, string: str, unescape: bool = True) -> "Target":
        split_string = string.strip().split(" ")

        if len(split_string) < 3:
            raise FieldParseError(
                "target",
                ("Too few fields in Target string. Expected 3 or 4, "
                 f"but got {len(split_string)}."),
                "attribute"
            )
        elif len(split_string) > 4:
            raise FieldParseError(
                "target",
                ("Too many fields in Target strand. Expected 3 or 4, "
                 f"but got {len(split_string)}."),
                "attribute"
            )

        if unescape:
            target_id = attr_target_id(
                attr_unescape(split_string[0])
                .replace("%20", " ")
            )
        else:
            target_id = attr_target_id(split_string[0])
        # We want 0 based exclusive
        start = attr_target_start(split_string[1]) - 1
        end = attr_target_end(split_string[2])

        if len(split_string) == 3:
            return cls(target_id, start, end)
        else:
            strand = TargetStrand.parse(attr_target_strand(split_string[3]))
            return cls(target_id, start, end, strand)


class GapCode(Enum):

    MATCH = 0
    INSERT = 1
    DELETE = 2
    FORWARD = 3
    REVERSE = 4

    def __str__(self) -> str:
        into_str_map = ["M", "I", "D", "F", "R"]
        return into_str_map[self.value]

    def __repr__(self) -> str:
        return f"GapCode.{self.name}"

    @classmethod
    def parse(cls, string: str) -> "GapCode":
        from_str_map = {
            "M": cls.MATCH,
            "I": cls.INSERT,
            "D": cls.DELETE,
            "F": cls.FORWARD,
            "R": cls.REVERSE,
        }

        try:
            return from_str_map[string]
        except KeyError:
            valid = list(from_str_map.keys())
            raise ValueError(f"Invalid option. Must be one of {valid}")
        return


class GapElement(object):

    def __init__(self, code: GapCode, length: int) -> None:
        """ """
        assert length >= 0

        self.code = code
        self.length = length
        return

    def __str__(self) -> str:
        return "{}{}".format(self.code, self.length)

    def __repr__(self) -> str:
        code = repr(self.code)
        return f"GapElement({code}, {self.length})"

    @classmethod
    def parse(cls, string: str) -> "GapElement":
        string = string.strip()
        code = GapCode.parse(attr_gap_code(string[:1]))
        length = attr_gap_len(string[1:])
        return cls(code, length)


class Gap(object):

    def __init__(self, elements: Sequence[GapElement]) -> None:
        self.elements = list(elements)
        return

    def __str__(self) -> str:
        return " ".join(map(str, self.elements))

    def __repr__(self) -> str:
        elems = repr(list(self.elements))
        return f"Gap({elems})"

    @classmethod
    def parse(cls, string: str) -> "Gap":
        split_string = string.strip().split(" ")
        elements = [GapElement.parse(s) for s in split_string if s != '']
        return cls(elements)


class GFFRecord(object):

    columns: List[str] = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes"
    ]

    def __init__(
        self,
        seqid: str,
        source: str,
        type: str,
        start: int,
        end: int,
        score: Optional[float] = None,
        strand: Strand = Strand.UNSTRANDED,
        phase: Phase = Phase.NOT_CDS,
        attributes: Optional["GFFAttributes"] = None,
        parents: Optional[Sequence["GFFRecord"]] = None,
        children: Optional[Sequence["GFFRecord"]] = None,
        derives_froms: Optional[Sequence["GFFRecord"]] = None,
        derivatives: Optional[Sequence["GFFRecord"]] = None,
    ) -> None:
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase

        if attributes is None:
            self.attributes = GFFAttributes()
        else:
            self.attributes = attributes

        self.parents: List[GFFRecord] = []
        if parents is not None:
            self.add_parents(parents)

        self.children: List[GFFRecord] = []
        if children is not None:
            self.add_children(children)

        self.derives_from: List[GFFRecord] = []
        if derives_froms is not None:
            self.add_derives_froms(derives_froms)

        self.derivatives: List[GFFRecord] = []
        if derivatives is not None:
            self.add_derivatives(derivatives)
        return

    def __str__(self) -> str:
        return self.as_str()

    def as_str(self, escape: bool = True) -> str:
        values = []
        for name in self.columns:
            value = getattr(self, name)

            if value is None:
                values.append(".")
            # Convert back to 1-based inclusive indexing.
            elif name == "start":
                values.append(str(value + 1))
            elif name == "attributes" and not escape:
                values.append(value.as_str(escape=False))
            else:
                values.append(str(value))

        return "\t".join(values)

    def __repr__(self) -> str:
        parameters = []
        for col in self.columns:
            val = repr(getattr(self, col))
            parameters.append(f"{val}")

        joined_parameters = ", ".join(parameters)
        return f"GFFRecord({joined_parameters})"

    def length(self) -> int:
        """ Returns the distance between start and end. """
        return self.end - self.start

    def __len__(self) -> int:
        return self.length()

    def add_child(self, child: "GFFRecord") -> None:
        if child not in self.children:
            self.children.append(child)
        if self not in child.parents:
            child.parents.append(self)
        return

    def add_parent(self, parent: "GFFRecord") -> None:
        if parent not in self.parents:
            self.parents.append(parent)

        if self not in parent.children:
            parent.children.append(self)
        return

    def add_children(self, children: Sequence["GFFRecord"]) -> None:
        for child in children:
            self.add_child(child)
        return

    def add_parents(self, parents: Sequence["GFFRecord"]) -> None:
        for parent in parents:
            self.add_parent(parent)
        return

    def add_derivative(self, derivative: "GFFRecord") -> None:
        if derivative not in self.derivatives:
            self.derivatives.append(derivative)
        if self not in derivative.parents:
            derivative.derives_from.append(self)
        return

    def add_derives_from(self, derives_from: "GFFRecord") -> None:
        if derives_from not in self.derives_from:
            self.derives_from.append(derives_from)

        if self not in derives_from.derivatives:
            derives_from.derivatives.append(self)
        return

    def add_derivatives(self, derivatives: Sequence["GFFRecord"]) -> None:
        for derivative in derivatives:
            self.add_derivative(derivative)
        return

    def add_derives_froms(self, derives_froms: Sequence["GFFRecord"]) -> None:
        for derives_from in derives_froms:
            self.add_derives_from(derives_from)
        return

    def update_parents(self) -> None:
        parent_ids = []
        for parent in self.parents:
            parent_id = parent.attributes.id
            assert parent_id is not None
            parent_ids.append(parent_id)

        self.attributes.parent = parent_ids
        return

    def update_derives_from(self) -> None:
        df_ids = []
        for df in self.derives_from:
            df_id = df.attributes.id
            assert df_id is not None
            df_ids.append(df_id)

        self.attributes.derives_from = df_ids
        return

    def traverse_children(
        self,
        sort: bool = False,
        breadth: bool = False,
    ) -> Iterator["GFFRecord"]:
        """ A graph traversal of this `GFFRecord`s children.

        Keyword arguments:
        sort -- Sort the children so that they appear in increasing order.
        breadth -- Perform a breadth first search rather than the default
            depth first search.

        Returns:
        A generator yielding `GFFRecord`s.
        """

        should_reverse = not breadth

        seen: Set[GFFRecord] = set()

        to_visit = deque([self])

        while len(to_visit) > 0:

            if breadth:
                node = to_visit.popleft()
            else:
                node = to_visit.pop()

            # NB uses id() for hashing
            if node in seen:
                continue
            else:
                seen.add(node)

            children = list(node.children)
            if sort:
                children.sort(
                    key=lambda f: (f.seqid, f.start, f.end, f.type),
                    reverse=should_reverse
                )

            to_visit.extend(children)
            yield node

        return None

    def traverse_parents(
        self,
        sort: bool = False,
        breadth: bool = False,
    ) -> Iterator["GFFRecord"]:
        """ A graph traversal of this `GFFRecord`s parents.

        Keyword arguments:
        sort -- Sort the parents so that they appear in increasing order.
        breadth -- Perform a breadth first search rather than the default
            depth first search.

        Returns:
        A generator yielding `GFFRecord`s.
        """

        should_reverse = not breadth

        seen: Set[GFFRecord] = set()
        to_visit = deque([self])

        while len(to_visit) > 0:

            if breadth:
                node = to_visit.popleft()
            else:
                node = to_visit.pop()

            if node in seen:
                continue
            else:
                seen.add(node)

            parents = list(node.parents)
            if sort:
                parents.sort(
                    key=lambda f: (f.seqid, f.start, f.end, f.type),
                    reverse=should_reverse
                )

            to_visit.extend(parents)
            yield node

        return None

    @classmethod
    def parse(
        cls,
        string: str,
        strip_quote: bool = False,
        unescape: bool = True,
    ) -> "GFFRecord":
        """ Parse a gff line string as a `GFFRecord`.

        Keyword arguments:
        string -- The gff line to parse.
        format -- What format the gff file is in.
            Currently only GFF3 is supported.
        strip_quote -- Strip quotes from attributes values. The specification
            says that they should not be stripped, so we don't by default.
        unescape -- Unescape reserved characters in the attributes to their
            original values. I.E. some commas, semicolons, newlines etc.

        Returns:
        A `GFFRecord`
        """

        sline = string.strip().split("\t")
        sline_len = len(sline)
        columns_len = len(cls.columns)

        if sline_len < (columns_len - 1):
            raise ValueError((
                "Line has too few columns. "
                f"Expected: {columns_len}, Encountered: {sline_len}"
            ))

        fields: Dict[str, str] = dict(zip(cls.columns, sline))
        if sline_len == columns_len - 1:
            fields["attributes"] = ""

        # 0-based indexing exclusive
        start = rec_start(fields["start"]) - 1
        end = rec_end(fields["end"])

        if start > end:
            tmp = start
            start = end
            end = tmp
            del tmp

        score = rec_score(fields["score"])
        strand = Strand.parse(rec_strand(fields["strand"]))
        phase = Phase.parse(rec_phase(fields["phase"]))

        attributes = GFFAttributes.parse(
            fields["attributes"],
            strip_quote=strip_quote,
            unescape=unescape,
        )

        return cls(
            rec_seqid(fields["seqid"]),
            rec_source(fields["source"]),
            rec_type(fields["type"]),
            start,
            end,
            score,
            strand,
            phase,
            attributes
        )

    def trim_ends(self, length: int) -> None:
        """ Trim a specified number of bases from each end of the feature. """

        from math import ceil

        if self.length() <= 2:
            length = 0
        elif self.length() < (2 * length):
            length = ceil(self.length() / 4)

        self.start += length
        self.end -= length
        return

    def pad_ends(self, length: int, max_value: Optional[int] = None) -> None:

        if (self.start - length) < 0:
            self.start = 0
        else:
            self.start -= length

        if max_value is not None and (self.end + length) > max_value:
            self.end = max_value
        else:
            self.end += length
        return

    def expand_to_children(self) -> None:
        if len(self.children) == 0:
            return

        min_ = min(c.start for c in self.children)
        max_ = max(c.end for c in self.children)

        if min_ < self.start:
            self.start = min_

        if max_ > self.end:
            self.end = max_
        return

    def shrink_to_children(self) -> None:
        if len(self.children) == 0:
            return

        min_ = min(c.start for c in self.children)
        max_ = max(c.end for c in self.children)

        if min_ > self.start:
            self.start = min_

        if max_ < self.end:
            self.end = max_
        return

    def copy(self) -> "GFFRecord":
        """ You'll still need to update the ID """
        from copy import copy, deepcopy

        node_copy = copy(self)
        node_copy.attributes = deepcopy(self.attributes)

        if node_copy.attributes is not None:
            node_copy.attributes.parent = []

        node_copy.children = []
        node_copy.parents = []
        return node_copy

    @classmethod  # noqa
    def from_file(
        cls,
        handle: TextIO,
        strip_quote: bool = False,
        unescape: bool = True,
    ) -> Iterator["GFFRecord"]:
        """
        Yes yes I need to make this more modular... not time right now.
        """

        from collections import defaultdict

        id_to_record = defaultdict(list)

        to_yield = list()
        lonely_children = set()
        lonely_derivatives = set()

        # Children that are encountered before their parents
        undefined_parents: Dict[str, List[GFFRecord]] = defaultdict(list)

        # derivatives encountered before their parents.
        undefined_derives_froms: Dict[str, List[GFFRecord]] = defaultdict(list)

        for i, line in enumerate(handle):
            if line.startswith("#"):
                continue
            elif line.strip() == "":
                continue

            try:
                record = GFFRecord.parse(line, strip_quote, unescape)
            except (LineParseError, FieldParseError) as e:
                raise e.as_parse_error(line=i).add_filename_from_handle(handle)

            id_ = record.attributes.id

            if id_ is not None:
                id_to_record[id_].append(record)

                uparents = undefined_parents.pop(id_, [])
                uderives = undefined_derives_froms.pop(id_, [])
                record.add_children(uparents)
                record.add_derivatives(uderives)

                for r in uparents:
                    lonely_children.discard(r)
                    if r not in lonely_derivatives:
                        to_yield.append(r)

                for r in uderives:
                    lonely_derivatives.discard(r)
                    if r not in lonely_children:
                        to_yield.append(r)

            is_missing_parent = False
            for parent in record.attributes.parent:
                if parent not in id_to_record:
                    undefined_parents[parent].append(record)
                    lonely_children.add(record)
                    is_missing_parent = True
                else:
                    record.add_parents(id_to_record.get(parent, []))

            is_missing_derives_from = False
            for derives_from in record.attributes.derives_from:
                if derives_from not in id_to_record:
                    undefined_derives_froms[derives_from].append(record)
                    lonely_derivatives.add(record)
                    is_missing_derives_from = True
                else:
                    record.add_derives_froms(
                        id_to_record.get(derives_from, []))

            if not (is_missing_parent or is_missing_derives_from):
                to_yield.append(record)

            while len(to_yield) > 0:
                yield to_yield.pop()

        if (len(undefined_parents) > 0) or (len(undefined_derives_froms) > 0):
            upar = set(undefined_parents.keys())
            uder = set(undefined_derives_froms.keys())

            message = [
                "Reached the end of GFF file with undefined references."
            ]

            if len(upar) > 0:
                message.append(f"Expected to find parent: {repr(upar)}.")

            if len(uder) > 0:
                message.append(f"Expected to find derives_from: {repr(uder)}.")

            raise (
                ParseError(filename=None, line=i, message=" ".join(message))
                .add_filename_from_handle(handle)
            )
        return


class GFFAttributes(object):

    def __init__(
        self,
        id: Optional[str] = None,
        name: Optional[str] = None,
        alias: Optional[Sequence[str]] = None,
        parent: Optional[Sequence[str]] = None,
        target: Optional[Target] = None,
        gap: Optional[Gap] = None,
        derives_from: Optional[Sequence[str]] = None,
        note: Optional[Sequence[str]] = None,
        dbxref: Optional[Sequence[str]] = None,
        ontology_term: Optional[Sequence[str]] = None,
        is_circular: bool = False,
        custom: Optional[Mapping[str, Sequence[str]]] = None,
    ) -> None:
        self.id = id
        self.name = name

        if alias is not None:
            self.alias: List[str] = list(alias)
        else:
            self.alias = []

        if parent is not None:
            self.parent: List[str] = list(parent)
        else:
            self.parent = []

        self.target = target
        self.gap = gap

        if derives_from is not None:
            self.derives_from: List[str] = list(derives_from)
        else:
            self.derives_from = []

        if note is not None:
            self.note: List[str] = list(note)
        else:
            self.note = []

        if dbxref is not None:
            self.dbxref: List[str] = list(dbxref)
        else:
            self.dbxref = []

        if ontology_term is not None:
            self.ontology_term: List[str] = list(ontology_term)
        else:
            self.ontology_term = []

        self.is_circular = is_circular

        self.custom: Dict[str, Union[str, List[str]]] = {}
        if custom is not None:
            for k, v in custom.items():
                if isinstance(v, str):
                    self.custom[k] = v
                else:
                    self.custom[k] = list(v)

        return

    @classmethod
    def _parse_list_of_strings(
        cls,
        value: str,
        strip_quote: bool = False,
        unescape: bool = True
    ) -> List[str]:
        """ Parses a gff attribute list of strings. """
        if value == "":
            return []

        if strip_quote:
            value = value.strip("\"' ")

        if unescape:
            return [attr_unescape(v) for v in parse_attr_list(value)]
        else:
            return parse_attr_list(value)

    @classmethod
    def parse(
        cls,
        string: str,
        strip_quote: bool = False,
        unescape: bool = True,
    ) -> "GFFAttributes":
        if string.strip() in (".", ""):
            return cls()

        fields = (
            f.split("=", maxsplit=1)
            for f
            in string.strip(" ;").split(";")
        )

        if strip_quote:
            kvpairs: Dict[str, str] = {
                k.strip(): v.strip(" '\"")
                for k, v
                in fields
            }
        else:
            kvpairs = {
                k.strip(): v.strip()
                for k, v
                in fields
            }

        if unescape:
            id = fmap(attr_unescape, kvpairs.pop("ID", None))
        else:
            id = kvpairs.pop("ID", None)

        if id == "":
            id = None

        if unescape:
            name = fmap(attr_unescape, kvpairs.pop("Name", None))
        else:
            name = kvpairs.pop("Name", None)

        if name == "":
            name = None

        alias = cls._parse_list_of_strings(
            kvpairs.pop("Alias", ""),
            strip_quote,
            unescape
        )

        parent = cls._parse_list_of_strings(
            kvpairs.pop("Parent", ""),
            strip_quote,
            unescape
        )

        target: Optional[Target] = fmap(
            lambda x: Target.parse(x, unescape),
            kvpairs.pop("Target", None)
        )

        gap = fmap(Gap.parse, kvpairs.pop("Gap", None))

        derives_from = cls._parse_list_of_strings(
            kvpairs.pop("Derives_from", ""),
            strip_quote,
            unescape
        )

        note = cls._parse_list_of_strings(
            kvpairs.pop("Note", ""),
            strip_quote,
            unescape
        )

        dbxref = cls._parse_list_of_strings(
            kvpairs.pop("Dbxref", ""),
            strip_quote,
            unescape
        )

        ontology_term = cls._parse_list_of_strings(
            kvpairs.pop("Ontology_term", ""),
            strip_quote,
            unescape
        )

        is_circular = attr_is_circular(kvpairs.pop("Is_circular", "false"))

        custom: Dict[str, Union[str, List[str]]] = dict()
        for k, v in kvpairs.items():
            if "," in v:
                custom[k] = cls._parse_list_of_strings(
                    v, strip_quote, unescape)
            elif v != "":
                custom[k] = v

        return cls(
            id,
            name,
            alias,
            parent,
            target,
            gap,
            derives_from,
            note,
            dbxref,
            ontology_term,
            is_circular,
            custom
        )

    def is_empty(self) -> bool:  # noqa
        # Yes, this could be written as single boolean comparison.
        # But it gets so long that it's hard to understand.
        if len(self.custom) > 0:
            return False
        elif self.id is not None:
            return False
        elif self.name is not None:
            return False
        elif len(self.alias) > 0:
            return False
        elif len(self.parent) > 0:
            return False
        elif self.target is not None:
            return False
        elif self.gap is not None:
            return False
        elif len(self.derives_from) > 0:
            return False
        elif len(self.note) > 0:
            return False
        elif len(self.dbxref) > 0:
            return False
        elif len(self.ontology_term) > 0:
            return False
        elif self.is_circular:
            return False
        else:
            return True

    def __str__(self) -> str:
        return self.as_str()

    def as_str(self, escape: bool = True) -> str:
        # Avoid having an empty string at the end.
        if self.is_empty():
            return "."

        keys = []
        keys.extend(GFF3_WRITE_ORDER)
        keys.extend(self.custom.keys())

        kvpairs = []
        for key in keys:
            value = self[key]
            if value is None or value == []:
                continue
            elif key == "Is_circular" and not value:
                continue

            if escape:
                key = attr_escape(key)

            if isinstance(value, list) and escape:
                value = ",".join(attr_escape(str(v)) for v in value)

            elif isinstance(value, list):
                value = ",".join(str(v) for v in value)

            elif isinstance(value, bool):
                value = "true" if value else "false"

            else:
                if escape:
                    value = attr_escape(str(value))
                else:
                    value = str(value)

            kvpairs.append((key, value))

        return ";".join(f"{k}={v}" for k, v in kvpairs)

    def __repr__(self) -> str:
        param_names = [GFF3_KEY_TO_ATTR[k] for k in GFF3_WRITE_ORDER]
        param_names.append("custom")

        parameters = []
        for param in param_names:
            value = getattr(self, param)

            if value is None or value == []:
                continue

            if isinstance(value, list):
                value = repr(list(value))
            else:
                value = repr(value)

            parameters.append(f"{param}={value}")

        joined_parameters = ", ".join(parameters)
        return f"GFF3Attributes({joined_parameters})"

    def __getitem__(
        self,
        key: str,
    ) -> Union[str, Sequence[str], Target, Gap, bool, None]:
        if key in GFF3_KEY_TO_ATTR:
            return getattr(self, GFF3_KEY_TO_ATTR[key])
        else:
            return self.custom[key]

    def __setitem__(
        self,
        key: str,
        value: Union[str, Sequence[str], Target, Gap, bool, None],
    ) -> None:
        """ If the key is an attr set it, otherwise update the custom dict."""

        if key in GFF3_KEY_TO_ATTR:
            setattr(self, GFF3_KEY_TO_ATTR[key], value)
        else:
            self.custom[key] = cast(str, value)  # Cast is for mypy
        return

    def __delitem__(self, key: str) -> None:
        """ Removes an item from the custom dictionary or resets
        attribute to default """

        if key in ("ID", "Name", "Target", "Gap"):
            setattr(self, GFF3_KEY_TO_ATTR[key], None)

        elif key in ("Alias", "Parent", "Derives_from", "Note",
                     "Dbxref", "Ontology_term"):
            setattr(self, GFF3_KEY_TO_ATTR[key], [])

        elif key == "Is_circular":
            setattr(self, GFF3_KEY_TO_ATTR[key], False)

        else:
            del self.custom[key]

        return

    def get(
        self,
        key: str,
        default: Union[str, Sequence[str], Target, Gap, bool, None] = None,
    ) -> Union[str, Sequence[str], Target, Gap, bool, None]:
        """ Gets an atrribute or element from the custom dict. """

        if key in GFF3_KEY_TO_ATTR:
            return getattr(self, GFF3_KEY_TO_ATTR[key])
        else:
            return self.custom.get(key, default)

    def pop(
        self,
        key: str,
        default: Union[str, Sequence[str], Target, Gap, bool, None] = None,
    ) -> Union[str, Sequence[str], Target, Gap, bool, None]:
        """ Removes an item from the attributes and returns its value.

        If the item is one of the reserved GFF3 terms, the
        value is reset to the default.
        """

        if key in GFF3_KEY_TO_ATTR:
            value = getattr(self, GFF3_KEY_TO_ATTR[key])
            del self[GFF3_KEY_TO_ATTR[key]]
            return value

        else:
            return self.custom.pop(key, default)
