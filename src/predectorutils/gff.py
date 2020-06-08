#!/usr/bin/env python3

import logging
from typing import cast
from typing import Optional, Union
from typing import Set, List, Dict
from typing import Sequence, Mapping
from typing import Iterator
# from typing import Hashable, Any
# from typing import TextIO
# from typing import Callable
from typing import TypeVar

from enum import Enum
from collections import deque

from predectorutils.higher import fmap
from predectorutils.parsers import (
    parse_float,
    parse_str,
    parse_int,
    parse_or_none,
    is_one_of,
)

logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)

T = TypeVar("T")


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
        target_id = self.target_id.replace(" ", "%20")

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
    def parse(cls, string: str) -> "Target":
        split_string = string.strip().split(" ")

        if len(split_string) < 3:
            raise ValueError("Too few fields in Target string.")

        elif len(split_string) == 3:
            target_id = split_string[0]
            start = int(split_string[1]) - 1  # We want 0 based exclusive
            end = int(split_string[2])
            return cls(target_id, start, end)

        elif len(split_string) == 4:
            target_id = split_string[0]
            start = int(split_string[1]) - 1  # We want 0 based exclusive
            end = int(split_string[2])
            strand = TargetStrand.parse(split_string[3])
            return cls(target_id, start, end, strand)

        else:
            raise ValueError("Too many fields in Target string.")


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
        code = GapCode.parse(string[:1])
        length = int(string[1:])
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


def parse_attr_list(string: str) -> List[str]:
    return list(f.strip() for f in string.strip(", ").split(","))


def parse_bool(string: str) -> bool:
    if string in ("true", "True", "TRUE"):
        return True
    elif string in ("false", "False", "FALSE"):
        return False
    else:
        raise ValueError(f"Invalid boolean encountered: {string}")


class GFF3Attributes(object):

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
        is_circular: Optional[bool] = None,
        custom: Optional[Mapping[str, str]] = None,
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

        if custom is not None:
            custom_not_none: Dict[str, str] = dict(custom)
        else:
            custom_not_none = {}

        self.custom: Dict[str, str] = custom_not_none
        return

    @classmethod
    def _parse_list_of_strings(
        cls,
        value: str,
        strip_quote: bool = False,
        unescape: bool = False
    ) -> List[str]:
        """ Parses a gff attribute list of strings. """
        if value == "":
            return []

        if strip_quote:
            value = value.strip("\"' ")

        if unescape:
            return [cls._attr_unescape(v) for v in parse_attr_list(value)]
        else:
            return parse_attr_list(value)

    @classmethod
    def parse(
        cls,
        string: str,
        strip_quote: bool = False,
        unescape: bool = False,
    ) -> "GFF3Attributes":
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

        id = kvpairs.pop("ID", None)
        name = kvpairs.pop("Name", None)

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
            Target.parse,
            fmap(
                lambda x: cls._attr_unescape(x) if unescape else x,
                kvpairs.pop("Target", None)
            )
        )

        gap = fmap(
            Gap.parse,
            fmap(
                lambda x: cls._attr_unescape(x) if unescape else x,
                kvpairs.pop("Gap", None)
            )
        )

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

        is_circular = fmap(parse_bool, kvpairs.pop("Is_circular", None))

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
            kvpairs
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
        # Nothing will be printed is this is false, so we have to check for
        # both
        elif (self.is_circular is not None) and (self.is_circular):
            return False
        else:
            return True

        return ((self.gene_id is None) and
                (self.transcript_id is None) and
                (len(self.custom) == 0))

    def __str__(self) -> str:
        return self.as_str(escape=True)

    def as_str(self, escape: bool = False) -> str:
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

            if escape:
                key = self._attr_escape(key)

            if isinstance(value, list) and escape:
                value = ",".join(self._attr_escape(str(v)) for v in value)

            elif isinstance(value, list):
                value = ",".join(str(v) for v in value)

            elif isinstance(value, bool):
                value = "true" if value else "false"

            else:
                if escape:
                    value = self._attr_escape(str(value))
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

    @staticmethod
    def _attr_escape(string: str) -> str:
        return (
            string
            .replace(";", "%3B")
            .replace(",", "%2C")
            .replace("=", "%3D")
            .replace("\t", "%09")
            .replace("\n", "%0A")
            .replace("\r", "%0D")
        )

    @staticmethod
    def _attr_unescape(string: str) -> str:
        return (
            string
            .replace("%3B", ";")
            .replace("%2C", ",")
            .replace("%3D", "=")
            .replace("%09", "\t")
            .replace("%0A", "\n")
            .replace("%0D", "\r")
        )

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

    @classmethod
    def infer_from_many(cls, strands: Sequence["Strand"]) -> "Strand":
        """ Infer the strand from a sequence of strands. """

        if len(set(strands)) > 1:
            logger.warning("Ambiguous strand found during parent inference.")
            return Strand.UNKNOWN
        else:
            return list(strands)[0]


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


class GFF3Record(object):

    columns: List[str] = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "id",
        "parents",
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
        attributes: Optional[GFF3Attributes] = None,
        parents: Optional[Sequence["GFF3Record"]] = None,
        children: Optional[Sequence["GFF3Record"]] = None,
    ) -> None:
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes

        self.parents: List[GFF3Record] = []
        if parents is not None:
            self.add_parents(parents)

        self.children: List[GFF3Record] = []
        if children is not None:
            self.add_children(children)
        return

    def __str__(self) -> str:
        return self.as_str(escape=True)

    def as_str(self, escape: bool = False) -> str:
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

    def add_child(self, child: "GFF3Record") -> None:
        if child not in self.children:
            self.children.append(child)
        if self not in child.parents:
            child.parents.append(self)
        return

    def add_parent(self, parent: "GFF3Record") -> None:
        if parent not in self.parents:
            self.parents.append(parent)

        if self not in parent.children:
            parent.children.append(self)
        return

    def add_children(self, children: Sequence["GFF3Record"]) -> None:
        for child in children:
            self.add_child(child)
        return

    def add_parents(self, parents: Sequence["GFF3Record"]) -> None:
        for parent in parents:
            self.add_parent(parent)
        return

    def traverse_children(
        self,
        sort: bool = False,
        breadth: bool = False,
    ) -> Iterator["GFF3Record"]:
        """ A graph traversal of this `GFFRecord`s children.

        Keyword arguments:
        sort -- Sort the children so that they appear in increasing order.
        breadth -- Perform a breadth first search rather than the default
            depth first search.

        Returns:
        A generator yielding `GFFRecord`s.
        """

        should_reverse = not breadth

        seen: Set[GFF3Record] = set()

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
    ) -> Iterator["GFF3Record"]:
        """ A graph traversal of this `GFFRecord`s parents.

        Keyword arguments:
        sort -- Sort the parents so that they appear in increasing order.
        breadth -- Perform a breadth first search rather than the default
            depth first search.

        Returns:
        A generator yielding `GFFRecord`s.
        """

        should_reverse = not breadth

        seen: Set[GFF3Record] = set()
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
        unescape: bool = False,
    ) -> "GFF3Record":
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

        if sline_len == columns_len - 1:
            logger.warning(
                ("Line has has too few columns columns. "
                 "Probably it is missing the attributes"),
            )
        elif sline_len < columns_len:
            raise ValueError((
                "Line has too few columns. "
                f"Expected: {columns_len}, Encountered: {sline_len}"
            ))
        elif sline_len > columns_len:
            logger.warning(
                "Line has too many columns. Expected: %s, Encountered: %s",
                columns_len,
                sline_len
            )

        fields: Dict[str, str] = dict(zip(cls.columns, sline))
        if sline_len == columns_len - 1:
            fields["attributes"] = ""

        # 0-based indexing exclusive
        start = parse_int(fields["start"], "start") - 1
        end = parse_int(fields["end"], "end")

        if start > end:
            tmp = start
            start = end
            end = tmp
            del tmp

        score = parse_or_none(fields["score"], "score", ".", parse_float)

        strand = Strand.parse(is_one_of(
            fields["strand"],
            ["-", "+", ".", "?"],
            "strand"
        ))

        phase = Phase.parse(is_one_of(
            fields["phase"],
            ["0", "1", "2", "."],
            "phase"
        ))

        attributes = GFF3Attributes.parse(
            fields["attributes"],
            strip_quote=strip_quote,
            unescape=unescape,
        )

        return cls(
            parse_string_not_empty(fields["seqid"], "seqid"),
            parse_string_not_empty(fields["source"], "source"),
            parse_string_not_empty(fields["type"], "type"),
            start,
            end,
            score,
            strand,
            phase,
            attributes
        )

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

    def update_parents(self) -> None:
        self.add_attributes_if_none()
        assert self.attributes is not None

        parent_ids = []
        for parent in self.parents:
            assert parent.attributes is not None
            assert parent.attributes.id is not None
            parent_id = parent.attributes.id
            parent_ids.append(parent_id)

        self.attributes.parent = parent_ids
        return

    def add_attributes_if_none(
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
        is_circular: Optional[bool] = None,
        custom: Optional[Mapping[str, str]] = None,
    ) -> None:
        if self.attributes is None:
            self.attributes = GFF3Attributes(
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
        return

    def copy(self) -> "GFF3Record":
        """ You'll still need to update the ID """
        from copy import copy, deepcopy

        node_copy = copy(self)
        node_copy.attributes = deepcopy(self.attributes)

        if node_copy.attributes is not None:
            node_copy.attributes.parent = []

        node_copy.children = []
        node_copy.parents = []
        return node_copy
