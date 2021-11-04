#!/usr/bin/env python3

from typing import TextIO
from typing import Iterator

from predectorutils.gff import GFFRecord, GFFAttributes, Strand
from predectorutils.analyses.base import Analysis, GFFAble
from predectorutils.parsers import (
    FieldParseError,
    LineParseError,
    parse_field,
    raise_it,
    parse_str,
    parse_int
)


re_name = raise_it(parse_field(parse_str, "name"))
re_kind = raise_it(parse_field(parse_str, "kind"))
re_pattern = raise_it(parse_field(parse_str, "pattern"))
re_match = raise_it(parse_field(parse_str, "match"))
re_start = raise_it(parse_field(parse_int, "start"))
re_end = raise_it(parse_field(parse_int, "end"))


class RegexAnalysis(Analysis, GFFAble):

    columns = [
        "name",
        "kind",
        "pattern",
        "match",
        "start",
        "end"
    ]

    types = [
        str,
        str,
        str,
        str,
        int,
        int
    ]

    analysis = "regex"
    software = "Predector"

    def __init__(
        self,
        name: str,
        kind: str,
        pattern: str,
        match: str,
        start: int,
        end: int,
    ) -> None:
        self.name = name
        self.kind = kind
        self.pattern = pattern
        self.match = match
        self.start = start
        self.end = end
        return None

    @classmethod
    def from_line(cls, line: str) -> "RegexAnalysis":
        """ Parse a table line as an object """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = [c.strip() for c in line.strip().split("\t")]

        if len(sline) != 6:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 6 but got {len(sline)}"
            )

        return cls(
            re_name(sline[0]),
            re_kind(sline[1]),
            re_pattern(sline[2]),
            re_match(sline[3]),
            re_start(sline[4]),
            re_end(sline[5])
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["RegexAnalysis"]:
        for i, line in enumerate(handle):
            sline = line.strip()

            if sline.startswith("#"):
                continue
            elif sline == "":
                continue
            elif sline == "\t".join(cls.columns):
                continue

            try:
                yield cls.from_line(sline)
            except (LineParseError, FieldParseError) as e:
                raise e.as_parse_error(line=i).add_filename_from_handle(handle)
        return

    def as_gff(
        self,
        keep_all: bool = True,
        id_index: int = 1,
    ) -> Iterator[GFFRecord]:
        attr = GFFAttributes(custom={
            "kind": self.kind,
            "pattern": self.pattern,
            "match": self.match,
        })

        if self.kind == "Kex2":
            type_ = "propeptide_cleavage_site"
        else:
            type_ = "polypeptide_motif"

        yield GFFRecord(
            seqid=self.name,
            source=self.analysis,
            type=type_,
            start=self.start,
            end=self.end,
            score=None,
            strand=Strand.PLUS,
            attributes=attr
        )
        return