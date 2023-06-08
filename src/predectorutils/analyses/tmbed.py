#!/usr/bin/env python3

from typing import TypeVar
from typing import TextIO
from collections.abc import Iterator, Sequence

from ..gff import (
    GFFRecord,
    Strand
)

from ..parsers import (
    FieldParseError,
    BlockParseError,
    parse_field,
    raise_it,
    parse_str,
    parse_sequence,
)

from .base import Analysis, GFFAble

__all__ = ["TMBed"]


tm_name = raise_it(parse_field(parse_str, "name"))
tm_topology = raise_it(parse_field(
    parse_sequence(["S", "B", "b", "H", "h", "."]),
    "topology"
))


T = TypeVar("T")


def parse_topology(s: str):

    if len(s) == 0:
        return s

    current_type = s[0]
    start = 0

    i = 1
    while i < len(s):
        if s[i] != current_type:
            yield (current_type, start, i)

            current_type = s[i]
            start = i

        i += 1

    yield (current_type, start, i)

    return


class TMBed(Analysis, GFFAble):

    """ .
    """

    columns = ["name", "has_sp", "has_tm", "topology"]
    types = [str, bool, bool, str]
    analysis = "tmbed"
    software = "TMBed"

    def __init__(
        self,
        name: str,
        has_sp: bool,
        has_tm: bool,
        topology: str,
    ) -> None:
        self.name = name
        self.has_sp = has_sp
        self.has_tm = has_tm
        self.topology = topology
        return

    @classmethod
    def from_block(cls, lines: Sequence[str]) -> "TMBed":
        """ Parse a tmbed line as an object. """

        ilines = list(lines)

        try:
            name = tm_name(ilines[0])
        except FieldParseError as e:
            raise e.as_block_error(0)

        name = name.lstrip(">")

        try:
            top = tm_topology(ilines[2].strip())
        except FieldParseError as e:
            raise e.as_block_error(2)

        # Slice is because model doesn't seem to have proper statemodel
        # so i think SPs could _potentially_ happen in middle of protein.
        has_sp = top.lower()[:30].count("s") > 10

        longest_run = {"B": 0, "b": 0, "H": 0, "h": 0, ".": 0, "S": 0}
        prev = top[:1]  # Slice prevents error if empty
        current_run = 1
        for this in top[1:]:
            if this == prev:
                current_run += 1
            elif (this != prev) and (current_run > longest_run[prev]):
                longest_run[prev] = current_run
                prev = this
                current_run = 1
            else:
                prev = this
                current_run = 1

        if current_run > longest_run[prev]:
            longest_run[prev] = current_run

        del longest_run["."]
        del longest_run["S"]

        has_tm = any([v > 10 for v in longest_run.values()])

        return cls(name, has_sp, has_tm, top)

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["TMBed"]:
        block: list[str] = []

        # Avoid case where handle is empty and we raise BlockParseError
        i = 0

        for i, line in enumerate(handle):
            sline = line.strip()
            if sline.startswith("#"):
                continue
            elif sline == "":
                continue
            elif sline.startswith(">") and len(block) > 0:
                try:
                    yield cls.from_block(block)
                except BlockParseError as e:
                    raise (
                        e.as_parse_error(line=i - len(block))
                        .add_filename_from_handle(handle)
                    )
                block = [sline]
            else:
                block.append(sline)

        if len(block) > 0:
            try:
                yield cls.from_block(block)
            except BlockParseError as e:
                raise (
                    e.as_parse_error(line=i - len(block))
                    .add_filename_from_handle(handle)
                )

        return

    def as_gff(
        self,
        software_version: str | None = None,
        database_version: str | None = None,
        keep_all: bool = False,
        id_index: int = 1
    ) -> Iterator[GFFRecord]:
        for (type_, start, end) in parse_topology(self.topology):
            mapp = {
                "B": "transmembrane_polypeptide_region",
                "b": "transmembrane_polypeptide_region",
                "H": "transmembrane_polypeptide_region",
                "h": "transmembrane_polypeptide_region",
                "S": "signal_peptide",
            }
            if type_ == ".":
                continue

            yield GFFRecord(
                seqid=self.name,
                source=self.gen_source(software_version, database_version),
                type=mapp[type_],
                start=start,
                end=end,
                strand=Strand.UNSTRANDED,
                attributes=None
            )
        return
