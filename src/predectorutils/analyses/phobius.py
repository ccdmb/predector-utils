#!/usr/bin/env python3

from typing import TextIO
from typing import Iterator

from predectorutils.analyses.base import Analysis
from predectorutils.parsers import (
    FieldParseError,
    LineParseError,
    parse_field,
    raise_it,
    parse_str,
    parse_int,
    parse_bool,
    MULTISPACE_REGEX
)

pb_name = raise_it(parse_field(parse_str, "name"))
pb_tm = raise_it(parse_field(parse_int, "tm"))
pb_sp = raise_it(parse_field(parse_bool("Y", "0"), "sp"))
pb_topology = raise_it(parse_field(parse_str, "topology"))


class Phobius(Analysis):

    """ .
    """

    columns = ["name", "tm", "sp", "topology"]
    types = [str, int, bool, str]
    analysis = "phobius"
    software = "Phobius"

    def __init__(self, name: str, tm: int, sp: bool, topology: str) -> None:
        self.name = name
        self.tm = tm
        self.sp = sp
        self.topology = topology
        return

    @classmethod
    def from_line(cls, line: str) -> "Phobius":
        """ Parse a phobius line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = MULTISPACE_REGEX.split(line.strip())

        if len(sline) != 4:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 4 but got {len(sline)}"
            )

        # Sequence is mis-spelled in the output
        if sline == ["SEQENCE", "ID", "TM", "SP", "PREDICTION"]:
            raise LineParseError("The line appears to be the header line")

        return cls(
            pb_name(sline[0]),
            pb_tm(sline[1]),
            pb_sp(sline[2]),
            pb_topology(sline[3]),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["Phobius"]:
        for i, line in enumerate(handle):
            sline = line.strip()
            if sline.startswith("#"):
                continue
            # Sequence is mis-spelled in the output
            elif i == 0 and sline.startswith("SEQENCE"):
                continue
            elif sline == "":
                continue

            try:
                yield cls.from_line(sline)
            except (LineParseError, FieldParseError) as e:
                raise e.as_parse_error(line=i).add_filename_from_handle(handle)
        return
