#!/usr/bin/env python3

from typing import Optional
from typing import TextIO
from typing import Iterator

from predectorutils.analyses.base import Analysis
from predectorutils.analyses.parsers import ParseError, LineParseError
from predectorutils.analyses.parsers import (
    parse_string_not_empty,
    parse_int,
    parse_bool,
    MULTISPACE_REGEX
)


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
            parse_string_not_empty(sline[0], "name"),
            parse_int(sline[1], "tm"),
            parse_bool(sline[2], "sp", "Y", "0"),
            parse_string_not_empty(sline[3], "topology")
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

            except LineParseError as e:
                if hasattr(handle, "name"):
                    filename: Optional[str] = handle.name
                else:
                    filename = None

                raise ParseError(
                    filename,
                    i,
                    e.message
                )
        return
