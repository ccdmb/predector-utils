#!/usr/bin/env python3

from typing import Optional
from typing import TextIO
from typing import Iterator

from predectorutils.analyses.base import Analysis
from predectorutils.analyses.parsers import ParseError, LineParseError
from predectorutils.analyses.parsers import (
    parse_string_not_empty,
    parse_float,
    parse_int,
    split_at_eq,
)


class TMHMM(Analysis):

    """ .
    """

    columns = ["name", "length", "exp_aa", "first_60", "pred_hel", "topology"]
    types = [str, int, float, float, int, str]
    analysis = "tmhmm"
    software = "TMHMM"

    def __init__(
        self,
        name: str,
        length: int,
        exp_aa: float,
        first_60: float,
        pred_hel: int,
        topology: str,
    ) -> None:
        self.name = name
        self.length = length
        self.exp_aa = exp_aa
        self.first_60 = first_60
        self.pred_hel = pred_hel
        self.topology = topology
        return

    @classmethod
    def from_line(cls, line: str) -> "TMHMM":
        """ Parse a tmhmm line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) != 6:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 6 but got {len(sline)}"
            )

        return cls(
            parse_string_not_empty(sline[0], "name"),
            parse_int(split_at_eq(sline[1], "length", "len"), "length"),
            parse_float(split_at_eq(sline[2], "exp_aa", "ExpAA"), "exp_aa"),
            parse_float(
                split_at_eq(sline[3], "first_60", "First60"),
                "first_60"
            ),
            parse_int(
                split_at_eq(sline[4], "pred_hel", "PredHel"),
                "pred_hel"
            ),
            parse_string_not_empty(
                split_at_eq(sline[5], "topology", "Topology"),
                "topology"
            )
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["TMHMM"]:
        for i, line in enumerate(handle):
            sline = line.strip()
            if sline.startswith("#"):
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
