#!/usr/bin/env python3

from typing import Optional
from typing import TextIO
from typing import Iterator

from predector.analyses import Analysis
from predector.analyses.parsers import ParseError, LineParseError
from predector.analyses.parsers import (
    parse_string_not_empty,
    parse_float,
    parse_int,
    is_one_of
)


class DeepSig(Analysis):

    """     """

    def __init__(
        self,
        name: str,
        prediction: str,
        prob: float,
        cs_pos: Optional[int]
    ) -> None:
        self.name = name
        self.prediction = prediction
        self.prob = prob
        self.cs_pos = cs_pos
        return

    @classmethod
    def from_line(cls, line: str) -> "DeepSig":
        """ Parse a deepsig line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) != 4:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 4 but got {len(sline)}"
            )

        if sline[3] == "-":
            cs_pos: Optional[int] = None
        else:
            cs_pos = parse_int(sline[3], "cs_pos")

        return cls(
            parse_string_not_empty(sline[0], "name"),
            is_one_of(
                sline[1],
                ["SignalPeptide", "Transmembrane", "Other"],
                "prediction"
            ),
            parse_float(sline[2], "prob"),
            cs_pos,
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["DeepSig"]:
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
