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
    parse_float,
    parse_int,
    split_at_eq,
)

tm_name = raise_it(parse_field(parse_str, "name"))
tm_length = raise_it(parse_field(split_at_eq(parse_int, "len"), "length"))
tm_exp_aa = raise_it(parse_field(split_at_eq(parse_float, "ExpAA"), "exp_aa"))
tm_first_60 = raise_it(parse_field(
    split_at_eq(parse_float, "First60"),
    "first_60"
))
tm_pred_hel = raise_it(parse_field(
    split_at_eq(parse_int, "PredHel"),
    "pred_hel"
))
tm_topology = raise_it(parse_field(
    split_at_eq(parse_str, "Topology"),
    "topology"
))


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
            tm_name(sline[0]),
            tm_length(sline[1]),
            tm_exp_aa(sline[2]),
            tm_first_60(sline[3]),
            tm_pred_hel(sline[4]),
            tm_topology(sline[5]),
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
            except (LineParseError, FieldParseError) as e:
                raise e.as_parse_error(line=i).add_filename_from_handle(handle)

        return
