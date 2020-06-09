#!/usr/bin/env python3

from typing import Optional
from typing import TextIO
from typing import Iterator

from predectorutils.analyses.base import Analysis
from predectorutils.analyses.base import float_or_none, str_or_none
from predectorutils.parsers import (
    FieldParseError,
    LineParseError,
    parse_field,
    raise_it,
    parse_str,
    parse_float,
    is_one_of
)


tp_name = raise_it(parse_field(parse_str, "name"))
tp_prediction = raise_it(parse_field(
    is_one_of(["noTP", "SP", "mTP"]),
    "prediction"
))
tp_other = raise_it(parse_field(parse_float, "OTHER"))
tp_sp = raise_it(parse_field(parse_float, "SP"))
tp_mtp = raise_it(parse_field(parse_float, "mTP"))


pl_prediction = raise_it(parse_field(
    is_one_of(["OTHER", "SP", "mTP", "cTP", "luTP"]),
    "prediction"
))
pl_ctp = raise_it(parse_field(parse_float, "cTP"))
pl_lutp = raise_it(parse_field(parse_float, "luTP"))


class TargetPNonPlant(Analysis):

    """ Doesn't have output format documentation yet
    """

    columns = ["name", "prediction", "other", "sp", "mtp", "cs_pos"]
    types = [str, str, float, float, float, str_or_none]
    analysis = "targetp_nonplant"
    software = "TargetP"

    def __init__(
        self,
        name: str,
        prediction: str,
        other: float,
        sp: float,
        mtp: float,
        cs_pos: Optional[str],
    ) -> None:
        self.name = name
        self.prediction = prediction
        self.other = other
        self.sp = sp
        self.mtp = mtp
        self.cs_pos = cs_pos
        return

    @classmethod
    def from_line(cls, line: str) -> "TargetPNonPlant":
        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) == 6:
            cs_pos: Optional[str] = str(sline[5])
        elif len(sline) == 5:
            cs_pos = None
        else:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 5 or 6 but got {len(sline)}"
            )

        prediction = tp_prediction(sline[1])
        if prediction == "noTP":
            prediction = "OTHER"

        return cls(
            tp_name(sline[0]),
            prediction,
            tp_other(sline[2]),
            tp_sp(sline[3]),
            tp_mtp(sline[4]),
            cs_pos=cs_pos,
        )

    @classmethod
    def from_file(
        cls,
        handle: TextIO,
    ) -> Iterator["TargetPNonPlant"]:
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


class TargetPPlant(Analysis):

    """ Doesn't have output format documentation yet
    """

    columns = ["name", "prediction", "other", "sp",
               "mtp", "ctp", "lutp", "cs_pos"]
    types = [str, str, float, float, float,
             float_or_none, float_or_none, str_or_none]
    analysis = "targetp_plant"
    software = "TargetP"

    def __init__(
        self,
        name: str,
        prediction: str,
        other: float,
        sp: float,
        mtp: float,
        ctp: Optional[float],
        lutp: Optional[float],
        cs_pos: Optional[str],
    ) -> None:
        self.name = name
        self.prediction = prediction
        self.other = other
        self.sp = sp
        self.mtp = mtp
        self.ctp = ctp
        self.lutp = lutp
        self.cs_pos = cs_pos
        return

    @classmethod
    def from_line(cls, line: str) -> "TargetPPlant":
        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) == 8:
            cs_pos: Optional[str] = str(sline[7])
        elif len(sline) == 7:
            cs_pos = None
        else:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 7 or 8 but got {len(sline)}"
            )

        return cls(
            tp_name(sline[0]),
            pl_prediction(sline[1]),
            tp_other(sline[2]),
            tp_sp(sline[3]),
            tp_mtp(sline[4]),
            pl_ctp(sline[5]),
            pl_lutp(sline[6]),
            cs_pos,
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["TargetPPlant"]:
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
