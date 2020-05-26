#!/usr/bin/env python3

from typing import Optional
from typing import TextIO
from typing import Iterator

from predectorutils.analyses.base import Analysis
from predectorutils.analyses.base import float_or_none, str_or_none
from predectorutils.analyses.parsers import ParseError, LineParseError
from predectorutils.analyses.parsers import (
    parse_string_not_empty,
    parse_float,
    is_one_of
)


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

        prediction = is_one_of(
            sline[1],
            ["noTP", "SP", "mTP"],
            "prediction"
        )

        if prediction == "noTP":
            prediction = "OTHER"

        return cls(
            parse_string_not_empty(sline[0], "name"),
            prediction,
            parse_float(sline[2], "OTHER"),
            parse_float(sline[3], "SP"),
            parse_float(sline[4], "mTP"),
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
            parse_string_not_empty(sline[0], "name"),
            is_one_of(
                sline[1],
                ["OTHER", "SP", "mTP", "cTP", "luTP"],
                "prediction"
            ),
            parse_float(sline[2], "OTHER"),
            parse_float(sline[3], "SP"),
            parse_float(sline[4], "mTP"),
            parse_float(sline[5], "cTP"),
            parse_float(sline[6], "luTP"),
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
