#!/usr/bin/env python3

import re
from typing import Optional
from typing import TextIO
from typing import Iterator
from typing import Tuple

from predectorutils.analyses.base import Analysis
from predectorutils.analyses.base import (
    int_or_none,
    float_or_none,
    str_or_none
)
from predectorutils.analyses.parsers import ParseError, LineParseError
from predectorutils.analyses.parsers import (
    parse_string_not_empty,
)

# This matches strings of the form "Y (0.962| 25-49)"
TP_REGEX = re.compile(
    r"^Y \((?P<prob>\d?\.?\d+)\s*\|\s*(?P<start>\d+)-(?P<end>\d+)\)$")
NUC_REGEX = re.compile(r"^Y \((?P<sigs>.+)\)$")
HEADER_REGEX = re.compile(
    r"^Identifier\s+Chloroplast\s+Mitochondria\s+Nucleus$")


def parse_tp_field(
    field: str,
    field_name: str,
) -> Tuple[bool, Optional[float], Optional[int], Optional[int]]:
    field = field.strip()

    if field == "-":
        return (False, None, None, None)
    else:
        res = TP_REGEX.match(field)
        if res is None:
            raise LineParseError(
                f"Invalid value: '{field}' in the column: "
                f"'{field_name}'. "
                "Must have the form '-' or 'Y (0.50 | 25-39)'"
            )
        else:
            dres = res.groupdict()
            return (True, float(dres["prob"]),
                    int(dres["start"]), int(dres["end"]))


def parse_nuc_field(
    field: str,
) -> Tuple[bool, Optional[str]]:
    field = field.strip()

    if field == "-":
        return (False, None)
    else:
        res = NUC_REGEX.match(field)
        if res is None:
            raise LineParseError(f"Invalid value: '{field}' in the column: "
                                 f"nucleus'. "
                                 "Must have the form '-' or 'Y (LGEV,PKPS)'")
        else:
            dres = res.groupdict()
            return (True, dres["sigs"])


class LOCALIZER(Analysis):

    """     """

    columns = ["name", "chloroplast_decision", "chloroplast_prob",
               "chloroplast_start", "chloroplast_end",
               "mitochondria_decision", "mitochondria_prob",
               "mitochondria_start", "mitochondria_end",
               "nucleus_decision", "nucleus_signals"]
    types = [str,
             bool, float_or_none, int_or_none, int_or_none,
             bool, float_or_none, int_or_none, int_or_none,
             bool, str_or_none]
    analysis = "localizer"
    software = "LOCALIZER"

    def __init__(
        self,
        name: str,
        chloroplast_decision: bool,
        chloroplast_prob: Optional[float],
        chloroplast_start: Optional[int],
        chloroplast_end: Optional[int],
        mitochondria_decision: bool,
        mitochondria_prob: Optional[float],
        mitochondria_start: Optional[int],
        mitochondria_end: Optional[int],
        nucleus_decision: bool,
        nucleus_signals: Optional[str],
    ) -> None:
        self.name = name
        self.chloroplast_decision = chloroplast_decision
        self.chloroplast_prob = chloroplast_prob
        self.chloroplast_start = chloroplast_start
        self.chloroplast_end = chloroplast_end
        self.mitochondria_decision = mitochondria_decision
        self.mitochondria_prob = mitochondria_prob
        self.mitochondria_start = mitochondria_start
        self.mitochondria_end = mitochondria_end
        self.nucleus_decision = nucleus_decision
        self.nucleus_signals = nucleus_signals
        return

    @classmethod
    def from_line(cls, line: str) -> "LOCALIZER":
        """ Parse an ApoplastP line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = [c.strip() for c in line.strip().split("\t")]

        if len(sline) != 4:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 4 but got {len(sline)}"
            )

        (cp, cp_prob, cp_start, cp_end) = parse_tp_field(
            sline[1],
            "chloroplast"
        )

        (mt, mt_prob, mt_start, mt_end) = parse_tp_field(
            sline[2],
            "mitochondria"
        )

        (nuc, nuc_sigs) = parse_nuc_field(sline[3])

        return cls(
            parse_string_not_empty(sline[0], "name"),
            cp,
            cp_prob,
            cp_start,
            cp_end,
            mt,
            mt_prob,
            mt_start,
            mt_end,
            nuc,
            nuc_sigs
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["LOCALIZER"]:
        for i, line in enumerate(handle):
            sline = line.strip()
            if sline.startswith("#"):
                continue
            elif sline == "":
                continue
            elif HEADER_REGEX.match(sline) is not None:
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
