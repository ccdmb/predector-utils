#!/usr/bin/env python3

from typing import TextIO
from typing import Iterator

from predectorutils.analyses.base import Analysis
from predectorutils.parsers import (
    FieldParseError,
    LineParseError,
    raise_it,
    parse_field,
    parse_str,
    parse_float,
    is_one_of
)

__all__ = ["ApoplastP"]


apo_name = raise_it(parse_field(parse_str, "name"))
apo_prediction = raise_it(parse_field(
    is_one_of(["Apoplastic", "Non-apoplastic"]),
    "prediction"
))
apo_prob = raise_it(parse_field(parse_float, "prob"))


class ApoplastP(Analysis):

    """     """

    columns = ["name", "prediction", "prob"]
    types = [str, str, float]
    analysis = "apoplastp"
    software = "ApoplastP"

    def __init__(self, name: str, prediction: str, prob: float) -> None:
        self.name = name
        self.prediction = prediction
        self.prob = prob
        return

    @classmethod
    def from_line(cls, line: str) -> "ApoplastP":
        """ Parse an ApoplastP line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) != 3:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 3 but got {len(sline)}"
            )

        return cls(
            apo_name(sline[0]),
            apo_prediction(sline[1]),
            apo_prob(sline[2]),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["ApoplastP"]:
        comment = False
        for i, line in enumerate(handle):
            sline = line.strip()
            if comment and sline.startswith("---------"):
                comment = False
                continue
            elif comment:
                continue
            elif (i == 0) and sline.startswith("---------"):
                comment = True
                continue

            if sline.startswith("#"):
                continue
            elif sline == "":
                continue

            try:
                yield cls.from_line(sline)

            except (LineParseError, FieldParseError) as e:
                raise e.as_parse_error(line=i).add_filename_from_handle(handle)
        return
