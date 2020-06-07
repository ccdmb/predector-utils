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


e1_name = raise_it(parse_field(parse_str, "name"))
e1_prediction = raise_it(parse_field(
    is_one_of(["Effector", "Non-effector"]),
    "prediction"
))
e1_prob = raise_it(parse_field(parse_float, "prob"))


class EffectorP1(Analysis):

    """ """

    columns = ["name", "prediction", "prob"]
    types = [str, str, float]
    analysis = "effectorp1"
    software = "EffectorP"

    def __init__(self, name: str, prediction: str, prob: float) -> None:
        self.name = name
        self.prediction = prediction
        self.prob = prob
        return

    @classmethod
    def from_line(cls, line: str) -> "EffectorP1":
        """ Parse an EffectorP1 line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) != 3:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 3 but got {len(sline)}."
            )

        return cls(
            e1_name(sline[0]),
            e1_prediction(sline[1]),
            e1_prob(sline[2]),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["EffectorP1"]:
        comment = False
        for i, line in enumerate(handle):
            sline = line.strip()
            if comment and sline.startswith("---------"):
                comment = False
                continue
            elif comment and sline.startswith("# Identifier"):
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


e2_name = raise_it(parse_field(parse_str, "name"))
e2_prediction = raise_it(parse_field(
    is_one_of(["Effector", "Unlikely effector", "Non-effector"]),
    "prediction"
))
e2_prob = raise_it(parse_field(parse_float, "prob"))


class EffectorP2(Analysis):

    """ """

    columns = ["name", "prediction", "prob"]
    types = [str, str, float]
    analysis = "effectorp2"
    software = "EffectorP"

    def __init__(self, name: str, prediction: str, prob: float) -> None:
        self.name = name
        self.prediction = prediction
        self.prob = prob
        return

    @classmethod
    def from_line(cls, line: str) -> "EffectorP2":
        """ Parse an EffectorP2 line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) != 3:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 3 but got {len(sline)}."
            )

        return cls(
            e2_name(sline[0]),
            e2_prediction(sline[1]),
            e2_prob(sline[2]),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["EffectorP2"]:
        comment = False
        for i, line in enumerate(handle):
            sline = line.strip()
            if comment and sline.startswith("---------"):
                comment = False
                continue
            elif comment and sline.startswith("# Identifier"):
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
