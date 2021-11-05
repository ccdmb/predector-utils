#!/usr/bin/env python3

import re

from typing import Optional
from typing import TextIO
from typing import Iterator

from predectorutils.analyses.base import Analysis
from predectorutils.parsers import (
    FieldParseError,
    LineParseError,
    raise_it,
    parse_field,
    parse_regex,
    parse_str,
    parse_float,
    is_one_of
)
from predectorutils.analyses.base import float_or_none


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


e3_name = raise_it(parse_field(parse_str, "name"))
e3_prediction = raise_it(parse_field(
    is_one_of([
        "Cytoplasmic effector", "Apoplastic effector",
        "Cytoplasmic/apoplastic effector", "Apoplastic/cytoplasmic effector",
        "Non-effector"
    ]),
    "prediction"
))

E3_REGEX = re.compile(r"^Y \((?P<prob>\d?\.?\d+)\)$")


def e3_parse_field(field: str, field_name: str) -> Optional[float]:
    field = field.strip()
    if field == "-":
        return None

    res = parse_field(parse_regex(E3_REGEX), field_name)(field)
    if isinstance(res, FieldParseError):
        raise res

    return float(res["prob"])


class EffectorP3(Analysis):

    """ """

    columns = [
        "name",
        "prediction",
        "cytoplasmic_prob",
        "apoplastic_prob",
        "noneffector_prob",
    ]
    types = [str, str, float_or_none, float_or_none, float_or_none]
    analysis = "effectorp3"
    software = "EffectorP"

    def __init__(
        self,
        name: str,
        prediction: str,
        cytoplasmic_prob: Optional[float],
        apoplastic_prob: Optional[float],
        noneffector_prob: Optional[float],
    ) -> None:
        self.name = name
        self.prediction = prediction
        self.cytoplasmic_prob = cytoplasmic_prob
        self.apoplastic_prob = apoplastic_prob
        self.noneffector_prob = noneffector_prob
        return

    @classmethod
    def from_line(cls, line: str) -> "EffectorP3":
        """ Parse an EffectorP3 line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) != 5:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 5 but got {len(sline)}."
            )

        return cls(
            e3_name(sline[0]),
            e3_prediction(sline[4]),
            e3_parse_field(sline[1], "cytoplasmic_prob"),
            e3_parse_field(sline[2], "apoplastic_prob"),
            e3_parse_field(sline[3], "noneffector_prob"),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["EffectorP3"]:
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


class EffectorP3Fungal(EffectorP3):
    analysis = "effectorp3_fungal"
