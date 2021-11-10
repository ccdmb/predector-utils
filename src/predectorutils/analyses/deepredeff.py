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
    is_one_of
)

dre_name = raise_it(parse_field(parse_str, "name"))
dre_s_score = raise_it(parse_field(parse_float, "s_score"))
dre_prediction = raise_it(parse_field(
    is_one_of(["effector", "non-effector"]),
    "prediction"
))


class Deepredeff(Analysis):

    """ """
    columns = [
        "name",
        "s_score",
        "prediction",
    ]

    types = [
        str,
        float,
        str,
    ]

    software = "deepredeff"
    analysis = "deepredeff"
    name_column = "name"

    def __init__(
        self,
        name: str,
        s_score: float,
        prediction: str,
    ):
        self.name = name
        self.s_score = s_score
        self.prediction = prediction
        return

    @classmethod
    def from_line(cls, line: str) -> "Deepredeff":
        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t", maxsplit=3)
        if len(sline) != 3:
            # Technically because of the max_split this should be impossible.
            # the description line is allowed to have spaces.
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 3 but got {len(sline)}"
            )

        return cls(
            dre_name(sline[0]),
            dre_s_score(sline[1]),
            dre_prediction(sline[2]),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["Deepredeff"]:
        for i, line in enumerate(handle):
            sline = line.strip()

            if sline.startswith("#"):
                continue
            elif sline == "":
                continue
            elif sline == "\t".join(cls.columns):
                continue

            try:
                yield cls.from_line(sline)
            except (LineParseError, FieldParseError) as e:
                raise e.as_parse_error(line=i).add_filename_from_handle(handle)
        return


class DeepredeffFungi(Deepredeff):
    analysis = "deepredeff_fungi"


class DeepredeffOomycete(Deepredeff):
    analysis = "deepredeff_oomycete"


class DeepredeffBacteria(Deepredeff):
    analysis = "deepredeff_bacteria"
