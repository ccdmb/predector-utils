#!/usr/bin/env python3

from typing import Dict
from typing import Union, Optional
from typing import TextIO
from typing import Iterator

from predector.analyses import Analysis
from predector.analyses.parsers import ParseError, LineParseError
from predector.analyses.parsers import (
    parse_string_not_empty,
    parse_float,
    is_one_of
)


class ApoplastP(Analysis):

    """     """

    columns = ["name", "prediction", "prob"]
    types = [str, str, float]

    def __init__(self, name: str, prediction: str, prob: float) -> None:
        self.name = name
        self.prediction = prediction
        self.prob = prob

    def as_dict(self) -> Dict[str, Union[str, int, float, bool, None]]:
        return {k: getattr(self, k) for k in self.columns}

    @classmethod
    def from_dict(
        cls,
        d: Dict[str, Union[str, int, float, bool, None]]
    ) -> "ApoplastP":
        # assertions appease the typechecker
        name = d["name"]
        assert isinstance(name, str)

        prediction = d["prediction"]
        assert isinstance(prediction, str)

        prob = d["prob"]
        assert isinstance(prob, float)
        return cls(name, prediction, prob)

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
            parse_string_not_empty(sline[0], "name"),
            is_one_of(
                sline[1],
                ["Apoplastic", "Non-apoplastic"],
                "prediction"
            ),
            parse_float(sline[2], "prob"),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["ApoplastP"]:
        comment = False
        for i, line in enumerate(handle):
            sline = line.strip()
            if comment and line.startswith("---------"):
                comment = False
                continue
            elif comment:
                continue
            elif (i == 0) and line.startswith("---------"):
                comment = True
                continue

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
