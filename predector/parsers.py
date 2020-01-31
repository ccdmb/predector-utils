#!/usr/bin/env python3

import re
from typing import NamedTuple
from typing import Optional, Union
from typing import Sequence, Iterator
from typing import Dict
from typing import TextIO


MULTISPACE_REGEX = re.compile(r"\s+")


class ParseError(Exception):
    """ Some aspect of parsing failed. """

    def __init__(
        self,
        filename: Optional[str],
        line: Optional[int],
        message: str
    ):
        self.filename = filename
        self.line = line
        self.message = message
        return


class LineParseError(Exception):

    def __init__(self, message: str):
        self.message = message
        return


class SignalP3NN(NamedTuple):

    """ For each organism class in SignalP; Eukaryote, Gram-negative and
    Gram-positive, two different neural networks are used, one for
    predicting the actual signal peptide and one for predicting the
    position of the signal peptidase I (SPase I) cleavage site.
    The S-score for the signal peptide prediction is reported for
    every single amino acid position in the submitted sequence,
    with high scores indicating that the corresponding amino acid is part
    of a signal peptide, and low scores indicating that the amino acid is
    part of a mature protein.

    The C-score is the 'cleavage site' score. For each position in the
    submitted sequence, a C-score is reported, which should only be
    significantly high at the cleavage site. Confusion is often seen
    with the position numbering of the cleavage site. When a cleavage
    site position is referred to by a single number, the number indicates
    the first residue in the mature protein, meaning that a reported
    cleavage site between amino acid 26-27 corresponds to that the mature
    protein starts at (and include) position 27.

    Y-max is a derivative of the C-score combined with the S-score
    resulting in a better cleavage site prediction than the raw C-score alone.
    This is due to the fact that multiple high-peaking C-scores can be found
    in one sequence, where only one is the true cleavage site.
    The cleavage site is assigned from the Y-score where the slope of the
    S-score is steep and a significant C-score is found.

    The S-mean is the average of the S-score, ranging from the N-terminal
    amino acid to the amino acid assigned with the highest Y-max score, thus
    the S-mean score is calculated for the length of the predicted signal
    peptide. The S-mean score was in SignalP version 2.0 used as the criteria
    for discrimination of secretory and non-secretory proteins.

    The D-score is introduced in SignalP version 3.0 and is a simple average
    of the S-mean and Y-max score. The score shows superior discrimination
    performance of secretory and non-secretory proteins to that of the S-mean
    score which was used in SignalP version 1 and 2.

    For non-secretory proteins all the scores represented in the SignalP3-NN
    output should ideally be very low.

    The hidden Markov model calculates the probability of whether the
    submitted sequence contains a signal peptide or not. The eukaryotic
    HMM model also reports the probability of a signal anchor, previously
    named uncleaved signal peptides. Furthermore, the cleavage site is
    assigned by a probability score together with scores for the n-region,
    h-region, and c-region of the signal peptide, if such one is found.
    """

    name: str
    cmax: float
    cmax_pos: int
    cmax_decision: bool
    ymax: float
    ymax_pos: int
    ymax_decision: bool
    smax: float
    smax_pos: int
    smax_decision: bool
    smean: float
    smean_decision: bool
    d: float
    d_decision: bool

    def as_dict(self) -> Dict[str, Union[str, int, float, bool]]:
        return {k: getattr(self, k) for k in self._fields}

    @classmethod
    def from_short_line(cls, line: str) -> "SignalP3NN":
        """ Parse a short-format NN line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = MULTISPACE_REGEX.split(line)

        if len(sline) != 14:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 14 but got {len(sline)}"
            )

        return cls(
            parse_string_not_empty(line[0], "name"),
            parse_float(line[1], "cmax"),
            parse_int(line[2], "cmax_pos"),
            parse_bool(line[3], "cmax_decision", "Y", "N"),
            parse_float(line[4], "ymax"),
            parse_int(line[5], "ymax_pos"),
            parse_bool(line[6], "ymax_decision", "Y", "N"),
            parse_float(line[7], "smax"),
            parse_int(line[8], "smax_pos"),
            parse_bool(line[9], "smax_decision", "Y", "N"),
            parse_float(line[10], "smean"),
            parse_bool(line[11], "smean_decision", "Y", "N"),
            parse_float(line[12], "d"),
            parse_bool(line[13], "d_decision", "Y", "N"),
        )

    @classmethod
    def from_short_file(cls, handle: TextIO) -> Iterator["SignalP3NN"]:
        for i, line in enumerate(handle):
            sline = line.strip()
            if sline.startswith("#"):
                continue

            try:
                yield cls.from_short_line(sline)

            except LineParseError as e:
                if hasattr(handle.name):
                    filename: Optional[str] = handle.name
                else:
                    filename = None

                raise ParseError(
                    filename,
                    i,
                    e.message
                )
        return


def parse_int(field: str, field_name: str) -> int:
    """ Parse a string as an integer, raising a custom error if fails. """

    try:
        return int(field)
    except ValueError:
        raise LineParseError(
            f"Could not parse value in {field_name} column as an integer. "
            f"The offending value was: '{field}'."
        )


def parse_float(field: str, field_name: str) -> float:
    """ Parse a string as a float, raising a custom error if fails. """

    try:
        return float(field)
    except ValueError:
        raise LineParseError(
            f"Could not parse value in {field_name} column as a float. "
            f"The offending value was: '{field}'."
        )


def parse_bool(
    field: str,
    field_name: str,
    true_value: str,
    false_value: str
) -> bool:
    """ """

    if field == true_value:
        return True
    elif field == false_value:
        return False
    else:
        raise LineParseError(
            f"Invalid value: '{field}' in the column: '{field_name}'. "
            f"Must be either {true_value} or {false_value}."
        )


def parse_string_not_empty(field: str, field_name: str) -> str:
    """ """

    if field.strip() == "":
        raise LineParseError(f"The value in column: '{field_name}' was empty.")
    else:
        return field


def is_one_of(field: str, options: Sequence[str], field_name: str) -> str:
    """ """

    if field in options:
        return field
    else:
        raise LineParseError(
            f"Invalid value: '{field}' in the column: '{field_name}'. "
            f"Must be one of {options}."
        )
