#!/usr/bin/env python3

import sys

from typing import Optional
from typing import TextIO
from typing import Iterator

from predectorutils.analyses.base import Analysis
from predectorutils.analyses.base import str_or_none
from predectorutils.analyses.parsers import ParseError, LineParseError
from predectorutils.analyses.parsers import (
    parse_string_not_empty,
    parse_float,
    parse_int,
    parse_bool,
    MULTISPACE_REGEX,
    is_one_of,
    is_value
)


class SignalP3NN(Analysis):

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

    columns = [
        "name", "cmax", "cmax_pos", "cmax_decision", "ymax", "ymax_pos",
        "ymax_decision", "smax", "smax_pos", "smax_decision", "smean",
        "smean_decision", "d", "d_decision"
    ]

    types = [str, float, int, bool, float, int, bool, float, int,
             bool, float, bool, float, bool]
    analysis = "signalp3_nn"
    software = "SignalP"

    def __init__(
        self,
        name: str,
        cmax: float,
        cmax_pos: int,
        cmax_decision: bool,
        ymax: float,
        ymax_pos: int,
        ymax_decision: bool,
        smax: float,
        smax_pos: int,
        smax_decision: bool,
        smean: float,
        smean_decision: bool,
        d: float,
        d_decision: bool,
    ):
        self.name = name
        self.cmax = cmax
        self.cmax_pos = cmax_pos
        self.cmax_decision = cmax_decision
        self.ymax = ymax
        self.ymax_pos = ymax_pos
        self.ymax_decision = ymax_decision
        self.smax = smax
        self.smax_pos = smax_pos
        self.smax_decision = smax_decision
        self.smean = smean
        self.smean_decision = smean_decision
        self.d = d
        self.d_decision = d_decision
        return

    @classmethod
    def from_line(cls, line: str) -> "SignalP3NN":
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
            parse_string_not_empty(sline[0], "name"),
            parse_float(sline[1], "cmax"),
            parse_int(sline[2], "cmax_pos"),
            parse_bool(sline[3], "cmax_decision", "Y", "N"),
            parse_float(sline[4], "ymax"),
            parse_int(sline[5], "ymax_pos"),
            parse_bool(sline[6], "ymax_decision", "Y", "N"),
            parse_float(sline[7], "smax"),
            parse_int(sline[8], "smax_pos"),
            parse_bool(sline[9], "smax_decision", "Y", "N"),
            parse_float(sline[10], "smean"),
            parse_bool(sline[11], "smean_decision", "Y", "N"),
            parse_float(sline[12], "d"),
            parse_bool(sline[13], "d_decision", "Y", "N"),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["SignalP3NN"]:
        for i, line in enumerate(handle):
            sline = line.strip()
            if sline.startswith("#"):
                continue
            elif sline == "":
                continue

            if sline == "error running HOW":
                print(
                    f"Encountered an error message on line {i}.",
                    file=sys.stderr
                )
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


class SignalP3HMM(Analysis):

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

    columns = ["name", "is_secreted", "cmax", "cmax_pos", "cmax_decision",
               "sprob", "sprob_decision"]
    types = [str, bool, float, int, bool, float, int]
    analysis = "signalp3_hmm"
    software = "SignalP"

    def __init__(
        self,
        name: str,
        is_secreted: bool,
        cmax: float,
        cmax_pos: int,
        cmax_decision: bool,
        sprob: float,
        sprob_decision: int,
    ) -> None:
        self.name = name
        self.is_secreted = is_secreted
        self.cmax = cmax
        self.cmax_pos = cmax_pos
        self.cmax_decision = cmax_decision
        self.sprob = sprob
        self.sprob_decision = sprob_decision
        return

    @classmethod
    def from_line(cls, line: str) -> "SignalP3HMM":
        """ Parse a short-format HMM line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = MULTISPACE_REGEX.split(line)

        if len(sline) != 7:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 7 but got {len(sline)}"
            )

        # in column !.
        # Q is non-secreted, A is something, possibly long signalpeptide?
        return cls(
            parse_string_not_empty(sline[0], "name"),
            is_value(sline[1], "is_secreted (!)", "S"),
            parse_float(sline[2], "cmax"),
            parse_int(sline[3], "cmax_pos"),
            parse_bool(sline[4], "cmax_decision", "Y", "N"),
            parse_float(sline[5], "sprob"),
            parse_bool(sline[6], "sprob_decision", "Y", "N"),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["SignalP3HMM"]:
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


class SignalP4(Analysis):

    """ The graphical output from SignalP (neural network) comprises
    three different scores, C, S and Y. Two additional scores are reported
    in the SignalP output, namely the S-mean and the D-score, but these
    are only reported as numerical values.

    For each organism class in SignalP; Eukaryote, Gram-negative and
    Gram-positive, two different neural networks are used, one for
    predicting the actual signal peptide and one for predicting the
    position of the signal peptidase I (SPase I) cleavage site.
    The S-score for the signal peptide prediction is reported for every
    single amino acid position in the submitted sequence, with high
    scores indicating that the corresponding amino acid is part of a
    signal peptide, and low scores indicating that the amino acid is
    part of a mature protein.

    The C-score is the ``cleavage site'' score. For each position in the
    submitted sequence, a C-score is reported, which should only be
    significantly high at the cleavage site. Confusion is often seen
    with the position numbering of the cleavage site. When a cleavage
    site position is referred to by a single number, the number indicates
    the first residue in the mature protein, meaning that a reported
    cleavage site between amino acid 26-27 corresponds to that the mature
    protein starts at (and include) position 27.

    Y-max is a derivative of the C-score combined with the S-score
    resulting in a better cleavage site prediction than the raw C-score
    alone. This is due to the fact that multiple high-peaking C-scores can
    be found in one sequence, where only one is the true cleavage site.
    The cleavage site is assigned from the Y-score where the slope of the
    S-score is steep and a significant C-score is found.

    The S-mean is the average of the S-score, ranging from the N-terminal
    amino acid to the amino acid assigned with the highest Y-max score,
    thus the S-mean score is calculated for the length of the predicted
    signal peptide. The S-mean score was in SignalP version 2.0 used as
    the criteria for discrimination of secretory and non-secretory proteins.

    The D-score is introduced in SignalP version 3.0. In version 4.0 this
    score is implemented as a weighted average of the S-mean and the
    Y-max scores. The score shows superior discrimination performance of
    secretory and non-secretory proteins to that of the S-mean score which
    was used in SignalP version 1 and 2.

    For non-secretory proteins all the scores represented in the SignalP
    output should ideally be very low.
   """

    types = [str, float, int, float, int, float, int,
             float, float, bool, float, str, str]
    columns = ["name", "cmax", "cmax_pos", "ymax", "ymax_pos",
               "smax", "smax_pos", "smean", "d", "decision", "dmax_cut",
               "networks_used"]
    analysis = "signalp4"
    software = "SignalP"

    def __init__(
        self,
        name: str,
        cmax: float,
        cmax_pos: int,
        ymax: float,
        ymax_pos: int,
        smax: float,
        smax_pos: int,
        smean: float,
        d: float,
        decision: bool,
        dmax_cut: float,
        networks_used: str,
    ) -> None:
        self.name = name
        self.cmax = cmax
        self.cmax_pos = cmax_pos
        self.ymax = ymax
        self.ymax_pos = ymax_pos
        self.smax = smax
        self.smax_pos = smax_pos
        self.smean = smean
        self.d = d
        self.decision = decision
        self.dmax_cut = dmax_cut
        self.networks_used = networks_used
        return

    @classmethod
    def from_line(cls, line: str) -> "SignalP4":
        """ Parse a short-format signalp4 line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = MULTISPACE_REGEX.split(line)

        if len(sline) != 12:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 12 but got {len(sline)}"
            )

        return cls(
            parse_string_not_empty(sline[0], "name"),
            parse_float(sline[1], "cmax"),
            parse_int(sline[2], "cmax_pos"),
            parse_float(sline[3], "ymax"),
            parse_int(sline[4], "ymax_pos"),
            parse_float(sline[5], "smax"),
            parse_int(sline[6], "smax_pos"),
            parse_float(sline[7], "smean"),
            parse_float(sline[8], "d"),
            parse_bool(sline[9], "decision", "Y", "N"),
            parse_float(sline[10], "dmax_cut"),
            is_one_of(
                sline[11],
                ["SignalP-noTM", "SignalP-TM"],
                "networks_used"
            ),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["SignalP4"]:
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


class SignalP5(Analysis):

    """ One annotation is attributed to each protein, the one that has
    the highest probability. The protein can have a Sec signal peptide
    (Sec/SPI), a Lipoprotein signal peptide (Sec/SPII), a Tat signal
    peptide (Tat/SPI) or No signal peptide at all (Other).

    If a signal peptide is predicted, the cleavage site position is
    reported as well. On the plot, three marginal probabilities are
    reported, i.e. SP(Sec/SPI) / LIPO(Sec/SPII) / TAT(Tat/SPI)
    (depending on what type of signal peptide is predicted), CS (the
    cleavage site) and OTHER (the probability that the sequence does
    not have any kind of signal peptide).
    """

    name: str
    prediction: str
    prob_signal: float
    prob_other: float
    cs_pos: Optional[str]

    columns = ["name", "prediction", "prob_signal", "prob_other", "cs_pos"]
    types = [str, str, float, float, str_or_none]
    analysis = "signalp5"
    software = "SignalP"

    def __init__(
        self,
        name: str,
        prediction: str,
        prob_signal: float,
        prob_other: float,
        cs_pos: Optional[str],
    ) -> None:
        self.name = name
        self.prediction = prediction
        self.prob_signal = prob_signal
        self.prob_other = prob_other
        self.cs_pos = cs_pos
        return

    @classmethod
    def from_line(cls, line: str) -> "SignalP5":
        """ Parse a short-format signalp5 line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) == 5:
            cs_pos: Optional[str] = str(sline[4])
        elif len(sline) == 4:
            cs_pos = None
        else:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 4 or 5 but got {len(sline)}"
            )

        return cls(
            parse_string_not_empty(sline[0], "name"),
            is_one_of(
                sline[1],
                ["SP(Sec/SPI)", "LIPO(Sec/SPII)", "TAT(Tat/SPI)", "OTHER"],
                "prediction"
            ),
            parse_float(sline[2], "prob_signal"),
            parse_float(sline[3], "prob_other"),
            cs_pos,
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["SignalP5"]:
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
