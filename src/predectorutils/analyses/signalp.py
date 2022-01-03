#!/usr/bin/env python3

import re
import sys

from typing import Optional
from typing import TextIO
from typing import Iterator

from predectorutils.gff import (
    GFFRecord,
    GFFAttributes,
    Strand,
)
from predectorutils.analyses.base import Analysis, GFFAble
from predectorutils.analyses.base import str_or_none
from predectorutils.parsers import (
    FieldParseError,
    LineParseError,
    raise_it,
    parse_field,
    parse_str,
    parse_float,
    parse_int,
    parse_bool,
    parse_regex,
    MULTISPACE_REGEX,
    is_one_of,
    is_value
)


__all__ = ["SignalP3NN", "SignalP3HMM", "SignalP4", "SignalP5"]


s3nn_name = raise_it(parse_field(parse_str, "name"))
s3nn_cmax = raise_it(parse_field(parse_float, "cmax"))
s3nn_cmax_pos = raise_it(parse_field(parse_int, "cmax_pos"))
s3nn_cmax_decision = raise_it(parse_field(
    parse_bool("Y", "N"),
    "cmax_decision"
))
s3nn_ymax = raise_it(parse_field(parse_float, "ymax"))
s3nn_ymax_pos = raise_it(parse_field(parse_int, "ymax_pos"))
s3nn_ymax_decision = raise_it(parse_field(
    parse_bool("Y", "N"),
    "ymax_decision"
))
s3nn_smax = raise_it(parse_field(parse_float, "smax"))
s3nn_smax_pos = raise_it(parse_field(parse_int, "smax_pos"))
s3nn_smax_decision = raise_it(parse_field(
    parse_bool("Y", "N"),
    "smax_decision"
))
s3nn_smean = raise_it(parse_field(parse_float, "smean"))
s3nn_smean_decision = raise_it(parse_field(
    parse_bool("Y", "N"),
    "smean_decision"
))
s3nn_d = raise_it(parse_field(parse_float, "d"))
s3nn_d_decision = raise_it(parse_field(parse_bool("Y", "N"), "d_decision"))


class SignalP3NN(Analysis, GFFAble):

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
            s3nn_name(sline[0]),
            s3nn_cmax(sline[1]),
            s3nn_cmax_pos(sline[2]),
            s3nn_cmax_decision(sline[3]),
            s3nn_ymax(sline[4]),
            s3nn_ymax_pos(sline[5]),
            s3nn_ymax_decision(sline[6]),
            s3nn_smax(sline[7]),
            s3nn_smax_pos(sline[8]),
            s3nn_smax_decision(sline[9]),
            s3nn_smean(sline[10]),
            s3nn_smean_decision(sline[11]),
            s3nn_d(sline[12]),
            s3nn_d_decision(sline[13]),
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

            except (LineParseError, FieldParseError) as e:
                raise e.as_parse_error(line=i).add_filename_from_handle(handle)
        return

    def as_gff(
        self,
        software_version: Optional[str] = None,
        database_version: Optional[str] = None,
        keep_all: bool = False,
        id_index: int = 1,
    ) -> Iterator[GFFRecord]:
        if not self.d_decision:
            return

        # d_decision = prediction of issecreted.
        # ymax = first aa of mature peptide
        attr = GFFAttributes(custom={
            "cmax": str(self.cmax),
            "cmax_pos": str(self.cmax_pos),
            "cmax_decision": "true" if self.cmax_decision else "false",
            "ymax": str(self.ymax),
            "ymax_pos": str(self.ymax_pos),
            "ymax_decision": "true" if self.ymax_decision else "false",
            "smax": str(self.smax),
            "smax_pos": str(self.smax_pos),
            "smax_decision": "true" if self.smax_decision else "false",
            "smean": str(self.smean),
            "smean_decision": "true" if self.smean_decision else "false",
            "d": str(self.d),
            "d_decision": "true" if self.d_decision else "false",
        })

        yield GFFRecord(
            seqid=self.name,
            source=self.gen_source(software_version, database_version),
            type="signal_peptide",
            start=0,
            end=self.ymax_pos - 1,
            score=self.d,
            strand=Strand.PLUS,
            attributes=attr
        )
        return


s3hmm_name = raise_it(parse_field(parse_str, "name"))
s3hmm_is_secreted = raise_it(parse_field(is_value("S"), "is_secreted (!)"))
s3hmm_cmax = raise_it(parse_field(parse_float, "cmax"))
s3hmm_cmax_pos = raise_it(parse_field(parse_int, "cmax_pos"))
s3hmm_cmax_decision = raise_it(parse_field(
    parse_bool("Y", "N"),
    "cmax_decision"
))
s3hmm_sprob = raise_it(parse_field(parse_float, "sprob"))
s3hmm_sprob_decision = raise_it(parse_field(
    parse_bool("Y", "N"),
    "sprob_decision"
))


class SignalP3HMM(Analysis, GFFAble):

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
    types = [str, bool, float, int, bool, float, bool]
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
        sprob_decision: bool,
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
            s3hmm_name(sline[0]),
            s3hmm_is_secreted(sline[1]),
            s3hmm_cmax(sline[2]),
            s3hmm_cmax_pos(sline[3]),
            s3hmm_cmax_decision(sline[4]),
            s3hmm_sprob(sline[5]),
            s3hmm_sprob_decision(sline[6]),
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

            except (FieldParseError, LineParseError) as e:
                raise e.as_parse_error(line=i).add_filename_from_handle(handle)
        return

    def as_gff(
        self,
        software_version: Optional[str] = None,
        database_version: Optional[str] = None,
        keep_all: bool = False,
        id_index: int = 1,
    ) -> Iterator[GFFRecord]:
        if not self.is_secreted:
            return

        # d_decision = prediction of issecreted.
        # ymax = first aa of mature peptide
        attr = GFFAttributes(custom={
            "is_secreted": "true" if self.is_secreted else "false",
            "cmax": str(self.cmax),
            "cmax_pos": str(self.cmax_pos),
            "cmax_decision": "true" if self.cmax_decision else "false",
            "sprob": str(self.sprob),
            "sprob_decision": "true" if self.sprob_decision else "false",
        })

        yield GFFRecord(
            seqid=self.name,
            source=self.gen_source(software_version, database_version),
            type="signal_peptide",
            start=0,
            end=self.cmax_pos - 1,
            score=self.sprob,
            strand=Strand.PLUS,
            attributes=attr
        )
        return


s4_name = raise_it(parse_field(parse_str, "name"))
s4_cmax = raise_it(parse_field(parse_float, "cmax"))
s4_cmax_pos = raise_it(parse_field(parse_int, "cmax_pos"))
s4_ymax = raise_it(parse_field(parse_float, "ymax"))
s4_ymax_pos = raise_it(parse_field(parse_int, "ymax_pos"))
s4_smax = raise_it(parse_field(parse_float, "smax"))
s4_smax_pos = raise_it(parse_field(parse_int, "smax_pos"))
s4_smean = raise_it(parse_field(parse_float, "smean"))
s4_d = raise_it(parse_field(parse_float, "d"))
s4_decision = raise_it(parse_field(parse_bool("Y", "N"), "decision"))
s4_dmax_cut = raise_it(parse_field(parse_float, "dmax_cut"))
s4_networks_used = raise_it(parse_field(
    is_one_of(["SignalP-noTM", "SignalP-TM"]),
    "networks_used"
))


class SignalP4(Analysis, GFFAble):

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
            s4_name(sline[0]),
            s4_cmax(sline[1]),
            s4_cmax_pos(sline[2]),
            s4_ymax(sline[3]),
            s4_ymax_pos(sline[4]),
            s4_smax(sline[5]),
            s4_smax_pos(sline[6]),
            s4_smean(sline[7]),
            s4_d(sline[8]),
            s4_decision(sline[9]),
            s4_dmax_cut(sline[10]),
            s4_networks_used(sline[11]),
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

            except (LineParseError, FieldParseError) as e:
                raise e.as_parse_error(line=i).add_filename_from_handle(handle)
        return

    def as_gff(
        self,
        software_version: Optional[str] = None,
        database_version: Optional[str] = None,
        keep_all: bool = False,
        id_index: int = 1,
    ) -> Iterator[GFFRecord]:
        if not self.decision:
            return

        # d_decision = prediction of issecreted.
        # ymax = first aa of mature peptide
        attr = GFFAttributes(custom={
            "cmax": str(self.cmax),
            "cmax_pos": str(self.cmax_pos),
            "ymax": str(self.ymax),
            "ymax_pos": str(self.ymax_pos),
            "smax": str(self.smax),
            "smax_pos": str(self.smax_pos),
            "smean": str(self.smean),
            "d": str(self.d),
            "decision": "true" if self.decision else "false",
            "dmax_cut": str(self.dmax_cut),
            "networks_used": str(self.networks_used),
        })

        yield GFFRecord(
            seqid=self.name,
            source=self.gen_source(software_version, database_version),
            type="signal_peptide",
            start=0,
            end=self.ymax_pos - 1,
            score=self.d,
            strand=Strand.PLUS,
            attributes=attr
        )
        return


s5_name = raise_it(parse_field(parse_str, "name"))
s5_prediction = raise_it(parse_field(
    is_one_of(["SP(Sec/SPI)", "LIPO(Sec/SPII)", "TAT(Tat/SPI)", "OTHER"]),
    "prediction"
))
s5_prob_signal = raise_it(parse_field(parse_float, "prob_signal"))
s5_prob_other = raise_it(parse_field(parse_float, "prob_other"))
s5_cs_pos = raise_it(parse_field(parse_str, "cs_pos"))

SP5_CS_POS_REGEX = re.compile(
    r"CS\s+pos:\s+\d+-(?P<cs>\d+)\.?\s+"
    r"[A-Za-z]+-[A-Za-z]+\.?\s+"
    r"Pr: (?P<cs_prob>[-+]?\d*\.?\d+)"
)
s5_cs_actual_pos = raise_it(parse_field(
    parse_regex(SP5_CS_POS_REGEX),
    "cs_pos"
))


class SignalP5(Analysis, GFFAble):

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
            cs_pos: Optional[str] = s5_cs_pos(sline[4])
        elif len(sline) == 4:
            cs_pos = None
        else:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 4 or 5 but got {len(sline)}"
            )

        return cls(
            s5_name(sline[0]),
            s5_prediction(sline[1]),
            s5_prob_signal(sline[2]),
            s5_prob_other(sline[3]),
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

            except (LineParseError, FieldParseError) as e:
                raise e.as_parse_error(line=i).add_filename_from_handle(handle)
        return

    def as_gff(
        self,
        software_version: Optional[str] = None,
        database_version: Optional[str] = None,
        keep_all: bool = False,
        id_index: int = 1,
    ) -> Iterator[GFFRecord]:
        if self.cs_pos is None:
            return

        # dict(cs, cs_prob)
        if self.cs_pos == "CS pos: ?. Probable protein fragment":
            return

        cs = s5_cs_actual_pos(self.cs_pos)

        # d_decision = prediction of issecreted.
        # ymax = first aa of mature peptide
        attr = GFFAttributes(custom={
            "prediction": str(self.prediction),
            "prob_signal": str(self.prob_signal),
            "prob_other": str(self.prob_other),
            "prob_cut_site": str(cs["cs_prob"]),
        })

        yield GFFRecord(
            seqid=self.name,
            source=self.gen_source(software_version, database_version),
            type="signal_peptide",
            start=0,
            end=int(cs["cs"]) - 1,
            score=self.prob_signal,
            strand=Strand.PLUS,
            attributes=attr
        )
        return


s6_name = raise_it(parse_field(parse_str, "name"))
s6_prediction = raise_it(parse_field(
    is_one_of(["OTHER", "SP", "LIPO", "TAT", "TATLIPO", "PILIN", "NO_SP"]),
    "prediction"
))
s6_prob_signal = raise_it(parse_field(parse_float, "prob_signal"))
s6_prob_other = raise_it(parse_field(parse_float, "prob_other"))
s6_cs_pos = raise_it(parse_field(parse_str, "cs_pos"))

SP6_CS_POS_REGEX = re.compile(
    r"CS\s+pos:\s+\d+-(?P<cs>\d+)\.?\s+"
    r"Pr: (?P<cs_prob>[-+]?\d*\.?\d+)"
)
s6_cs_actual_pos = raise_it(parse_field(
    parse_regex(SP6_CS_POS_REGEX),
    "cs_pos"
))


class SignalP6(Analysis, GFFAble):

    name: str
    prediction: str
    prob_signal: float
    prob_other: float
    cs_pos: Optional[str]

    columns = ["name", "prediction", "prob_signal", "prob_other", "cs_pos"]
    types = [str, str, float, float, str_or_none]
    analysis = "signalp6"
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
    def from_line(cls, line: str) -> "SignalP6":
        """ Parse a short-format signalp5 line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) == 5:
            cs_pos: Optional[str] = s6_cs_pos(sline[4])
        elif len(sline) == 4:
            cs_pos = None
        else:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 4 or 5 but got {len(sline)}"
            )

        return cls(
            s6_name(sline[0]),
            s6_prediction(sline[1]),
            s6_prob_signal(sline[3]),
            s6_prob_other(sline[2]),
            cs_pos,
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["SignalP6"]:
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

    def as_gff(
        self,
        software_version: Optional[str] = None,
        database_version: Optional[str] = None,
        keep_all: bool = False,
        id_index: int = 1,
    ) -> Iterator[GFFRecord]:

        if self.cs_pos is None:
            return

        # dict(cs, cs_prob)
        if self.cs_pos == "CS pos: ?. Probable protein fragment":
            return

        cs = s6_cs_actual_pos(self.cs_pos)

        # d_decision = prediction of issecreted.
        # ymax = first aa of mature peptide
        attr = GFFAttributes(custom={
            "prediction": str(self.prediction),
            "prob_signal": str(self.prob_signal),
            "prob_other": str(self.prob_other),
            "prob_cut_site": str(cs["cs_prob"]),
        })

        yield GFFRecord(
            seqid=self.name,
            source=self.gen_source(software_version, database_version),
            type="signal_peptide",
            start=0,
            end=int(cs["cs"]) - 1,
            score=self.prob_signal,
            strand=Strand.PLUS,
            attributes=attr
        )
        return
