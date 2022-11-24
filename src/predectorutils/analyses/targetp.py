#!/usr/bin/env python3

import re
from typing import Optional
from typing import TextIO
from typing import Iterator

from predectorutils.gff import (
    GFFRecord,
    GFFAttributes,
    Strand,
)
from predectorutils.analyses.base import Analysis, GFFAble
from predectorutils.analyses.base import float_or_none, str_or_none
from predectorutils.parsers import (
    FieldParseError,
    LineParseError,
    parse_field,
    raise_it,
    parse_str,
    parse_regex,
    parse_float,
    is_one_of
)


tp_name = raise_it(parse_field(parse_str, "name"))
tp_prediction = raise_it(parse_field(
    is_one_of(["OTHER", "noTP", "SP", "mTP", "cTP", "luTP"]),
    "prediction"
))
tp_other = raise_it(parse_field(parse_float, "OTHER"))
tp_sp = raise_it(parse_field(parse_float, "SP"))
tp_mtp = raise_it(parse_field(parse_float, "mTP"))


pl_prediction = raise_it(parse_field(
    is_one_of(["OTHER", "noTP", "SP", "mTP", "cTP", "luTP"]),
    "prediction"
))
pl_ctp = raise_it(parse_field(parse_float, "cTP"))
pl_lutp = raise_it(parse_field(parse_float, "luTP"))

CS_POS_REGEX = re.compile(
    r"CS\s+pos:\s+\d+-(?P<cs>\d+)\.?\s+"
    r"[A-Za-z]+-[A-Za-z]+\.?\s+"
    r"Pr: (?P<cs_prob>[-+]?\d*\.?\d+)"
)
cs_actual_pos = raise_it(parse_field(
    parse_regex(CS_POS_REGEX),
    "cs_pos"
))

# CS pos luTP: 90-91. ALA-SP. Pr: 0.8378; CS pos cTP: 41-42. Pr: 0.7827

PL_CS_POS_REGEX = re.compile(
    r"CS\s+"
    r"pos\s*(?P<kind>luTP|cTP|mTP|SP)?:\s+\d+-(?P<cs>\d+)\.?\s+"
    r"(?P<AA>[A-Za-z]+-[A-Za-z]+)?\.?\s*"
    r"Pr: (?P<cs_prob>[-+]?\d*\.?\d+)"
)
pl_cs_actual_pos = raise_it(parse_field(
    parse_regex(PL_CS_POS_REGEX),
    "pl_cs_pos"
))



class TargetPNonPlant(Analysis, GFFAble):

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

        prediction = tp_prediction(sline[1])
        if prediction == "noTP":
            prediction = "OTHER"

        return cls(
            tp_name(sline[0]),
            prediction,
            tp_other(sline[2]),
            tp_sp(sline[3]),
            tp_mtp(sline[4]),
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
        elif "Probable protein fragment" in self.cs_pos:
            return
        elif self.prediction not in ("SP", "mTP"):
            return

        # dict(cs, cs_prob)
        cs = cs_actual_pos(self.cs_pos)

        # d_decision = prediction of issecreted.
        # ymax = first aa of mature peptide
        attr = GFFAttributes(custom={
            "prediction": str(self.prediction),
            "prob_signal": str(self.sp),
            "prob_mitochondrial": str(self.mtp),
            "prob_other": str(self.other),
            "prob_cut_site": str(cs["cs_prob"]),
        })

        if self.prediction == "SP":
            type_ = "signal_peptide"
            prob = self.sp

        elif self.prediction == "mTP":
            type_ = "mitochondrial_targeting_signal"
            prob = self.mtp

        else:
            # Should happen
            return

        yield GFFRecord(
            seqid=self.name,
            source=self.gen_source(software_version, database_version),
            type=type_,
            start=0,
            end=int(cs["cs"]) - 1,
            score=prob,
            strand=Strand.UNSTRANDED,
            attributes=attr
        )
        return


class TargetPPlant(Analysis, GFFAble):

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
            tp_name(sline[0]),
            pl_prediction(sline[1]),
            tp_other(sline[2]),
            tp_sp(sline[3]),
            tp_mtp(sline[4]),
            pl_ctp(sline[5]),
            pl_lutp(sline[6]),
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
        elif "Probable protein fragment" in self.cs_pos:
            return

        # luTP predictions include both a cTP and luTP cutsite.
        # We need to parse both.
        cutsites = [
            pl_cs_actual_pos(c.strip())
            for c
            in self.cs_pos.split(";")
        ]

        for cs in cutsites:
            # d_decision = prediction of issecreted.
            # ymax = first aa of mature peptide
            attr = GFFAttributes(custom={
                "prediction": str(self.prediction),
                "prob_signal": str(self.sp),
                "prob_mitochondrial": str(self.mtp),
                "prob_chloroplast": str(self.ctp),
                "prob_lumen": str(self.lutp),
                "prob_other": str(self.other),
                "prob_cut_site": str(cs["cs_prob"]),
            })
            if cs.get("kind", None) is not None:
                attr.custom["kind"] = cs["kind"]

            if cs.get("kind", None) is not None:
                type_ = {
                    "CP": "signal_peptide",
                    "mTP": "mitochondrial_targeting_signal",
                    "cTP": "transit_peptide",
                    "luTP": "transit_peptide"
                }[cs["kind"]]
                prob: Optional[float] = cs.get("cs_prob", None)
            elif self.prediction == "SP":
                type_ = "signal_peptide"
                prob = self.sp

            elif self.prediction == "mTP":
                type_ = "mitochondrial_targeting_signal"
                prob = self.mtp

            elif self.prediction == "cTP":
                type_ = "transit_peptide"
                prob = self.ctp

            elif self.prediction == "luTP":
                type_ = "transit_peptide"
                prob = self.lutp

            else:
                # Shouldn't happen
                return

            yield GFFRecord(
                seqid=self.name,
                source=self.gen_source(software_version, database_version),
                type=type_,
                start=0,
                end=int(cs["cs"]) - 1,
                score=prob,
                strand=Strand.UNSTRANDED,
                attributes=attr
            )
        return
