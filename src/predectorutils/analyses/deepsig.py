#!/usr/bin/env python3

from typing import Optional
from typing import TextIO
from typing import Iterator

from predectorutils.gff import GFFRecord, GFFAttributes, Strand
from predectorutils.analyses.base import Analysis, GFFAble
from predectorutils.analyses.base import int_or_none
from predectorutils.parsers import (
    FieldParseError,
    LineParseError,
    raise_it,
    parse_field,
    parse_str,
    parse_float,
    parse_int,
    parse_or_none,
    is_one_of
)

__all__ = ["DeepSig"]


ds_name = raise_it(parse_field(parse_str, "name"))
ds_prediction = raise_it(parse_field(
    is_one_of(["SignalPeptide", "Transmembrane", "Other"]),
    "prediction"
))
ds_prob = raise_it(parse_field(parse_float, "prob"))
ds_cs_pos = raise_it(parse_field(parse_or_none(parse_int, "-"), "cs_pos"))


class DeepSig(Analysis, GFFAble):

    """     """

    columns = ["name", "prediction", "prob", "cs_pos"]
    types = [str, str, float, int_or_none]
    analysis = "deepsig"
    software = "DeepSig"

    def __init__(
        self,
        name: str,
        prediction: str,
        prob: float,
        cs_pos: Optional[int]
    ) -> None:
        self.name = name
        self.prediction = prediction
        self.prob = prob
        self.cs_pos = cs_pos
        return

    @classmethod
    def from_line(cls, line: str) -> "DeepSig":
        """ Parse a deepsig line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) != 4:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 4 but got {len(sline)}"
            )

        return cls(
            ds_name(sline[0]),
            ds_prediction(sline[1]),
            ds_prob(sline[2]),
            ds_cs_pos(sline[3]),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["DeepSig"]:
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

        if not self.prediction == "SignalPeptide":
            return

        # THis should always be true for signal peptides i think.
        assert self.cs_pos is not None

        # d_decision = prediction of issecreted.
        # ymax = first aa of mature peptide
        attr = GFFAttributes(custom={
            "prediction": str(self.prediction),
            "prob": str(self.prob),
            "cs_pos": str(self.cs_pos),
        })


        yield GFFRecord(
            seqid=self.name,
            source=self.gen_source(software_version, database_version),
            type="signal_peptide",
            start=0,
            end=self.cs_pos,
            score=self.prob,
            strand=Strand.UNSTRANDED,
            attributes=attr
        )
        return
