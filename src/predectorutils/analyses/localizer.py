#!/usr/bin/env python3

import re
from typing import Optional
from typing import TextIO
from typing import Iterator
from typing import Tuple

from predectorutils.higher import fmap
from predectorutils.gff import GFFRecord, GFFAttributes, Strand
from predectorutils.analyses.base import Analysis, GFFAble
from predectorutils.analyses.base import (
    int_or_none,
    float_or_none,
    str_or_none
)
from predectorutils.parsers import (
    FieldParseError,
    LineParseError,
    parse_field,
    raise_it,
    parse_str,
    parse_regex
)

# This matches strings of the form "Y (0.962| 25-49)"
TP_REGEX = re.compile(
    r"^Y \((?P<prob>\d?\.?\d+)\s*\|\s*(?P<start>\d+)-(?P<end>\d+)\)$")
NUC_REGEX = re.compile(r"^Y \((?P<sigs>.+)\)$")
HEADER_REGEX = re.compile(
    r"^Identifier\s+Chloroplast\s+Mitochondria\s+Nucleus$")


def parse_tp_field(
    field: str,
    field_name: str,
) -> Tuple[bool, Optional[float], Optional[int], Optional[int]]:
    field = field.strip()

    if field == "-":
        return (False, None, None, None)

    res = parse_field(parse_regex(TP_REGEX), field_name)(field)
    if isinstance(res, FieldParseError):
        raise res

    # This should be safe because we know the regex matched.
    return (True, float(res["prob"]), int(res["start"]), int(res["end"]))


def parse_nuc_field(
    field: str,
) -> Tuple[bool, Optional[str]]:
    field = field.strip()

    if field == "-":
        return (False, None)

    res = parse_field(parse_regex(NUC_REGEX), "nucleus")(field)
    if isinstance(res, FieldParseError):
        raise res

    # This should be safe because we know the regex matched.
    return (True, res["sigs"])


class LOCALIZER(Analysis, GFFAble):

    """     """

    columns = ["name", "chloroplast_decision", "chloroplast_prob",
               "chloroplast_start", "chloroplast_end",
               "mitochondria_decision", "mitochondria_prob",
               "mitochondria_start", "mitochondria_end",
               "nucleus_decision", "nucleus_signals"]
    types = [str,
             bool, float_or_none, int_or_none, int_or_none,
             bool, float_or_none, int_or_none, int_or_none,
             bool, str_or_none]
    analysis = "localizer"
    software = "LOCALIZER"

    def __init__(
        self,
        name: str,
        chloroplast_decision: bool,
        chloroplast_prob: Optional[float],
        chloroplast_start: Optional[int],
        chloroplast_end: Optional[int],
        mitochondria_decision: bool,
        mitochondria_prob: Optional[float],
        mitochondria_start: Optional[int],
        mitochondria_end: Optional[int],
        nucleus_decision: bool,
        nucleus_signals: Optional[str],
    ) -> None:
        self.name = name
        self.chloroplast_decision = chloroplast_decision
        self.chloroplast_prob = chloroplast_prob
        self.chloroplast_start = chloroplast_start
        self.chloroplast_end = chloroplast_end
        self.mitochondria_decision = mitochondria_decision
        self.mitochondria_prob = mitochondria_prob
        self.mitochondria_start = mitochondria_start
        self.mitochondria_end = mitochondria_end
        self.nucleus_decision = nucleus_decision
        self.nucleus_signals = nucleus_signals
        return

    @classmethod
    def from_line(cls, line: str) -> "LOCALIZER":
        """ Parse an ApoplastP line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = [c.strip() for c in line.strip().split("\t")]

        if len(sline) != 4:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 4 but got {len(sline)}"
            )

        (cp, cp_prob, cp_start, cp_end) = parse_tp_field(
            sline[1],
            "chloroplast"
        )

        (mt, mt_prob, mt_start, mt_end) = parse_tp_field(
            sline[2],
            "mitochondria"
        )

        (nuc, nuc_sigs) = parse_nuc_field(sline[3])

        return cls(
            raise_it(parse_field(parse_str, "name"))(sline[0]),
            cp,
            cp_prob,
            fmap(lambda x: x - 1 + 20, cp_start),
            cp_end,
            mt,
            mt_prob,
            fmap(lambda x: x - 1 + 20, mt_start),
            mt_end,
            nuc,
            nuc_sigs
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["LOCALIZER"]:
        for i, line in enumerate(handle):
            sline = line.strip()
            if sline.startswith("#"):
                continue
            elif sline == "":
                continue
            elif HEADER_REGEX.match(sline) is not None:
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
        keep_all: bool = True,
        id_index: int = 1,
    ) -> Iterator[GFFRecord]:
        source = self.gen_source(software_version, database_version)

        if self.chloroplast_decision:
            assert self.chloroplast_start is not None
            assert self.chloroplast_end is not None

            attr = GFFAttributes(
                note=["Putative internal chloroplast localization peptide"],
                custom={
                    "prob": str(self.chloroplast_prob),
                }
            )

            yield GFFRecord(
                seqid=self.name,
                source=source,
                type="peptide_localization_signal",
                start=self.chloroplast_start,
                end=self.chloroplast_end,
                score=self.chloroplast_prob,
                strand=Strand.PLUS,
                attributes=attr
            )

        if self.mitochondria_decision:
            assert self.mitochondria_start is not None
            assert self.mitochondria_end is not None
            attr = GFFAttributes(
                note=["Putative internal mitochondrial localization peptide"],
                custom={
                    "prob": str(self.mitochondria_prob),
                }
            )

            yield GFFRecord(
                seqid=self.name,
                source=source,
                type="mitochondrial_targeting_signal",
                start=self.mitochondria_start,
                end=self.mitochondria_end,
                score=self.mitochondria_prob,
                strand=Strand.UNSTRANDED,
                attributes=attr
            )
        return
