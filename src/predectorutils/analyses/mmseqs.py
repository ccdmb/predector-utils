#!/usr/bin/env python3

import re
from typing import TextIO
from typing import Iterator
from typing import Optional

from predectorutils.gff import (
    GFFRecord,
    GFFAttributes,
    Strand,
    Target,
    Gap,
    GapCode,
    GapElement
)
from predectorutils.analyses.base import Analysis, GFFAble
from predectorutils.parsers import (
    FieldParseError,
    LineParseError,
    parse_field,
    raise_it,
    parse_str,
    parse_float,
    parse_int,
)


mm_query = raise_it(parse_field(parse_str, "query"))
mm_target = raise_it(parse_field(parse_str, "target"))
mm_qstart = raise_it(parse_field(parse_int, "qstart"))
mm_qend = raise_it(parse_field(parse_int, "qend"))
mm_qlen = raise_it(parse_field(parse_int, "qlen"))
mm_tstart = raise_it(parse_field(parse_int, "tstart"))
mm_tend = raise_it(parse_field(parse_int, "tend"))
mm_tlen = raise_it(parse_field(parse_int, "tlen"))
mm_evalue = raise_it(parse_field(parse_float, "evalue"))
mm_gapopen = raise_it(parse_field(parse_int, "gapopen"))
mm_pident = raise_it(parse_field(parse_float, "pident"))
mm_alnlen = raise_it(parse_field(parse_int, "alnlen"))
mm_raw = raise_it(parse_field(parse_float, "raw"))
mm_bits = raise_it(parse_field(parse_float, "bits"))
mm_cigar = raise_it(parse_field(parse_str, "cigar"))
mm_mismatch = raise_it(parse_field(parse_int, "mismatch"))
mm_qcov = raise_it(parse_field(parse_float, "qcov"))
mm_tcov = raise_it(parse_field(parse_float, "tcov"))


def parse_cigar(cigar):
    output = []

    parts = re.findall(r"\d+[MIDFR]", cigar)
    for part in parts:
        num = int(part[:-1])
        code = GapCode.parse(part[-1])
        output.append(GapElement(code, num))
    return Gap(output)


class MMSeqs(Analysis, GFFAble):

    """ """
    columns = [
        "query",
        "target",
        "qstart",
        "qend",
        "qlen",
        "tstart",
        "tend",
        "tlen",
        "evalue",
        "gapopen",
        "pident",
        "alnlen",
        "raw",
        "bits",
        "cigar",
        "mismatch",
        "qcov",
        "tcov"
    ]

    types = [
        str,
        str,
        int,
        int,
        int,
        int,
        int,
        int,
        float,
        int,
        float,
        int,
        float,
        float,
        str,
        int,
        float,
        float
    ]

    software = "MMSeqs2"
    analysis = "mmseqs"
    name_column = "query"

    def __init__(
        self,
        query: str,
        target: str,
        qstart: int,
        qend: int,
        qlen: int,
        tstart: int,
        tend: int,
        tlen: int,
        evalue: float,
        gapopen: int,
        pident: float,
        alnlen: int,
        raw: float,
        bits: float,
        cigar: str,
        mismatch: int,
        qcov: float,
        tcov: float
    ):
        self.query = query
        self.target = target
        self.qstart = qstart
        self.qend = qend
        self.qlen = qlen
        self.tstart = tstart
        self.tend = tend
        self.tlen = tlen
        self.evalue = evalue
        self.gapopen = gapopen
        self.pident = pident
        self.alnlen = alnlen
        self.raw = raw
        self.bits = bits
        self.cigar = cigar
        self.mismatch = mismatch
        self.qcov = qcov
        self.tcov = tcov
        return

    @classmethod
    def from_line(cls, line: str) -> "MMSeqs":
        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t", maxsplit=17)
        if len(sline) != 18:
            # Technically because of the max_split this should be impossible.
            # the description line is allowed to have spaces.
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 18 but got {len(sline)}"
            )

        return cls(
            mm_query(sline[0]),
            mm_target(sline[1]),
            mm_qstart(sline[2]) - 1,
            mm_qend(sline[3]),
            mm_qlen(sline[4]),
            mm_tstart(sline[5]) - 1,
            mm_tend(sline[6]),
            mm_tlen(sline[7]),
            mm_evalue(sline[8]),
            mm_gapopen(sline[9]),
            mm_pident(sline[10]),
            mm_alnlen(sline[11]),
            mm_raw(sline[12]),
            mm_bits(sline[13]),
            mm_cigar(sline[14]),
            mm_mismatch(sline[15]),
            mm_qcov(sline[16]),
            mm_tcov(sline[17]),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["MMSeqs"]:
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

    def decide_significant(self, evalue=1e-5, tcov=0.5) -> bool:
        return (self.evalue <= evalue) and (self.tcov >= tcov)

    def as_gff(
        self,
        software_version: Optional[str] = None,
        database_version: Optional[str] = None,
        keep_all: bool = False,
        id_index: int = 1,
    ) -> Iterator[GFFRecord]:

        if not (keep_all or self.decide_significant()):
            return

        attr = GFFAttributes(
            target=Target(self.target, self.tstart, self.tend),
            gap=parse_cigar(self.cigar),
            custom={
                "tlen": str(self.tlen),
                "evalue": str(self.evalue),
                "gapopen": str(self.gapopen),
                "pident": str(self.pident),
                "alnlen": str(self.alnlen),
                "raw": str(self.raw),
                "bits": str(self.bits),
                "mismatch": str(self.mismatch),
                "qcov": str(self.qcov),
                "tcov": str(self.tcov),
            }
        )

        yield GFFRecord(
            seqid=self.query,
            source=self.gen_source(software_version, database_version),
            type="protein_match",
            start=self.qstart,
            end=self.qend,
            score=self.evalue,
            strand=Strand.UNSTRANDED,
            attributes=attr
        )
        return


class PHIBase(MMSeqs):
    analysis = "phibase"
    software = "MMSeqs2"
    database = "PHIBase"
