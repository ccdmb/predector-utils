#!/usr/bin/env python3

from typing import Optional
from typing import TextIO
from typing import Iterator

from predectorutils.analyses.base import Analysis
from predectorutils.analyses.parsers import ParseError, LineParseError
from predectorutils.analyses.parsers import (
    parse_string_not_empty,
    parse_float,
    parse_int,
)


class MMSeqs(Analysis):

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
            parse_string_not_empty(sline[0], "query"),
            parse_string_not_empty(sline[1], "target"),
            parse_int(sline[2], "qstart"),
            parse_int(sline[3], "qend"),
            parse_int(sline[4], "qlen"),
            parse_int(sline[5], "tstart"),
            parse_int(sline[6], "tend"),
            parse_int(sline[7], "tlen"),
            parse_float(sline[8], "evalue"),
            parse_int(sline[9], "gapopen"),
            parse_float(sline[10], "pident"),
            parse_int(sline[11], "alnlen"),
            parse_float(sline[12], "raw"),
            parse_float(sline[13], "bits"),
            parse_string_not_empty(sline[14], "cigar"),
            parse_int(sline[15], "mismatch"),
            parse_float(sline[16], "qcov"),
            parse_float(sline[17], "tcov"),
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


class PHIBase(MMSeqs):
    analysis = "phibase"
    software = "MMSeqs2"
    database = "PHIBase"
