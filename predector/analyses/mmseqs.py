#!/usr/bin/env python3

import re
from typing import Optional
from typing import TextIO
from typing import Iterator

from predector.analyses.base import Analysis
from predector.analyses.base import str_or_none
from predector.analyses.parsers import ParseError, LineParseError
from predector.analyses.parsers import (
    parse_string_not_empty,
    parse_float,
    parse_int,
    parse_bool,
    MULTISPACE_REGEX,
)

ACT_SITE_REGEX = re.compile(r"predicted_active_site\[(?P<sites>[\d,\s]+)\]$")


class MMSeqs(Analysis):

    """ """
    columns = [
        "target",
        "query",
        "tstart",
        "tend",
        "tlen",
        "qstart",
        "qend",
        "qlen",
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

    def __init__(
        self,
        target: str,
        query: str,
        tstart: int,
        tend: int,
        tlen: int,
        qstart: int,
        qend: int,
        qlen: int,
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
        self.target = target
        self.query = query
        self.tstart = tstart
        self.tend = tend
        self.tlen = tlen
        self.qstart = qstart
        self.qend = qend
        self.qlen = qlen
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

        sline = line.strip().split("\t", maxsplit=16)
        if len(sline) != 15 and len(sline) != 16:
            # Technically because of the max_split this should be impossible.
            # the description line is allowed to have spaces.
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 15 or 16 but got {len(sline)}"
            )

        if len(sline) == 15:
            active_sites: Optional[str] = None
        else:
            active_sites = parse_predicted_active_site(sline[15])

        if sline[14] == "No_clan":
            clan: Optional[str] = None
        else:
            clan = parse_string_not_empty(sline[14], "clan")

        return cls(
            parse_string_not_empty(sline[0], "name"),
            parse_int(sline[1], "ali_start"),
            parse_int(sline[2], "ali_end"),
            parse_int(sline[3], "env_start"),
            parse_int(sline[4], "env_end"),
            parse_string_not_empty(sline[5], "hmm"),
            parse_string_not_empty(sline[6], "hmm_name"),
            parse_string_not_empty(sline[7], "hmm_type"),
            parse_int(sline[8], "hmm_start"),
            parse_int(sline[9], "hmm_end"),
            parse_int(sline[10], "hmm_len"),
            parse_float(sline[11], "bitscore"),
            parse_float(sline[12], "evalue"),
            parse_bool(sline[13], "is_significant", "1", "0"),
            clan,
            active_sites,
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
