#!/usr/bin/env python3

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
    MULTISPACE_REGEX,
)


class DomTbl(Analysis):

    """ """

    columns = [
        "query",
        "hmm",
        "hmm_len",
        "query_len",
        "full_evalue",
        "full_score",
        "full_bias",
        "nmatches",
        "domain_c_evalue",
        "domain_i_evalue",
        "domain_score",
        "domain_bias",
        "hmm_from",
        "hmm_to",
        "query_from",
        "query_to",
        "acc",
        "description"
    ]

    types = [
        str,
        str,
        int,
        int,
        float,
        float,
        float,
        int,
        float,
        float,
        float,
        float,
        int,
        int,
        int,
        int,
        float,
        str_or_none
    ]
    analysis = "hmmer"
    software = "HMMER"
    name_column = "query"

    def __init__(
        self,
        query: str,
        hmm: str,
        hmm_len: int,
        query_len: int,
        full_evalue: float,
        full_score: float,
        full_bias: float,
        nmatches: int,
        domain_c_evalue: float,
        domain_i_evalue: float,
        domain_score: float,
        domain_bias: float,
        hmm_from: int,
        hmm_to: int,
        query_from: int,
        query_to: int,
        acc: float,
        description: Optional[str]
    ) -> None:
        self.query = query
        self.hmm = hmm
        self.hmm_len = hmm_len
        self.query_len = query_len
        self.full_evalue = full_evalue
        self.full_score = full_score
        self.full_bias = full_bias
        self.nmatches = nmatches
        self.domain_c_evalue = domain_c_evalue
        self.domain_i_evalue = domain_i_evalue
        self.domain_score = domain_score
        self.domain_bias = domain_bias
        self.hmm_from = hmm_from
        self.hmm_to = hmm_to
        self.query_from = query_from
        self.query_to = query_to
        self.acc = acc
        self.description = description
        return

    @classmethod
    def from_line(cls, line: str) -> "DomTbl":
        if line == "":
            raise LineParseError("The line was empty.")

        sline = MULTISPACE_REGEX.split(line.strip(), maxsplit=22)
        if len(sline) != 22 and len(sline) != 23:
            # Technically because of the max_split this should be impossible.
            # the description line is allowed to have spaces.
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 22 or 23 but got {len(sline)}"
            )

        if len(sline) == 22:
            description: Optional[str] = None
        elif sline[22] == "-" or sline[22] == "":
            description = None
        else:
            description = sline[22]

        return cls(
            parse_string_not_empty(sline[3], "name"),  # query name
            split_hmm(parse_string_not_empty(sline[0], "hmm")),  # target name
            parse_int(sline[2], "hmm_len"),  # tlen
            parse_int(sline[5], "query_len"),  # qlen
            parse_float(sline[6], "full_evalue"),
            parse_float(sline[7], "full_score"),
            parse_float(sline[8], "full_bias"),
            parse_int(sline[10], "nmatches"),
            parse_float(sline[11], "domain_c_evalue"),
            parse_float(sline[12], "domain_i_evalue"),
            parse_float(sline[13], "domain_score"),
            parse_float(sline[14], "domain_bias"),
            parse_int(sline[15], "hmm_from"),
            parse_int(sline[16], "hmm_to"),
            parse_int(sline[17], "query_from"),
            parse_int(sline[18], "query_to"),
            parse_float(sline[21], "acc"),
            description
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["DomTbl"]:
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

    def coverage(self) -> float:
        return (self.hmm_to - self.hmm_from) / self.hmm_len

    def decide_significant(
        self,
        short_threshold=1e-3,
        long_threshold=1e-5,
        coverage_threshold=0.3,
        length_threshold=80,
    ) -> bool:
        """ These are roughly the criteria from dbcan """

        if self.coverage() < coverage_threshold:
            return False
        elif (self.query_to - self.query_from) > length_threshold:
            return self.domain_i_evalue < 1e-5
        else:
            return self.domain_i_evalue < 1e-3


class DBCAN(DomTbl):
    analysis = "dbcan"
    database = "DBCan"


def split_hmm(s: str) -> str:
    return s.rsplit(".hmm", maxsplit=1)[0]
