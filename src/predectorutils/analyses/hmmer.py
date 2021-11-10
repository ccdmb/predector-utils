#!/usr/bin/env python3

from typing import Optional, Union
from typing import TextIO
from typing import Iterator

from predectorutils.gff import (
    GFFRecord,
    GFFAttributes,
    Target,
    Strand,
)

from predectorutils.analyses.base import Analysis, GFFAble
from predectorutils.analyses.base import str_or_none
from predectorutils.parsers import (
    FieldParseError,
    LineParseError,
    ValueParseError,
    raise_it,
    parse_field,
    parse_str,
    parse_float,
    parse_int,
    MULTISPACE_REGEX,
)


def split_hmm(s: str) -> Union[ValueParseError, str]:
    s1 = parse_str(s)
    if isinstance(s1, ValueParseError):
        return s1
    else:
        return s.rsplit(".hmm", maxsplit=1)[0]


hm_name = raise_it(parse_field(parse_str, "name"))  # query name
hm_hmm = raise_it(parse_field(split_hmm, "hmm"))  # target name
hm_hmm_len = raise_it(parse_field(parse_int, "hmm_len"))  # tlen
hm_query_len = raise_it(parse_field(parse_int, "query_len"))  # qlen
hm_full_evalue = raise_it(parse_field(parse_float, "full_evalue"))
hm_full_score = raise_it(parse_field(parse_float, "full_score"))
hm_full_bias = raise_it(parse_field(parse_float, "full_bias"))
hm_nmatches = raise_it(parse_field(parse_int, "nmatches"))
hm_domain_c_evalue = raise_it(parse_field(parse_float, "domain_c_evalue"))
hm_domain_i_evalue = raise_it(parse_field(parse_float, "domain_i_evalue"))
hm_domain_score = raise_it(parse_field(parse_float, "domain_score"))
hm_domain_bias = raise_it(parse_field(parse_float, "domain_bias"))
hm_hmm_from = raise_it(parse_field(parse_int, "hmm_from"))
hm_hmm_to = raise_it(parse_field(parse_int, "hmm_to"))
hm_query_from = raise_it(parse_field(parse_int, "query_from"))
hm_query_to = raise_it(parse_field(parse_int, "query_to"))
hm_acc = raise_it(parse_field(parse_float, "acc"))


class DomTbl(Analysis, GFFAble):

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
            hm_name(sline[3]),
            hm_hmm(sline[0]),
            hm_hmm_len(sline[2]),
            hm_query_len(sline[5]),
            hm_full_evalue(sline[6]),
            hm_full_score(sline[7]),
            hm_full_bias(sline[8]),
            hm_nmatches(sline[10]),
            hm_domain_c_evalue(sline[11]),
            hm_domain_i_evalue(sline[12]),
            hm_domain_score(sline[13]),
            hm_domain_bias(sline[14]),
            hm_hmm_from(sline[15]) - 1,
            hm_hmm_to(sline[16]),
            hm_query_from(sline[17]) - 1,
            hm_query_to(sline[18]),
            hm_acc(sline[21]),
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

            except (FieldParseError, LineParseError) as e:
                raise e.as_parse_error(line=i).add_filename_from_handle(handle)
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
            target=Target(self.hmm, self.hmm_from, self.hmm_to),
            custom={
                "hmm_len": str(self.hmm_len),
                "query_len": str(self.query_len),
                "full_evalue": str(self.full_evalue),
                "full_score": str(self.full_score),
                "full_bias": str(self.full_bias),
                "nmatches": str(self.nmatches),
                "domain_c_evalue": str(self.domain_c_evalue),
                "domain_i_evalue": str(self.domain_i_evalue),
                "domain_score": str(self.domain_score),
                "domain_bias": str(self.domain_bias),
                "acc": str(self.acc),
                "description": str(self.description),
            }
        )

        yield GFFRecord(
            seqid=self.query,
            source=self.gen_source(software_version, database_version),
            type="protein_hmm_match",
            start=self.query_from,
            end=self.query_to,
            score=self.domain_i_evalue,
            strand=Strand.UNSTRANDED,
            attributes=attr
        )
        return


class DBCAN(DomTbl):
    analysis = "dbcan"
    database = "DBCan"


class EffectorDB(DomTbl):
    analysis = "effectorsearch"
    database = "effectordb"
