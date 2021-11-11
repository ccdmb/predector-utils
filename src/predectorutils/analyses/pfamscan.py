#!/usr/bin/env python3

import re
from typing import Optional
from typing import TextIO
from typing import Iterator

from predectorutils.gff import GFFRecord, GFFAttributes, Strand, Target
from predectorutils.analyses.base import Analysis, GFFAble
from predectorutils.analyses.base import str_or_none
from predectorutils.parsers import (
    FieldParseError,
    LineParseError,
    parse_field,
    raise_it,
    parse_str,
    parse_float,
    parse_int,
    parse_bool,
    parse_or_none,
    MULTISPACE_REGEX,
)

ACT_SITE_REGEX = re.compile(r"predicted_active_site\[(?P<sites>[\d,\s]+)\]$")


def parse_predicted_active_site(
    field: str,
    field_name: str = "active_site",
) -> str:
    """ """

    field = field.strip()
    if not field.startswith("predicted_active_site"):
        raise LineParseError(
            f"Invalid value: '{field}' in the column: '{field_name}'. "
            "Must have the form 'predicted_active_site[1,2,3]'."
        )

    field = field[len("predicted_active_site"):]
    sfield = (f.strip("[],; ") for f in field.split('['))
    return ';'.join(f.replace(' ', '') for f in sfield if len(f) > 0)


ps_name = raise_it(parse_field(parse_str, "name"))
ps_ali_start = raise_it(parse_field(parse_int, "ali_start"))
ps_ali_end = raise_it(parse_field(parse_int, "ali_end"))
ps_env_start = raise_it(parse_field(parse_int, "env_start"))
ps_env_end = raise_it(parse_field(parse_int, "env_end"))
ps_hmm = raise_it(parse_field(parse_str, "hmm"))
ps_hmm_name = raise_it(parse_field(parse_str, "hmm_name"))
ps_hmm_type = raise_it(parse_field(parse_str, "hmm_type"))
ps_hmm_start = raise_it(parse_field(parse_int, "hmm_start"))
ps_hmm_end = raise_it(parse_field(parse_int, "hmm_end"))
ps_hmm_len = raise_it(parse_field(parse_int, "hmm_len"))
ps_bitscore = raise_it(parse_field(parse_float, "bitscore"))
ps_evalue = raise_it(parse_field(parse_float, "evalue"))
ps_is_significant = raise_it(parse_field(
    parse_bool("1", "0"),
    "is_significant"
))

ps_clan = raise_it(parse_field(parse_or_none(parse_str, "No_clan"), "clan"))


class PfamScan(Analysis, GFFAble):

    """ """

    columns = [
        "name",
        "ali_start",
        "ali_end",
        "env_start",
        "env_end",
        "hmm",
        "hmm_name",
        "hmm_type",
        "hmm_start",
        "hmm_end",
        "hmm_len",
        "bitscore",
        "evalue",
        "is_significant",
        "clan",
        "active_sites"
    ]

    types = [
        str,
        int,
        int,
        int,
        int,
        str,
        str,
        str,
        int,
        int,
        int,
        float,
        float,
        bool,
        str_or_none,
        str_or_none
    ]
    analysis = "pfamscan"
    software = "Pfam-scan"
    database = "Pfam"

    def __init__(
        self,
        name: str,
        ali_start: int,
        ali_end: int,
        env_start: int,
        env_end: int,
        hmm: str,
        hmm_name: str,
        hmm_type: str,
        hmm_start: int,
        hmm_end: int,
        hmm_len: int,
        bitscore: float,
        evalue: float,
        is_significant: bool,
        clan: Optional[str],
        active_sites: Optional[str]
    ):
        self.name = name
        self.ali_start = ali_start
        self.ali_end = ali_end
        self.env_start = env_start
        self.env_end = env_end
        self.hmm = hmm
        self.hmm_name = hmm_name
        self.hmm_type = hmm_type
        self.hmm_start = hmm_start
        self.hmm_end = hmm_end
        self.hmm_len = hmm_len
        self.bitscore = bitscore
        self.evalue = evalue
        self.is_significant = is_significant
        self.clan = clan
        self.active_sites = active_sites
        return

    @classmethod
    def from_line(cls, line: str) -> "PfamScan":
        if line == "":
            raise LineParseError("The line was empty.")

        sline = MULTISPACE_REGEX.split(line.strip(), maxsplit=16)
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

        return cls(
            ps_name(sline[0]),
            ps_ali_start(sline[1]) - 1,
            ps_ali_end(sline[2]),
            ps_env_start(sline[3]) - 1,
            ps_env_end(sline[4]),
            ps_hmm(sline[5]),
            ps_hmm_name(sline[6]),
            ps_hmm_type(sline[7]),
            ps_hmm_start(sline[8]) - 1,
            ps_hmm_end(sline[9]),
            ps_hmm_len(sline[10]),
            ps_bitscore(sline[11]),
            ps_evalue(sline[12]),
            ps_is_significant(sline[13]),
            ps_clan(sline[14]),
            active_sites,
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["PfamScan"]:
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
        if not (keep_all or self.is_significant):
            return

        attr = GFFAttributes(
            target=Target(self.hmm, self.hmm_start, self.hmm_end),
            custom={
                "env_start": str(self.env_start),
                "env_end": str(self.env_end),
                "hmm_name": str(self.hmm_name),
                "hmm_type": str(self.hmm_type),
                "hmm_len": str(self.hmm_len),
                "bitscore": str(self.bitscore),
                "evalue": str(self.evalue),
                "is_significant": "true" if self.is_significant else "false",
            }
        )

        if self.clan is not None:
            attr.custom["clan"] = str(self.clan)

        if self.active_sites is not None:
            attr.custom["active_sites"] = list(
                a.strip()
                for a
                in self.active_sites
                .replace("]", ",")
                .replace("[", ",")
                .split(",")
                if a != ""
            )

        yield GFFRecord(
            seqid=self.name,
            source=self.gen_source(software_version, database_version),
            type="protein_hmm_match",
            start=self.ali_start,
            end=self.ali_end,
            score=self.evalue,
            strand=Strand.UNSTRANDED,
            attributes=attr
        )
        return
