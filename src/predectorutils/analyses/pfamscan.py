#!/usr/bin/env python3

import re
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
    parse_bool,
    MULTISPACE_REGEX,
)

ACT_SITE_REGEX = re.compile(r"predicted_active_site\[(?P<sites>[\d,\s]+)\]$")


class PfamScan(Analysis):

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
    def from_file(cls, handle: TextIO) -> Iterator["PfamScan"]:
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
