#!/usr/bin/env python3

import re
from typing import TextIO
from typing import Iterator
from typing import List, Tuple
from typing import Optional

from predectorutils.gff import (
    GFFRecord,
    GFFAttributes,
    Strand,
)
from predectorutils.analyses.base import Analysis, GFFAble
from predectorutils.parsers import (
    FieldParseError,
    LineParseError,
    parse_field,
    raise_it,
    parse_str,
    parse_int,
    parse_bool,
    MULTISPACE_REGEX
)

pb_name = raise_it(parse_field(parse_str, "name"))
pb_tm = raise_it(parse_field(parse_int, "tm"))
pb_sp = raise_it(parse_field(parse_bool("Y", "0"), "sp"))
pb_topology = raise_it(parse_field(parse_str, "topology"))


def parse_topology(string: str) -> List[Tuple[str, int, int]]:
    parts = re.findall(
        r"(?P<tag>[ncio])(?P<start>\d+)[-/](?P<end>\d+)",
        string
    )
    out = []
    for tag, start, end in parts:
        if tag == "n":
            out.append(("n_terminal_region", 0, int(start) - 1))
            out.append((
                "central_hydrophobic_region_of_signal_peptide",
                int(start) - 1,
                int(end)
            ))
        elif tag == "c":
            assert out[-1][0] == "central_hydrophobic_region_of_signal_peptide"
            lend = int(out[-1][2])
            out.append(("c_terminal_region", lend, int(start)))
        else:
            assert tag in ("i", "o"), string
            out.append((
                "transmembrane_polypeptide_region",
                int(start) - 1,
                int(end)
            ))

    return out


class Phobius(Analysis, GFFAble):

    """ .
    """

    columns = ["name", "tm", "sp", "topology"]
    types = [str, int, bool, str]
    analysis = "phobius"
    software = "Phobius"

    def __init__(self, name: str, tm: int, sp: bool, topology: str) -> None:
        self.name = name
        self.tm = tm
        self.sp = sp
        self.topology = topology
        return

    @classmethod
    def from_line(cls, line: str) -> "Phobius":
        """ Parse a phobius line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = MULTISPACE_REGEX.split(line.strip())

        if len(sline) != 4:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 4 but got {len(sline)}"
            )

        # Sequence is mis-spelled in the output
        if sline == ["SEQENCE", "ID", "TM", "SP", "PREDICTION"]:
            raise LineParseError("The line appears to be the header line")

        return cls(
            pb_name(sline[0]),
            pb_tm(sline[1]),
            pb_sp(sline[2]),
            pb_topology(sline[3]),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["Phobius"]:
        for i, line in enumerate(handle):
            sline = line.strip()
            if sline.startswith("#"):
                continue
            # Sequence is mis-spelled in the output
            elif i == 0 and sline.startswith("SEQENCE"):
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
        id_index: int = 1
    ) -> Iterator[GFFRecord]:

        records = []
        source = self.gen_source(software_version, database_version)

        for (type_, start, end) in parse_topology(self.topology):
            records.append(GFFRecord(
                seqid=self.name,
                source=source,
                type=type_,
                start=start,
                end=end,
                strand=Strand.UNSTRANDED,
            ))

        if self.sp:
            gff_record = GFFRecord(
                seqid=self.name,
                source=source,
                type="signal_peptide",
                start=0,
                end=1,
                strand=Strand.UNSTRANDED,
                attributes=GFFAttributes(id=f"signal_peptide{id_index}"),
                children=[c
                          for c
                          in records
                          if c.type != "transmembrane_polypeptide_region"]
            )
            gff_record.expand_to_children()
            yield gff_record

        for r in records:
            r.update_parents()
            yield r
        return
