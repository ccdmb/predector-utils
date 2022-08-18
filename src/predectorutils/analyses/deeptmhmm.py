#!/usr/bin/env python3

import re
from typing import TypeVar
from typing import Sequence, List, Tuple
from typing import Dict
from typing import TextIO
from typing import Iterable, Iterator
from typing import Optional
from typing import Callable
from typing import Pattern

from predectorutils.gff import (
    GFFRecord,
    GFFAttributes,
    Strand
)

from predectorutils.analyses.base import Analysis, GFFAble
from predectorutils.parsers import (
    LineParseError,
    BlockParseError,
    ValueParseError,
    parse_field,
    raise_it,
    parse_str,
    parse_regex,
)

__all__ = ["DeepTMHMM"]


tm_name = raise_it(parse_field(parse_str, "name"))


T = TypeVar("T")


def parse_topology(s: str):

    if len(s) == 0:
        return s

    current_type = s[0]
    start = 0

    i = 1
    while i < len(s):
        if s[i] != current_type:
            yield (current_type, start, i)

            current_type = s[i]
            start = i

        i += 1

    yield (current_type, start, i)

    return


def convert_line_err(
    lineno: int,
    field: str,
    parser: Callable[[str], T]
) -> T:
    try:
        return parser(field)
    except LineParseError as e:
        raise e.as_block_error(lineno)


def get_line(lines: Iterator[Tuple[int, str]]) -> Tuple[int, str]:
    i, line = next(lines)

    while line.strip() == "":
        i, line = next(lines)

    return i, line.strip()


def parse_regex_line(
    lines: Iterable[Tuple[int, str]],
    regex: Pattern,
    record: Dict[str, str]
) -> None:
    i, line = get_line(iter(lines))
    rline = parse_regex(regex)(line)

    if isinstance(rline, ValueParseError):
        raise rline.as_block_error(i)
    else:
        record.update(rline)
    return


NAME_REGEX = re.compile(
    r"^>(?P<name>[^\s]+)\s*\|\s*"
    r"(?P<prediction>GLOB|SP|TM|SP+TM|TM+SP|BETA)"
)


class DeepTMHMM(Analysis, GFFAble):

    """ .
    """

    columns = ["name", "prediction", "topology"]
    types = [str, str, str]
    analysis = "deeptmhmm"
    software = "DeepTMHMM"

    def __init__(
        self,
        name: str,
        prediction: str,
        topology: str,
    ) -> None:
        self.name = name
        self.prediction = prediction
        self.topology = topology
        return

    @classmethod
    def from_block(cls, lines: Sequence[str]) -> "DeepTMHMM":
        """ Parse a deeptmhmm line as an object. """

        if not isinstance(lines, Iterable):
            ilines: Iterator[Tuple[int, str]] = enumerate(iter(lines))
        else:
            ilines = enumerate(lines)

        record: Dict[str, str] = dict()
        parse_regex_line(ilines, NAME_REGEX, record)
        next(ilines)
        i, li = next(ilines)
        record["topology"] = li.strip()
        return cls(**record)

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["DeepTMHMM"]:
        block: List[str] = []

        for i, line in enumerate(handle):
            sline = line.strip()
            if sline.startswith("#"):
                continue
            elif sline == "":
                continue
            elif sline.startswith(">") and len(block) > 0:
                try:
                    yield cls.from_block(block)
                except BlockParseError as e:
                    raise (
                        e.as_parse_error(line=i - len(block))
                        .add_filename_from_handle(handle)
                    )
                block = [sline]
            else:
                block.append(sline)

        if len(block) > 0:
            try:
                yield cls.from_block(block)
            except BlockParseError as e:
                raise (
                    e.as_parse_error(line=i - len(block))
                    .add_filename_from_handle(handle)
                )

        return

    def as_gff(
        self,
        software_version: Optional[str] = None,
        database_version: Optional[str] = None,
        keep_all: bool = False,
        id_index: int = 1
    ) -> Iterator[GFFRecord]:
        for (type_, start, end) in parse_topology(self.topology):
            if type_ in ("I", "O", "P"):
                continue

            mapp = {
                "M": "transmembrane_polypeptide_region",
                "B": "transmembrane_polypeptide_region",  # TODO
                "S": "signal_peptide",
            }
            mapp2 = {
                "M": "alpha_helix",
                "B": "beta_barrel",
            }

            if type_ in "MB":
                attr = GFFAttributes(custom={
                    "kind": mapp2[type_],
                })
            else:
                attr = None

            yield GFFRecord(
                seqid=self.name,
                source=self.gen_source(software_version, database_version),
                type=mapp[type_],
                start=start,
                end=end,
                strand=Strand.UNSTRANDED,
                attributes=attr
            )
        return
