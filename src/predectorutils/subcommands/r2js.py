#!/usr/bin/env python3

import sys
import argparse
import json
import hashlib

from typing import Any, Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid


from ..analyses import Analysis, Analyses


def cli(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "format",
        type=Analyses.from_string,
        choices=list(Analyses),
        help="The file results to parse into a line delimited JSON format."
    )

    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help="The text file to parse as input. Use '-' for stdin."
    )

    parser.add_argument(
        "infasta",
        type=argparse.FileType('r'),
        help="The fasta file used to calculate results"
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Where to write the output to. Default: stdout"
    )

    parser.add_argument(
        "--pipeline-version",
        dest="pipeline_version",
        type=str,
        default=None,
        help="The version of predector that you're running."
    )

    parser.add_argument(
        "--software-version",
        dest="software_version",
        type=str,
        default=None,
        help=(
            "The version of the software that you're using. "
            "Note that the software itself is determined from the format arg."
        ),
    )

    parser.add_argument(
        "--database-version",
        dest="database_version",
        type=str,
        default=None,
        help=(
            "The version of the database that you're searching. "
            "Note that the database itself is determined from the format arg."
        ),
    )

    return


def get_line(
    pipeline_version: Optional[str],
    software_version: Optional[str],
    database_version: Optional[str],
    analysis_type: Analyses,
    analysis: Analysis,
    checksums: dict[str, str],
    md5sums: dict[str, str],
) -> dict[Any, Any]:
    name = getattr(analysis, analysis.name_column)
    out = {
        "software": analysis.software,
        "database": analysis.database,
        "analysis": str(analysis_type),
        "checksum": checksums.get(name, None),
        "md5sum": md5sums.get(name, None),
        "data": analysis.as_dict()
    }

    if pipeline_version is not None:
        out["pipeline_version"] = pipeline_version

    if software_version is not None:
        out["software_version"] = software_version

    if database_version is not None:
        out["database_version"] = database_version

    return out


def get_checksum(seq: SeqRecord) -> tuple[str, str]:
    checksum = seguid(str(seq.seq))
    return seq.id, checksum


def get_md5sum(seq: SeqRecord) -> tuple[str, str]:
    md = hashlib.md5(str(seq.seq).encode()).hexdigest()
    return seq.id, md


def get_checksums(seqs: list[SeqRecord]) -> dict[str, str]:
    out: dict[str, str] = {}

    for seq in seqs:
        id_, checksum = get_checksum(seq)
        out[id_] = checksum
    return out


def get_md5sums(seqs: list[SeqRecord]) -> dict[str, str]:
    out: dict[str, str] = {}

    for seq in seqs:
        id_, checksum = get_md5sum(seq)
        out[id_] = checksum
    return out


def runner(args: argparse.Namespace) -> None:
    seqs = list(SeqIO.parse(args.infasta, "fasta"))
    checksums = get_checksums(seqs)
    md5sums = get_md5sums(seqs)
    analysis = args.format.get_analysis()
    for line in analysis.from_file(args.infile):
        dline = get_line(
            args.pipeline_version,
            args.software_version,
            args.database_version,
            args.format,
            line,
            checksums,
            md5sums
        )
        print(json.dumps(dline), file=args.outfile)
    return
