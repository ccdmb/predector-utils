#!/usr/bin/env python3

import sys
import argparse
import json

from typing import Dict
from typing import Tuple
from typing import Any
from typing import Optional
from typing import TextIO

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid


from predectorutils import analyses


def cli(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "format",
        type=analyses.Analyses.from_string,
        choices=list(analyses.Analyses),
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
    analysis_type: analyses.Analyses,
    analysis: analyses.Analysis,
    checksums: Dict[str, str]
) -> Dict[Any, Any]:
    name = getattr(analysis, analysis.name_column)
    out = {
        "software": analysis.software,
        "database": analysis.database,
        "analysis": str(analysis_type),
        "checksum": checksums.get(name, None),
        "data": analysis.as_dict()
    }

    if pipeline_version is not None:
        out["pipeline_version"] = pipeline_version

    if software_version is not None:
        out["software_version"] = software_version

    if database_version is not None:
        out["database_version"] = database_version

    return out


def get_checksum(seq: SeqRecord) -> Tuple[str, str]:
    checksum = seguid(str(seq.seq))
    return seq.id, checksum


def get_checksums(handle: TextIO) -> Dict[str, str]:
    seqs = SeqIO.parse(handle, "fasta")
    out: Dict[str, str] = {}

    for seq in seqs:
        id_, checksum = get_checksum(seq)
        out[id_] = checksum
    return out


def runner(args: argparse.Namespace) -> None:
    checksums = get_checksums(args.infasta)
    analysis = args.format.get_analysis()
    for line in analysis.from_file(args.infile):
        dline = get_line(
            args.pipeline_version,
            args.software_version,
            args.database_version,
            args.format,
            line,
            checksums
        )
        print(json.dumps(dline), file=args.outfile)
    return
