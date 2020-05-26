#!/usr/bin/env python3

import sys
import argparse
import json

from typing import Dict
from typing import Any
from typing import Optional

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
) -> Dict[Any, Any]:
    out = {
        "protein_name": getattr(analysis, analysis.name_column),
        "software": analysis.software,
        "analysis": str(analysis_type),
        "data": analysis.as_dict()
    }

    if pipeline_version is not None:
        out["pipeline_version"] = pipeline_version

    if software_version is not None:
        out["software_version"] = software_version

    if database_version is not None:
        out["database_version"] = database_version

    return out


def runner(args: argparse.Namespace) -> None:
    analysis = args.format.get_analysis()
    for line in analysis.from_file(args.infile):
        dline = get_line(
            args.pipeline_version,
            args.software_version,
            args.database_version,
            args.format,
            line
        )
        print(json.dumps(dline), file=args.outfile)
    return
