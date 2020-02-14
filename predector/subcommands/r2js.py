#!/usr/bin/env python3

import sys
import argparse
import json
from datetime import datetime

from typing import Dict
from typing import Any
from typing import Optional

from predector import analyses


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
        "-r", "--run-name",
        dest="run_name",
        type=str,
        default=None,
        help="Add a run name to the output."
    )

    parser.add_argument(
        "-s", "--session-id",
        dest="session_id",
        type=str,
        default=None,
        help="Add a session id to the output."
    )

    # Should I convert this to a proper time object?
    parser.add_argument(
        "-t", "--start",
        dest="start",
        type=datetime.fromisoformat,
        default=None,
        help=("Add the run start-time to the output. Should be in the "
              "iso format i.e. '2011-11-04 00:05:23.283'")
    )

    return


def get_line(
    run_name: Optional[str],
    session_id: Optional[str],
    start: Optional[datetime],
    protein_name: str,
    analysis_type: analyses.Analyses,
    analysis: analyses.Analysis,
) -> Dict[Any, Any]:
    out = {
        "protein_name": protein_name,
        "analysis": str(analysis_type),
        "data": analysis.as_dict()
    }

    if run_name is not None:
        out["run_name"] = run_name

    if session_id is not None:
        out["session_id"] = session_id

    if start is not None:
        out["start"] = start.isoformat(sep=' ')

    return out


def runner(args: argparse.Namespace) -> None:
    analysis = args.format.get_analysis()
    for line in analysis.from_file(args.infile):
        dline = get_line(
            args.run_name,
            args.session_id,
            args.start,
            line.name,
            args.format,
            line
        )
        print(json.dumps(dline), file=args.outfile)
    return
