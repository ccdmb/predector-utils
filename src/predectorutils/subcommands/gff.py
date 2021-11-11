#!/usr/bin/env python3

import sys
import argparse
import json
from collections import defaultdict

from typing import Dict

from predectorutils.analyses import Analyses, GFFAble


def cli(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help="The ldjson file to parse as input. Use '-' for stdin."
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Where to write the output to. Default: stdout"
    )

    parser.add_argument(
        "-a", "--keep-all",
        action="store_true",
        default=False,
        help="Don't filter out insignificant hits."
    )

    return


def get_analysis(dline):
    cls = Analyses.from_string(dline["analysis"]).get_analysis()
    analysis = cls.from_dict(dline["data"])
    return analysis


def runner(args: argparse.Namespace) -> None:
    id_counter: Dict[str, int] = defaultdict(lambda: 1)

    for line in args.infile:
        sline = line.strip()
        if sline == "":
            continue

        dline = json.loads(sline)
        record = get_analysis(dline)

        if not isinstance(record, GFFAble):
            continue

        try:
            records = record.as_gff(
                software_version=dline.get("software_version", None),
                database_version=dline.get("database_version", None),
                keep_all=args.keep_all,
                id_index=id_counter[dline["analysis"]]
            )
            gline = "\n".join(map(str, records))
        except Exception as e:
            print(sline)
            raise e

        id_counter[dline["analysis"]] += 1

        if len(gline) > 0:
            args.outfile.write(gline + "\n")

    return
