#!/usr/bin/env python3

import os
import argparse
import json
from collections import defaultdict

import pandas as pd

from predectorutils.analyses import Analyses


def cli(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help="The ldjson file to parse as input. Use '-' for stdin."
    )

    parser.add_argument(
        "-t", "--template",
        type=str,
        default="{analysis}.tsv",
        help=(
            "A template for the output filenames. Can use python `.format` "
            "style variable analysis. Directories will be created."
        )
    )

    return


def get_analysis(dline):
    cls = Analyses.from_string(dline["analysis"]).get_analysis()
    analysis = cls.from_dict(dline["data"])
    return analysis


def runner(args: argparse.Namespace) -> None:
    records = defaultdict(list)

    for line in args.infile:
        sline = line.strip()
        if sline == "":
            continue

        dline = json.loads(sline)
        record = get_analysis(dline)
        records[dline["analysis"]].append(record)

    for key in list(records.keys()):
        df = pd.DataFrame(map(lambda x: x.as_series(), records.pop(key, [])))

        fname = args.template.format(analysis=key)
        dname = os.path.dirname(fname)
        if dname != '':
            os.makedirs(dname, exist_ok=True)

        df.to_csv(fname, sep="\t", index=False, na_rep=".")
    return
