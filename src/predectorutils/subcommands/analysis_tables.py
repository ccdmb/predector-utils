#!/usr/bin/env python3

import os
import argparse
import json

from typing import Iterator

import pandas as pd

from predectorutils.indexedresults import IndexedResults
from predectorutils.analyses import Analysis, Analyses


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


def get_analysis(results, an_key) -> Iterator[Analysis]:
    for line in results[an_key]:
        sline = line.strip()
        dline = json.loads(sline)
        cls = Analyses.from_string(dline["analysis"]).get_analysis()
        assert cls == an_key.analysis, cls
        analysis = cls.from_dict(dline["data"])
        yield analysis
    return


def runner(args: argparse.Namespace) -> None:
    # This thing just keeps the contents on disk.
    # Helps save memory.
    results = IndexedResults.parse(args.infile)

    for analysis in results.analyses():
        records = get_analysis(results, analysis)
        df = pd.DataFrame(map(lambda x: x.as_series(), records))

        fname = args.template.format(analysis=analysis.analysis)
        dname = os.path.dirname(fname)
        if dname != '':
            os.makedirs(dname, exist_ok=True)

        df.to_csv(fname, sep="\t", index=False, na_rep=".")
    return
