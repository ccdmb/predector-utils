#!/usr/bin/env python3

import os
import argparse
import json

from typing import Iterator
from typing import Set

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


def get_analysis(results, subdf) -> Iterator[Analysis]:
    for _, line in results.fetch_df(subdf):
        sline = line.strip()
        dline = json.loads(sline)
        cls = Analyses.from_string(dline["analysis"]).get_analysis()
        analysis = cls.from_dict(dline["data"])
        yield analysis
    return


def runner(args: argparse.Namespace) -> None:
    # This thing just keeps the contents on disk.
    # Helps save memory.
    results = IndexedResults.parse(args.infile)
    seen_analyses: Set[str] = set()

    for (analysis, _, _), subdf in results.index.groupby([
        "analysis", "software_version", "database_version"
    ]):
        if analysis in seen_analyses:
            raise ValueError(
                "You shouldn't reach this point"
            )
        else:
            seen_analyses.add(analysis)
        records = get_analysis(results, subdf)
        df = pd.DataFrame(map(lambda x: x.as_series(), records))

        fname = args.template.format(analysis=analysis.analysis)
        dname = os.path.dirname(fname)
        if dname != '':
            os.makedirs(dname, exist_ok=True)

        df.to_csv(fname, sep="\t", index=False, na_rep=".")
    return
