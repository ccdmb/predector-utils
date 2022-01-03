#!/usr/bin/env python3

import os
import argparse

from typing import Iterator
from typing import Set

import sqlite3

import pandas as pd

from predectorutils.database import (
    load_db,
    ResultsTable,
    ResultRow,
    TargetRow
)


def cli(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "db",
        type=str,
        help="Where to store the database"
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

    parser.add_argument(
        "--mem",
        type=float,
        default=1.0,
        help=(
            "The amount of RAM in gibibytes to let "
            "SQLite use for cache."
        )
    )

    return


def inner(
    con: sqlite3.Connection,
    cur: sqlite3.Cursor,
    args: argparse.Namespace
) -> None:
    from ..analyses import Analyses

    tab = ResultsTable(con, cur)
    targets = list(tab.fetch_targets())

    seen: Set[Analyses] = set()
    for target in targets:
        if target.analysis in seen:
            raise ValueError(
                "There are multiple versions of the same analysis."
            )
        else:
            seen.add(target.analysis)

        records = tab.select_target(target, checksums=False)
        df = pd.DataFrame(map(lambda x: x.as_analysis().as_series(), records))

        fname = args.template.format(analysis=str(target.analysis))
        dname = os.path.dirname(fname)
        if dname != '':
            os.makedirs(dname, exist_ok=True)

        df.to_csv(fname, sep="\t", index=False, na_rep=".")


def runner(args: argparse.Namespace) -> None:
    try:
        con, cur = load_db(args.db, args.mem)
        inner(con, cur, args)
    except Exception as e:
        raise e
    finally:
        con.commit()
        con.close()
    return
