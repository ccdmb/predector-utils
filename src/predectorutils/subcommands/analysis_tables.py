#!/usr/bin/env python3

import os
import argparse

from typing import Iterator
from typing import Set

import sqlite3

import pandas as pd

from predectorutils.indexedresults import ResultsTable, TargetRow, ResultRow


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

    parser.add_argument(
        "--db",
        type=str,
        default=":memory:",
        help="Where to store the database"
    )

    return


def fetch_targets(
    tab: ResultsTable,
    table: str = "results"
) -> Iterator[TargetRow]:
    assert tab.exists_table(table), f"table {table} does not exist"

    result = tab.cur.execute((
        "SELECT DISTINCT analysis, software_version, database_version "
        f"FROM {table}"
    ))
    for r in result:
        yield TargetRow(*r)
    return


def inner(con: sqlite3.Connection, args: argparse.Namespace) -> None:
    # This thing just keeps the contents on disk.
    # Helps save memory.
    cur = con.cursor()

    tab = ResultsTable(con, cur)
    tab.create_tables()
    tab.insert_results(ResultRow.from_file(args.infile))

    targets = list(fetch_targets(tab, "results"))

    seen: Set[str] = set()
    for target in targets:
        if target.analysis in seen:
            raise ValueError(
                "There are multiple versions of the same analysis."
            )
        else:
            seen.add(target.analysis)

        records = tab.select_target(target, "results")
        df = pd.DataFrame(map(lambda x: x.as_analysis().as_series(), records))

        fname = args.template.format(analysis=target.analysis)
        dname = os.path.dirname(fname)
        if dname != '':
            os.makedirs(dname, exist_ok=True)

        df.to_csv(fname, sep="\t", index=False, na_rep=".")


def runner(args: argparse.Namespace) -> None:
    con = sqlite3.connect(args.db)
    try:
        inner(con, args)
    except Exception as e:
        raise e
    finally:
        con.commit()
        con.close()
    return
