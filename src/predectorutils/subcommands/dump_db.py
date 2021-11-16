#!/usr/bin/env python3

import argparse
import sqlite3
import sys

from predectorutils.indexedresults import ResultsTable, ResultRow


def cli(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--replace-name",
        action="store_true",
        default=False,
        help="Replace the names of analyses with 'dummy'",
    )

    parser.add_argument(
        "--table",
        default="results_deduplicated",
        choices=["results", "results_deduplicated"],
        help="Which table to dump from."
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Where to write the output to. Default: stdout"
    )

    parser.add_argument(
        "db",
        type=str,
        help="The database to dump results from."
    )

    return


def inner(con: sqlite3.Connection, args: argparse.Namespace) -> None:
    cur = con.cursor()

    tab = ResultsTable(con, cur)
    for row in tab.select_all(args.table):
        if args.replace_name:
            row = row.replace_name()

        print(row.as_str(), file=args.outfile)
    return


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
