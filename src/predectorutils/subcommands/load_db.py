#!/usr/bin/env python3

import argparse
import sqlite3
import sys

from predectorutils.indexedresults import ResultsTable, ResultRow


def cli(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--no-replace-name",
        dest="replace_name",
        action="store_false",
        default=True,
        help="Don't replace the analysis names with 'dummy'"
    )

    parser.add_argument(
        "db",
        type=str,
        help="Where to store the sqlite database"
    )

    parser.add_argument(
        "results",
        type=argparse.FileType('r'),
        default=sys.stdin,
        help="The ldjson to insert."
    )

    return


def inner(con: sqlite3.Connection, args: argparse.Namespace) -> None:
    cur = con.cursor()

    tab = ResultsTable(con, cur)
    tab.create_tables()

    tab.insert_results(
        ResultRow.from_file(
            args.results,
            replace_name=args.replace_name
        )
    )
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
