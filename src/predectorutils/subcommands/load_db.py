#!/usr/bin/env python3

import argparse
import sqlite3
import sys

from predectorutils.database import load_db, ResultsTable, ResultRow


def cli(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-r", "--replace-name",
        dest="replace_name",
        action="store_true",
        default=False,
        help="Replace the analysis names with 'd'"
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


def inner(
    con: sqlite3.Connection,
    cur: sqlite3.Cursor,
    args: argparse.Namespace
) -> None:
    tab = ResultsTable(con, cur)
    tab.create_tables()

    tab.insert_results(
        ResultRow.from_file(
            args.results,
            replace_name=args.replace_name
        )
    )

    tab.index_results()
    return


def runner(args: argparse.Namespace) -> None:
    try:
        con, cur = load_db(args.db)
        inner(con, cur, args)
    except Exception as e:
        raise e
    finally:
        con.commit()
        con.close()
    return
