#!/usr/bin/env python3

import argparse
import sqlite3
import sys

from predectorutils.database import load_db, ResultsTable


def cli(parser: argparse.ArgumentParser) -> None:
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
    tab = ResultsTable(con, cur)

    for row in tab.select_all():
        print(row.as_str(), file=args.outfile)
    return


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
