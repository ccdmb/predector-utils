#!/usr/bin/env python3

import argparse
import sqlite3
import sys

from ..analyses import Analyses
from ..database import load_db, ResultsTable, ResultRow


def cli(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-r", "--drop-name",
        dest="drop_name",
        action="store_true",
        default=False,
        help="Don't insert the names"
    )

    parser.add_argument(
        "-d", "--drop-null-dbversion",
        dest="drop_null_dbversion",
        action="store_true",
        default=False,
        help="Filter out rows that require databases, but no db version set."
    )

    parser.add_argument(
        "--include",
        dest="include",
        metavar="ANALYSIS",
        nargs="+",
        type=Analyses.from_string,
        choices=list(Analyses),
        default=list(Analyses),
        help="Only include these analyses, specify multiple with spaces."
    )

    parser.add_argument(
        "--exclude",
        dest="exclude",
        metavar="ANALYSIS",
        nargs="+",
        type=Analyses.from_string,
        choices=list(Analyses),
        default=[],
        help=(
            "Exclude these analyses, specify multiple with spaces. "
            "Overrides analyses specified in --include."
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
    tab.create_result_tables()
    if not args.drop_name:
        tab.create_decoder_table()

    target_analyses = set(args.include).difference(args.exclude)

    tab.insert_results(
        ResultRow.from_file(
            args.results,
            drop_name=args.drop_name,
            drop_null_dbversion=args.drop_null_dbversion,
            target_analyses=target_analyses
        ),
        insert_names=not args.drop_name
    )

    tab.create_result_index()

    if not args.drop_name:
        tab.create_decoder_index()
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
