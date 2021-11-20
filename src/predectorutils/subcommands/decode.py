#!/usr/bin/env python3

import os
from os.path import basename, splitext, dirname

import argparse
import sqlite3

from typing import Iterator, TextIO

from predectorutils.database import (
    load_db,
    ResultsTable,
    ResultRow,
    DecoderRow
)


def cli(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "db",
        type=str,
        help="Where the sqlite database is"
    )

    parser.add_argument(
        "map",
        type=argparse.FileType("r"),
        help="Where to save the id mapping file."
    )

    parser.add_argument(
        "-t", "--template",
        type=str,
        default="{filename}.ldjson",
        help="What to name the output files."
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


def make_outdir(filename: str) -> None:
    dname = dirname(filename)
    if dname != "":
        os.makedirs(dname, exist_ok=True)

    return


def write_results(
    handle: TextIO,
    results: Iterator[ResultRow],
) -> None:
    buf = []
    for line in results:
        buf.append(line.as_str())
        if len(buf) > 10000:
            print("\n".join(buf), file=handle)
            buf = []

    if len(buf) > 0:
        print("\n".join(buf), file=handle)
    return


def inner(
    con: sqlite3.Connection,
    cur: sqlite3.Cursor,
    args: argparse.Namespace
) -> None:
    tab = ResultsTable(con, cur)
    tab.insert_decoder(DecoderRow.from_file(args.map))

    for fname, decoded in tab.decode():
        filename_noext = splitext(basename(fname))[0]
        filename = args.template.format(
            filename=fname,
            filename_noext=filename_noext,
        )
        make_outdir(filename)
        with open(filename, "w") as handle:
            write_results(handle, decoded)

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
