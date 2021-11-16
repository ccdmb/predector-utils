#!/usr/bin/env python3

import os
from os.path import basename, splitext, dirname

import argparse
import sqlite3

from predectorutils.indexedresults import ResultsTable, ResultRow, DecoderRow


def cli(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "map",
        type=argparse.FileType("r"),
        help="Where to save the id mapping file."
    )

    parser.add_argument(
        "infile",
        type=argparse.FileType("r"),
        help="The input ldjson file to decode.",
    )

    parser.add_argument(
        "-t", "--template",
        type=str,
        default="{filename}.ldjson",
        help="What to name the output files."
    )

    parser.add_argument(
        "--db",
        type=str,
        default=":memory:",
        help="Where to store the database"
    )

    return


def make_outdir(filename: str) -> None:
    dname = dirname(filename)
    if dname != "":
        os.makedirs(dname, exist_ok=True)

    return


def inner(con: sqlite3.Connection, args: argparse.Namespace) -> None:
    cur = con.cursor()

    tab = ResultsTable(con, cur)
    tab.create_tables()
    tab.insert_results(ResultRow.from_file(args.infile))
    tab.insert_decoder(DecoderRow.from_file(args.map))

    for fname, decoded in tab.decode():
        filename_noext = splitext(basename(fname))[0]
        filename = args.template.format(
            filename=fname,
            filename_noext=filename_noext,
        )

        first_chunk = True

        buf = []
        for line in decoded:
            buf.append(line)
            if len(buf) > 10000:
                if first_chunk:
                    make_outdir(filename)
                    mode = "w"
                    first_chunk = False
                else:
                    mode = "a"

                with open(filename, mode) as handle:
                    print("\n".join(buf), file=handle)
                buf = []

        if len(buf) > 0:
            if first_chunk:
                make_outdir(filename)
                mode = "w"
                first_chunk = False
            else:
                mode = "a"

            with open(filename, mode) as handle:
                print("\n".join(buf), file=handle)

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
