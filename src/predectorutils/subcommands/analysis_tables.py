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


def select_target(
    tab: ResultsTable,
    target: TargetRow,
    checksums: bool = False,
) -> Iterator[ResultRow]:
    an = target.analysis
    assert tab.exists_table("decoder"), "no decoder table"

    anstr = str(an)
    table_name = f"results_{anstr}"

    if not tab.exists_table(f"results_{anstr}"):
        return

    if checksums:
        assert tab.exists_table("checksums"), "no checksums table"
        chk_where = "AND a.checksum IN checksums"
    else:
        chk_where = ""

    target_dict = {
        "analysis": target.analysis,
        "software_version": target.software_version
    }
    if an.needs_database():
        db_cols = (
            """
            database,
            database_version,
            """
        )
        db_where = "AND a.database_version = IFNULL(:database_version, '')"
        target_dict["database_version"] = target.database_version
    else:
        db_cols = ""
        db_where = ""

    result = tab.cur.execute(
        f"""
        SELECT
            d.name as name,
            CAST(:analysis AS analyses) as analysis,
            a.software,
            a.software_version,
            {db_cols}
            a.pipeline_version,
            a.checksum,
            a.data
        FROM {table_name} a
        INNER JOIN decoder d ON a.checksum = d.checksum
        WHERE software_version = :software_version
        {db_where}
        {chk_where}
        """,
        target_dict
    )

    for r in result:
        yield ResultRow.from_rowfactory(r)

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

        records = select_target(tab, target, checksums=False)
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
