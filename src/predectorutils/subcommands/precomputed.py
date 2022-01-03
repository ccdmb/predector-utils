#!/usr/bin/env python3

import os
import argparse
import sqlite3

from typing import Iterator
from typing import Tuple
from typing import TextIO
from typing import List, Set, Dict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid

from predectorutils.database import load_db, TargetRow, ResultsTable, ResultRow


def cli(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        help="Where to write the precomputed ldjson results to.",
        default=None,
    )

    parser.add_argument(
        "-t", "--template",
        type=str,
        default="{analysis}.fasta",
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

    parser.add_argument(
        "db",
        type=str,
        help="Where the sqlite database is"
    )

    parser.add_argument(
        "analyses",
        type=argparse.FileType('r'),
        help=(
            "A 3 column tsv file, no header. "
            "'analysis<tab>software_version<tab>database_version'. "
            "database_version should be empty string if None."
        )
    )

    parser.add_argument(
        "infasta",
        type=argparse.FileType('r'),
        help="The fasta file to parse as input. Cannot be stdin."
    )

    return


def get_checksum(seq: SeqRecord) -> Tuple[str, str]:
    checksum = seguid(str(seq.seq))
    return seq.id, checksum


def get_checksum_to_ids(seqs: Dict[str, SeqRecord]) -> Dict[str, Set[str]]:
    d: Dict[str, Set[str]] = dict()

    for seq in seqs.values():
        id_, chk = get_checksum(seq)
        if chk in d:
            d[chk].add(id_)
        else:
            d[chk] = {id_}

    return d


def write_remaining_seqs(
    remaining: List[str],
    seqs: Dict[str, SeqRecord],
    target: TargetRow,
    checksum_to_ids: Dict[str, Set[str]],
    template: str
) -> None:
    if len(remaining) == 0:
        return

    fname = template.format(**target.as_dict())

    dname = os.path.dirname(fname)
    if dname != '':
        os.makedirs(dname, exist_ok=True)

    buf: List[SeqRecord] = []

    with open(fname, "w") as handle:
        for checksum in remaining:
            for id_ in checksum_to_ids[checksum]:
                buf.append(seqs[id_])

            if len(buf) > 10000:
                SeqIO.write(buf, handle, "fasta")
                buf = []

        if len(buf) > 0:
            SeqIO.write(buf, handle, "fasta")
    return


def write_results(
    results: Iterator[ResultRow],
    outfile: TextIO
) -> None:
    buf: List[str] = []
    for result in results:
        buf.append(result.as_str())

        if len(buf) > 50000:
            print("\n".join(buf), file=outfile)
            buf = []

    if len(buf) > 0:
        print("\n".join(buf), file=outfile)
    return


def inner(
    con: sqlite3.Connection,
    cur: sqlite3.Cursor,
    args: argparse.Namespace
) -> None:
    seqs: Dict[str, SeqRecord] = SeqIO.to_dict(
        SeqIO.parse(args.infasta, "fasta")
    )
    checksum_to_ids = get_checksum_to_ids(seqs)
    checksums = set(checksum_to_ids.keys())

    tab = ResultsTable(con, cur)
    tab.insert_checksums(checksums)

    for target in TargetRow.from_file(args.analyses):
        if (
            target.analysis.needs_database() and
            (target.database_version is None)
        ):
            continue

        if args.outfile is not None:
            local_results = tab.select_target(target, checksums=True)
            write_results(local_results, args.outfile)

        remaining_checksums = list(tab.find_remaining(target))

        write_remaining_seqs(
            remaining_checksums,
            seqs,
            target,
            checksum_to_ids,
            args.template
        )
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
