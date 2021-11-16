#!/usr/bin/env python3

import os
import argparse
import sqlite3

from typing import Iterator
from typing import Tuple
from typing import TextIO
from typing import List, Set, Dict
from typing import Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid

from predectorutils.analyses import Analyses
from predectorutils.indexedresults import TargetRow, ResultsTable, ResultRow


def cli(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        help="Where to write the precomputed ldjson results to.",
        default="-",
    )

    parser.add_argument(
        "-p", "--precomputed",
        type=argparse.FileType('r'),
        default=None,
        help="The ldjson to parse as precomputed input."
    )

    parser.add_argument(
        "--db",
        type=str,
        default=":memory:",
        help="Where to store the sqlite database"
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


def filter_seqs_by_done(
    seqs: Iterator[SeqRecord],
    an: Tuple[str, str, Optional[str]],
    checksums: Dict[str, str],
    done: Dict[Tuple[str, str, Optional[str]], Set[str]],
) -> Iterator[SeqRecord]:
    if an not in done:
        return seqs

    for seq in seqs:
        chk = checksums[seq.id]
        if chk not in done[an]:
            yield seq
    return


def fetch_local_precomputed(
    tab: ResultsTable,
    checksums: Set[str],
    outfile: TextIO
) -> None:
    buf: List[str] = []
    for precomp in tab.select_checksums(checksums, "subset"):
        buf.append(precomp.as_str())

        if len(buf) > 10000:
            print("\n".join(buf), file=outfile)
            buf = []

    if len(buf) > 0:
        print("\n".join(buf), file=outfile)
    return


def write_remaining_seqs(
    seqs: Dict[str, SeqRecord],
    tab: ResultsTable,
    checksums: Set[str],
    checksum_to_ids: Dict[str, Set[str]],
    template: str
) -> None:
    for an, chks in tab.find_remaining(checksums, "subset"):
        analysis, sversion, dversion = an

        fname = template.format(analysis=analysis)
        dname = os.path.dirname(fname)
        if dname != '':
            os.makedirs(dname, exist_ok=True)

        these = []
        for chk in chks:
            for id_ in checksum_to_ids[chk]:
                these.append(seqs[id_])

        if len(these) > 0:
            SeqIO.write(these, fname, "fasta")
    return


def filter_results_by_database_version(
    results: Iterator[ResultRow],
    requires_database: Set[str],
) -> Iterator[ResultRow]:
    for result in results:
        if (
            (result.analysis in requires_database) and
            (result.database_version is None)
        ):
            continue
        else:
            yield result
    return


def inner(con: sqlite3.Connection, args: argparse.Namespace) -> None:
    seqs: Dict[str, SeqRecord] = SeqIO.to_dict(
        SeqIO.parse(args.infasta, "fasta")
    )
    checksum_to_ids = get_checksum_to_ids(seqs)
    checksums = set(checksum_to_ids.keys())

    cur = con.cursor()

    tab = ResultsTable(con, cur)
    tab.create_tables()

    requires_database = {
        str(a)
        for a
        in Analyses
        if (a.get_analysis().database is not None)
    }

    if args.precomputed is not None:
        results = filter_results_by_database_version(
            ResultRow.from_file(
                args.precomputed,
                replace_name=True
            ),
            requires_database
        )
        tab.insert_results(results)

    tab.deduplicate_table()

    tab.insert_targets(TargetRow.from_file(args.analyses))

    fetch_local_precomputed(tab, checksums, args.outfile)

    # Write out remaining tasks to be done
    write_remaining_seqs(seqs, tab, checksums, checksum_to_ids, args.template)
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
