#!/usr/bin/env python3

import os
import argparse
import json

from typing import Iterator
from typing import Tuple
from typing import TextIO
from typing import List, Set, Dict
from typing import Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid

import pandas as pd

from predectorutils.indexedresults import IndexedResults
from predectorutils.analyses import Analyses

ANALYSIS_COLUMNS = ["analysis", "software_version", "database_version"]


def cli(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        help="Where to write the precomputed ldjson results to.",
        default="-",
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
        "-p", "--precomputed",
        type=argparse.FileType('r'),
        default=None,
        help="The ldjson file to parse as input. Use '-' for stdin."
    )

    parser.add_argument(
        "infasta",
        type=argparse.FileType('r'),
        help="The fasta file to parse as input. Cannot be stdin."
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

    return


def get_checksum(seq: SeqRecord) -> Tuple[str, str]:
    checksum = seguid(str(seq.seq))
    return seq.id, checksum


def get_id_to_checksum(infile: TextIO) -> Dict[str, str]:
    checksums: Dict[str, str] = {}

    infile.seek(0)
    seqs = SeqIO.parse(infile, "fasta")
    for seq in seqs:
        id_, checksum = get_checksum(seq)
        checksums[id_] = checksum

    infile.seek(0)
    return checksums


def get_checksum_to_ids(checksums: Dict[str, str]) -> Dict[str, List[str]]:
    d: Dict[str, List[str]] = dict()

    for id_, chk in checksums.items():
        if chk in d:
            d[chk].append(id_)
        else:
            d[chk] = [id_]

    return d


def update_analysis_line(s: str, name: str) -> str:
    s = s.strip()
    d = json.loads(s)
    cls = Analyses.from_string(d["analysis"]).get_analysis()
    analysis = cls.from_dict(d["data"])
    setattr(analysis, analysis.name_column, name)
    d["data"] = analysis.as_dict()
    return json.dumps(d)


def get_an_tuple(s: pd.Series) -> Tuple[str, str, Optional[str]]:
    analysis = s["analysis"]
    assert isinstance(analysis, str)

    software_version = s["software_version"]
    assert isinstance(software_version, str)

    database_version = s["database_version"]
    if database_version is not None:
        assert isinstance(database_version, str)
    return (analysis, software_version, database_version)


def fetch_local_precomputed(
    inhandle: TextIO,
    targets: pd.DataFrame,
    checksum_to_ids: Dict[str, List[str]],
    outhandle: TextIO
) -> Dict[Tuple[str, str, Optional[str]], Set[str]]:
    precomputed = IndexedResults.parse(inhandle)

    done_df = pd.merge(
        targets,
        precomputed,
        on=ANALYSIS_COLUMNS,
        how="inner"
    )

    done: Dict[Tuple[str, str, Optional[str]], Set[str]] = {}
    buf: List[str] = []

    for an, line in precomputed.fetch_df(done_df):
        an_tup = get_an_tuple(an)
        if an in done:
            done[an_tup].add(an["checksum"])
        else:
            done[an_tup] = {an["checksum"]}

        for id_ in checksum_to_ids.get(an["checksum"], []):
            buf_line = update_analysis_line(line, id_)
            buf.append(buf_line)

        if len(buf) > 10000:
            print("\n".join(buf), file=outhandle)
            buf = []

    if len(buf) > 0:
        print("\n".join(buf), file=outhandle)
    return done


def find_remaining(
    checksums: Set[str],
    targets: pd.DataFrame,
    done: Dict[Tuple[str, str, Optional[str]], Set[str]]
) -> Iterator[Dict[str, str]]:
    for i, row in targets.iterrows():
        an = get_an_tuple(row)
        sub_done = done.get(an, set())
        for checksum in checksums:
            if checksum in sub_done:
                continue
            else:
                d: Dict[str, str] = {}
                for k, v in zip(ANALYSIS_COLUMNS, an):
                    if v is not None:
                        d[k] = v

                d["checksum"] = checksum
                yield d
    return


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


def write_remaining_seqs(
    infasta: TextIO,
    targets: pd.DataFrame,
    id_to_checksum: Dict[str, str],
    done: Dict[Tuple[str, str, Optional[str]], Set[str]],
    template: str
):
    for _, an in targets.iterrows():
        infasta.seek(0)
        an_tup = get_an_tuple(an)
        filtered = filter_seqs_by_done(
            SeqIO.parse(infasta, "fasta"),
            an_tup,
            id_to_checksum,
            done,
        )

        fname = template.format(analysis=an.analysis)
        dname = os.path.dirname(fname)
        if dname != '':
            os.makedirs(dname, exist_ok=True)

        SeqIO.write(filtered, fname, "fasta")

    return


def runner(args: argparse.Namespace) -> None:
    id_to_checksum = get_id_to_checksum(args.infasta)
    checksum_to_ids = get_checksum_to_ids(id_to_checksum)

    targets = pd.read_csv(
        args.analyses,
        sep="\t",
        names=ANALYSIS_COLUMNS
    )

    # This thing just keeps the contents on disk.
    # Helps save memory.
    if args.precomputed is not None:
        done = fetch_local_precomputed(
            args.precomputed,
            targets,
            checksum_to_ids,
            args.outfile
        )
    else:
        done = {}

    write_remaining_seqs(
        args.infasta,
        targets,
        id_to_checksum,
        done,
        args.template
    )
    return
