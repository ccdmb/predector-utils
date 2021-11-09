#!/usr/bin/env python3

import os
import argparse
import json

from typing import Iterator
from typing import Tuple
from typing import TextIO
from typing import List, Set, Dict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid

from predectorutils.indexedresults import AnalysisTuple, IndexedResults
from predectorutils.analyses import Analysis, Analyses


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
        "inldjson",
        type=argparse.FileType('r'),
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


def get_targets(infile: TextIO) -> Iterator[AnalysisTuple]:
    for line in infile:
        sline = line.strip().split("\n")
        an = Analyses.from_string(sline[0])
        yield AnalysisTuple(an, sline[1], sline[2])
    return


def get_precomputed_analysis(results, an_key, checksums) -> Iterator[Analysis]:
    for line in results[(an_key, )]:
        sline = line.strip()
        dline = json.loads(sline)
        cls = Analyses.from_string(dline["analysis"]).get_analysis()
        analysis = cls.from_dict(dline["data"])
        yield analysis
    return


def update_analysis_line(s: str, name: str) -> str:
    s = s.strip()
    d = json.loads(s)
    cls = Analyses.from_string(d["analysis"]).get_analysis()
    analysis = cls.from_dict(d["data"])
    setattr(analysis, analysis.name_column, name)
    d["data"] = analysis.as_dict()
    return json.dumps(d)


def get_checksum(seq: SeqRecord) -> Tuple[str, str]:
    checksum = seguid(str(seq.seq))
    return seq.id, checksum


def get_checksums(infile: TextIO) -> Dict[str, str]:
    checksums: Dict[str, str] = {}

    infile.seek(0)
    seqs = SeqIO.parse(infile, "fasta")
    for seq in seqs:
        id_, checksum = get_checksum(seq)
        checksums[id_] = checksum

    infile.seek(0)
    return checksums


def filter_seqs_by_done(
    seqs: Iterator[SeqRecord],
    an: Analyses,
    checksums: Dict[str, str],
    done: Dict[Analyses, Set[str]],
) -> Iterator[SeqRecord]:
    if an not in done:
        return seqs

    for seq in seqs:
        chk = checksums[seq.id]
        if chk not in done[an]:
            yield seq
    return


def get_inv_checksums(checksums: Dict[str, str]) -> Dict[str, List[str]]:
    d: Dict[str, List[str]] = dict()

    for id_, chk in checksums.items():
        if chk in d:
            d[chk].append(id_)
        else:
            d[chk] = [id_]

    return d


def runner(args: argparse.Namespace) -> None:
    # This thing just keeps the contents on disk.
    # Helps save memory.
    precomputed = IndexedResults.parse(args.inldjson)
    checksums = get_checksums(args.infasta)
    inv_checksums = get_inv_checksums(checksums)

    targets = get_targets(args.analyses)

    done: Dict[Analyses, Set[str]] = {}
    buf: List[str] = []
    for an in targets:
        for chk in checksums.values():
            for line in precomputed[(targets, chk)]:
                if an.analysis in done:
                    done[an.analysis].add(chk)
                else:
                    done[an.analysis] = {chk}

                for id_ in inv_checksums.get(chk, []):
                    buf_line = update_analysis_line(line, id_)
                    buf.append(buf_line)

                if len(buf) > 5000:
                    print("\n".join(buf), file=args.outfile)
                    buf = []

        if len(buf) > 0:
            print("\n".join(buf), file=args.outfile)

        args.infasta.seek(0)
        filtered = filter_seqs_by_done(
            SeqIO.parse(args.infasta, "fasta"),
            an.analysis,
            checksums,
            done,
        )

        fname = args.template.format(analysis=an.analysis)
        dname = os.path.dirname(fname)
        if dname != '':
            os.makedirs(dname, exist_ok=True)

        SeqIO.write(filtered, fname, "fasta")
    return
