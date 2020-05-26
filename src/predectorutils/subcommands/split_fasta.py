#!/usr/bin/env python3

import os
import argparse

from typing import List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def cli(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "infile",
        metavar="INFILE",
        type=str,
        help="The input files to encode.",
    )

    parser.add_argument(
        "-t", "--template",
        default="{fname}.split/chunk{index:0>4}.fasta",
        type=str,
        help="Where to put the split fasta files. Will create directories.",
    )

    parser.add_argument(
        "-s", "--size",
        default=1000,
        type=int,
        help="How many seqs per chunk?."
    )

    return


def runner(args: argparse.Namespace) -> None:
    seqs = SeqIO.parse(args.infile, "fasta")

    index = 1
    chunk: List[SeqRecord] = []
    for seq in seqs:
        if len(chunk) >= args.size:
            fname = args.template.format(fname=args.infile, index=index)
            dname = os.path.dirname(fname)

            if dname not in ('', '.'):
                os.makedirs(dname, exist_ok=True)

            SeqIO.write(chunk, fname, "fasta")
            index += 1
            chunk = [seq]
        else:
            chunk.append(seq)

    if len(chunk) > 0:
        fname = args.template.format(fname=args.infile, index=index)
        dname = os.path.dirname(fname)

        if dname not in ('', '.'):
            os.makedirs(dname, exist_ok=True)

        SeqIO.write(chunk, fname, "fasta")

    return
