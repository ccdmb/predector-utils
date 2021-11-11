#!/usr/bin/env python3

import re
import sys
import argparse
import json

from typing import List

from Bio import SeqIO

from ..analyses import Analysis
from ..regex import (
    FLAGS,
    KEX2_CLASSIC,
    KEX2_OUTRAM_2020,
    KEX2_OUTRAM_2021,
    RXLRLIKE,
    RegexMatcher,
)


def cli(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Where to write the output to. Default: stdout"
    )

    parser.add_argument(
        "-l", "--ldjson",
        default=False,
        action="store_true",
        help="Write the output as ldjson rather than tsv."
    )

    parser.add_argument(
        "-k", "--kind",
        default="custom",
        choices=["kex2_cutsite", "rxlr_like_motif", "custom"],
        help=(
            "Which regular expressions to search for. "
            "If custom you must specify a regular expression to --regex. "
            "Default: custom."
        )
    )

    parser.add_argument(
        "-r", "--regex",
        default=None,
        help=(
            "The regular expression to search for. "
            "Ignored if --kind is not custom."
        )
    )

    parser.add_argument(
        "infile",
        metavar="INFILE",
        type=argparse.FileType("r"),
        help="The input fasta file.",
    )

    return


def format_ldjson(a: Analysis, first: bool) -> str:
    return json.dumps(a.as_dict()) + "\n"


def format_tsv(a: Analysis, first: bool) -> str:
    return (
        a
        .as_series()
        .to_frame()
        .T
        .to_csv(sep="\t", index=False, header=first)
    )


def runner(args: argparse.Namespace) -> None:
    seqs = SeqIO.parse(args.infile, "fasta")
    output: List[str] = []

    if args.kind == "custom":
        if args.regex is None:
            raise ValueError("If --kind is custom, you must specify a regex")

        regexes = [re.compile(args.regex, flags=FLAGS)]
    elif args.kind == "kex2_cutsite":
        regexes = [
            re.compile(p) for p
            in [KEX2_CLASSIC, KEX2_OUTRAM_2020, KEX2_OUTRAM_2021]
        ]
    elif args.kind == "rxlr_like_motif":
        regexes = [re.compile(RXLRLIKE)]
    else:
        raise ValueError("It shouldn't be possible to get to this point")

    matchers = [RegexMatcher(r, args.kind) for r in regexes]

    if args.ldjson:
        formatter = format_ldjson
    else:
        formatter = format_tsv

    i = 0
    for seq in seqs:
        if (len(output) / 5000) > 1:
            print("".join(output), file=args.outfile)
            output = []

        for matcher in matchers:
            for match in matcher.find(seq):
                output.append(formatter(match, i == 0))
                i += 1

    if len(output) > 0:
        print("".join(output), file=args.outfile)
    return
