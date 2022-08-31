#!/usr/bin/env python3

import re
import sys
import argparse

from Bio import SeqIO
from predectorutils.gff import GFFRecord, GFFAttributes, Strand, Phase

REGEX = re.compile(
    r"^\s*(?P<seqid>\S+)\s+"
    r"\[\s*(?P<start>\d*)\s*-\s*(?P<end>\d*)\s*\]\s*"
    r"(?P<strand>\(REVERSE SENSE\))?"
)


def cli(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "infile",
        metavar="INFILE",
        type=argparse.FileType("r"),
        help="The input fasta file from getorf.",
    )

    parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType("w"),
        help="Where to write the GFF3 file.",
    )

    return


def runner(args: argparse.Namespace) -> None:
    seqs = SeqIO.parse(args.infile, "fasta")
    records = []

    for seq in seqs:
        match = REGEX.match(seq.description)
        forward = match.group("strand") is None
        if forward:
            strand = Strand.PLUS
        else:
            strand = Strand.MINUS

        start = int(match.group("start"))
        end = int(match.group("end"))

        if not forward:
            assert start > end
            tmp = start
            start = end
            end = tmp
            del tmp

        start -= 1

        id_ = match.group("seqid")
        seqid, _ = id_.rsplit("_", maxsplit=1)
        record = GFFRecord(
            seqid=seqid,
            source="getorf",
            type="CDS",
            start=start,
            end=end,
            score=None,
            phase=Phase.FIRST,
            strand=strand,
            attributes=GFFAttributes(id=id_)
        )
        seq.id = id_
        seq.name = id_
        seq.description = id_
        records.append(record)

    records.sort(key=lambda f: (f.seqid, f.start))

    for line in records:
        print(line, file=args.outfile)

    return
