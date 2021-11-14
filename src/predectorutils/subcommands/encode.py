#!/usr/bin/env python3

from os.path import split as psplit

import re
import argparse

from typing import NamedTuple
from typing import Dict
from typing import Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid

from predectorutils.baseconv import IdConverter

INVALID_CHARS = re.compile(r"[^A-Z]", flags=re.ASCII | re.IGNORECASE)


def cli(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "outfasta",
        type=argparse.FileType('w'),
        help="Where to write the output to."
    )

    parser.add_argument(
        "outmap",
        type=argparse.FileType("w"),
        help="Where to save the id mapping file."
    )

    parser.add_argument(
        "infiles",
        metavar="INFILE",
        nargs="+",
        type=str,
        help="The input files to encode.",
    )

    parser.add_argument(
        "-l", "--length",
        default=5,
        type=int,
        help="The number of characters to use in the new id."
    )

    parser.add_argument(
        "-p", "--prefix",
        default="SR",
        type=str,
        help="The prefix to add to the beginning of the ids.",
    )

    return


class TableLine(NamedTuple):

    encoded: str
    filename: str
    id: str
    checksum: str


def get_checksum(seq: SeqRecord) -> Tuple[str, str]:
    checksum = seguid(str(seq.seq))
    return seq.id, checksum


def format_table_line(t: TableLine) -> str:
    return f"{t.encoded}\t{t.filename}\t{t.id}\t{t.checksum}"


def runner(args: argparse.Namespace) -> None:
    checksums: Dict[str, str] = dict()
    id_conv = IdConverter(prefix=args.prefix, length=args.length)

    i = 0
    j = 1
    seq_chunk = list()
    tab_chunk = list()
    for infile in args.infiles:

        seqs = SeqIO.parse(infile, "fasta")
        for seq in seqs:
            fixed_seq = (
                str(seq.seq)
                .replace("-", "")
                .rstrip("*")
                .upper()
                .replace("*", "X")
                .replace("J", "X")
                .replace("B", "X")
                .replace("Z", "X")
                .replace("U", "X")
                .replace("O", "X")
            )

            if INVALID_CHARS.match(fixed_seq) is not None:
                raise ValueError(
                    f"The sequence {seq.id} contains invalid characters."
                )

            seq.seq = Seq(fixed_seq)

            id_, checksum = get_checksum(seq)
            if checksum in checksums:
                encoded = checksums[checksum]
                new_seq = False
            else:
                encoded = id_conv.encode(i)
                checksums[checksum] = encoded
                i += 1
                new_seq = True

            line = TableLine(encoded, psplit(infile)[1], id_, checksum)
            tab_chunk.append(format_table_line(line))

            if new_seq:
                seq.id = encoded
                seq.name = encoded
                seq.description = encoded

                seq_chunk.append(seq.format("fasta").strip())

            if j % 10000 == 0:
                args.outfasta.write('\n'.join(seq_chunk) + '\n')
                args.outmap.write('\n'.join(tab_chunk) + '\n')
                seq_chunk = list()
                tab_chunk = list()

            j += 1

    if len(seq_chunk) > 0:
        args.outfasta.write('\n'.join(seq_chunk) + '\n')

    if len(tab_chunk) > 0:
        args.outmap.write('\n'.join(tab_chunk) + '\n')

    print(len(checksums))
    return
