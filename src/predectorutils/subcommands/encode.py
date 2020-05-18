#!/usr/bin/env python3

from os.path import split as psplit

import argparse

from typing import NamedTuple
from typing import Dict
from typing import Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid

from predectorutils.baseconv import IdConverter


def cli(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "outfasta",
        type=argparse.FileType('w'),
        help="Where to write the output to. Default: stdout"
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

    parser.add_argument(
        "-s", "--strip",
        default=None,
        type=str,
        help=(
            "Strip these characters from the end of the sequence. "
            "This is only respected for 'fasta' formats."
        ),
    )

    parser.add_argument(
        "-U", "--upper",
        default=False,
        action="store_true",
        help=(
            "Convert sequences to uppercase. "
            "This is only used for 'fasta' format."
        ),
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


def get_table_line(
    filename: str,
    seq: SeqRecord,
    id_conv: IdConverter,
    checksums: Dict[str, str]
) -> TableLine:
    id, checksum = get_checksum(seq)

    encoded = checksums.get(checksum, None)
    if encoded is None:
        encoded = next(id_conv)

    basename = psplit(filename)[1]

    return TableLine(encoded, basename, id, checksum)


def format_table_line(t: TableLine) -> str:
    return f"{t.encoded}\t{t.filename}\t{t.id}\t{t.checksum}"


def runner(args: argparse.Namespace) -> None:
    checksums: Dict[str, str] = dict()
    id_conv = IdConverter(prefix=args.prefix, length=args.length)

    seq_chunk = list()
    tab_chunk = list()

    for infile in args.infiles:
        seqs = SeqIO.parse(infile, "fasta")
        for i, seq in enumerate(seqs, 1):
            seq.seq = Seq(
                str(seq.seq)
                .rstrip("*")
                .upper()
                .replace("*", "X")
            )

            line = get_table_line(infile, seq, id_conv, checksums)
            tab_chunk.append(format_table_line(line))

            if line.checksum not in checksums:
                seq.id = line.encoded
                seq.name = line.encoded
                seq.description = line.encoded

                seq_chunk.append(seq.format("fasta"))

            if i % 10000 == 0:
                args.outfasta.write(''.join(seq_chunk))
                args.outmap.write('\n'.join(tab_chunk))

    if len(seq_chunk) > 0:
        args.outfasta.write(''.join(seq_chunk))

    if len(tab_chunk) > 0:
        args.outmap.write('\n'.join(tab_chunk))
    return
