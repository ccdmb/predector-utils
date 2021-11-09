#!/usr/bin/env python3

import sys
import os
from os.path import basename, splitext, dirname
from copy import deepcopy

import argparse
import json
from collections import defaultdict

from typing import NamedTuple
from typing import Any
from typing import Dict, Set
from typing import List
from typing import TextIO

from predectorutils.analyses import Analysis, Analyses


def cli(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "map",
        type=argparse.FileType("r"),
        help="Where to save the id mapping file."
    )

    parser.add_argument(
        "infile",
        type=argparse.FileType("r"),
        help="The input ldjson file to decode.",
    )

    parser.add_argument(
        "-t", "--template",
        type=str,
        default="{filename}.ldjson",
        help="What to name the output files."
    )

    return


class TableLine(NamedTuple):

    encoded: str
    filename: str
    id: str
    checksum: str


def read_table_line(line: str) -> TableLine:
    sline = line.strip().split("\t")
    return TableLine(sline[0], sline[1], sline[2], sline[3])


def parse_table(handle: TextIO) -> Dict[str, List[TableLine]]:
    out = defaultdict(list)
    for line in handle:
        tl = read_table_line(line)
        out[tl.checksum].append(tl)

    return out


def make_outdir(filename: str) -> None:
    dname = dirname(filename)
    if dname != "":
        os.makedirs(dname, exist_ok=True)

    return


def get_analysis(dline: Dict[Any, Any]) -> Analysis:
    cls = Analyses.from_string(dline["analysis"]).get_analysis()
    analysis = cls.from_dict(dline["data"])
    return analysis


def runner(args: argparse.Namespace) -> None:

    previously_written: Set[str] = set()
    outchunks: Dict[str, List[str]] = defaultdict(list)
    tls = parse_table(args.map)

    for i, line in enumerate(args.infile, 1):
        sline = line.strip()
        if sline == "":
            continue

        dline = json.loads(sline)
        record = get_analysis(dline)
        try:
            record_name = dline["checksum"]
            table_lines = tls[record_name]
        except KeyError:
            print(
                f"ERROR: one of the protein names {record_name} was "
                "not in the mapping file.",
                file=sys.stderr
            )
            sys.exit(1)

        for table_line in table_lines:
            this_dline = deepcopy(dline)
            this_record = deepcopy(record)
            filename_noext = splitext(basename(table_line.filename))[0]

            filename = args.template.format(
                filename=table_line.filename,
                filename_noext=filename_noext,
            )

            setattr(this_record, this_record.name_column, table_line.id)
            this_dline["data"] = this_record.as_dict()

            outchunks[filename].append(json.dumps(this_dline))

        if i % 10000 == 0:
            for filename, chunk in outchunks.items():
                if filename in previously_written:
                    mode = "a"
                else:
                    make_outdir(filename)
                    mode = "w"

                with open(filename, mode) as handle:
                    print("\n".join(chunk), file=handle)

                previously_written.add(filename)

            outchunks = defaultdict(list)

    for filename, chunk in outchunks.items():
        if len(chunk) == 0:
            continue

        if filename in previously_written:
            mode = "a"
        else:
            make_outdir(filename)
            mode = "w"

        with open(filename, mode) as handle:
            print("\n".join(chunk), file=handle)

    return
