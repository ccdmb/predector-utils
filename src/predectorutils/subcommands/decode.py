#!/usr/bin/env python3

import os
from os.path import basename, splitext, dirname
from copy import deepcopy

import argparse
import json
from collections import defaultdict

from typing import NamedTuple
from typing import Dict
from typing import List
from typing import TextIO
from typing import Iterator

from predectorutils.indexedresults import IndexedResults
from predectorutils.analyses import Analyses


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
        out[tl.filename].append(tl)

    return out


def make_outdir(filename: str) -> None:
    dname = dirname(filename)
    if dname != "":
        os.makedirs(dname, exist_ok=True)

    return


def get_inv_checksums(tlines: List[TableLine]) -> Dict[str, List[str]]:
    d: Dict[str, List[str]] = dict()

    for tline in tlines:
        chk = tline.checksum
        id_ = tline.id
        if chk in d:
            d[chk].append(id_)
        else:
            d[chk] = [id_]

    return d


def get_analysis(results, tlines) -> Iterator[str]:
    chk_to_id = get_inv_checksums(tlines)
    for line in results[(None, list(chk_to_id.keys()))]:
        sline = line.strip()
        dline = json.loads(sline)
        cls = Analyses.from_string(dline["analysis"]).get_analysis()
        analysis = cls.from_dict(dline["data"])

        for new_id in chk_to_id[dline["checksum"]]:
            new_analysis = deepcopy(analysis)
            new_dline = deepcopy(dline)
            setattr(new_analysis, new_analysis.name_column, new_id)
            new_dline["data"] = new_analysis.as_dict()
            yield json.dumps(new_dline)
    return


def runner(args: argparse.Namespace) -> None:

    encoded_index = IndexedResults.parse(args.infile)
    tls = parse_table(args.map)

    for fname, tlines in tls.items():
        filename_noext = splitext(basename(fname))[0]
        filename = args.template.format(
            filename=fname,
            filename_noext=filename_noext,
        )

        first_chunk = True

        buf = []
        for line in get_analysis(encoded_index, tlines):
            buf.append(line)
            if len(buf) > 10000:
                if first_chunk:
                    make_outdir(filename)
                    mode = "w"
                    first_chunk = False
                else:
                    mode = "a"

                with open(filename, mode) as handle:
                    print("\n".join(buf), file=handle)
                buf = []

        if len(buf) > 0:
            if first_chunk:
                make_outdir(filename)
                mode = "w"
                first_chunk = False
            else:
                mode = "a"

            with open(filename, mode) as handle:
                print("\n".join(buf), file=handle)

    return
