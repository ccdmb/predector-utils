#!/usr/bin/env python3
from __future__ import annotations

import sys
import argparse
from collections import defaultdict

from typing import TextIO
from typing import Optional
from typing import (
    List,
    DefaultDict
)

import pandas as pd
from intervaltree import IntervalTree, Interval

from predectorutils.gff import GFFRecord

NUMERIC_COLUMNS = [
    'effector_score',
    'manual_effector_score',
    'manual_secretion_score',
    'effectorp1',
    'effectorp2',
    'effectorp3_cytoplasmic',
    'effectorp3_apoplastic',
    'effectorp3_noneffector',
    'deepredeff_fungi',
    'deepredeff_oomycete',
    'apoplastp',
    'molecular_weight',
    'residue_number',
    'charge',
    'isoelectric_point',
    'aa_c_number',
    'aa_tiny_number',
    'aa_small_number',
    'aa_aliphatic_number',
    'aa_aromatic_number',
    'aa_nonpolar_number',
    'aa_charged_number',
    'aa_basic_number',
    'aa_acidic_number',
    'fykin_gap',
    'targetp_secreted_prob',
    'targetp_mitochondrial_prob',
    'deeploc_membrane',
    'deeploc_nucleus',
    'deeploc_cytoplasm',
    'deeploc_extracellular',
    'deeploc_mitochondrion',
    'deeploc_cell_membrane',
    'deeploc_endoplasmic_reticulum',
    'deeploc_plastid',
    'deeploc_golgi',
    'deeploc_lysosome',
    'deeploc_peroxisome',
    'signalp3_nn_d',
    'signalp3_hmm_s',
    'signalp4_d',
    'signalp5_prob',
    'signalp6_prob',
    'deepsig_signal_prob',
    'deepsig_transmembrane_prob',
    'deepsig_other_prob'
]


def cli(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "genes",
        type=argparse.FileType('r'),
        help="Gene GFF to use."
    )

    parser.add_argument(
        "annotations",
        type=argparse.FileType('r'),
        help="ranking table to use (i.e. the tsv output of predector)"
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Where to write the output to. Default: stdout"
    )

    parser.add_argument(
        "--target",
        type=str,
        default=None,
        nargs="+",
        help=(
            "Only output these columns into the bedgraph file. "
            "By default writes a multi bedgraph with all columns."
        )
    )

    parser.add_argument(
        "--id",
        dest="id_field",
        default="Parent",
        help=(
            "What GFF attribute field corresponds to your protein feature "
            "seqids? Default uses the Parent field. Because some fields "
            "(like Parent) can have multiple values, we'll raise an error "
            "if there is more than 1 unique value. Any CDSs missing the "
            "specified field (e.g. ID) will be skipped."
        )
    )

    parser.add_argument(
        "--reducer",
        default="max",
        choices=["min", "max", "mean", "median"],
        help=(
            "How should we combine scores if two features overlap?"
        )
    )

    return


def get_id(
    record: GFFRecord,
    id_field: str
) -> Optional[str]:

    id_ = record.attributes.get(id_field, None)
    if id_ is None:
        return None

    if isinstance(id_, list):
        id_ = list(set(id_))  # Remove duplicates
        if len(id_) > 1:
            raise ValueError(
                f"Feature had multiple ids for specified ID field "
                f"{id_field}. Got: {','.join(id_)}."
            )
        elif len(id_) == 0:
            return None
        else:
            id_ = id_.pop()

    assert isinstance(id_, str)
    return id_


def inner(  # noqa: C901
    genes: TextIO,
    annotations: TextIO,
    outfile: TextIO,
    id_field: str,
    target: Optional[List[str]],
    reducer: str
):
    gff = list(GFFRecord.from_file(genes))

    preds = pd.read_csv(
        annotations,
        sep="\t",
        na_values=".",
    )
    preds.set_index("name", inplace=True)
    if target is None:
        preds = preds[NUMERIC_COLUMNS]
    else:
        for t in target:
            if t not in NUMERIC_COLUMNS:
                raise ValueError(
                    f"ERROR: target {t} is not a numeric column."
                )
        preds = preds[target]

    # Gets the intervals of CDS entries.
    intervals: DefaultDict[str, IntervalTree] = defaultdict(IntervalTree)
    for feature in gff:
        if feature.type != "CDS":
            continue

        id_ = get_id(feature, id_field)
        intervals[feature.seqid].add(Interval(
            feature.start,
            feature.end,
            data={id_}
        ))

    # Splits any overlapping CDS entries
    # then merges them.
    for _, itree in intervals.items():
        itree.split_overlaps()
        itree.merge_equals(data_reducer=lambda a, b: a | b)

    rows = []
    for seqid, itree in intervals.items():
        for interval in itree:
            ids = [id_ for id_ in interval.data if id_ in preds.index]

            if len(ids) == 0:
                continue

            loc = pd.Series(
                [seqid, interval.begin, interval.end],
                index=["#seqid", "start", "end"]
            )
            vals = getattr(preds.loc[ids], reducer)(axis=0)
            row = pd.concat([loc, vals])
            rows.append(row)

    bedgraph = pd.DataFrame(rows)
    bedgraph.sort_values(by=["#seqid", "start"], inplace=True)
    bedgraph.to_csv(outfile, sep="\t", index=False, na_rep="NA")
    return


def runner(args: argparse.Namespace) -> None:
    try:
        inner(
            genes=args.genes,
            annotations=args.annotations,
            outfile=args.outfile,
            target=args.target,
            id_field=args.id_field,
            reducer=args.reducer,
        )
    except Exception as e:
        raise e
    return
