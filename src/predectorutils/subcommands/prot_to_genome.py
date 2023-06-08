#!/usr/bin/env python3
from __future__ import annotations

import sys
import argparse
from collections import defaultdict

from typing import TextIO
from collections.abc import (
    Sequence,
    Iterator, Iterable
)

from ..gff import (
    Strand,
    Phase,
    GFFRecord,
)


def cli(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "genes",
        type=argparse.FileType('r'),
        help="Gene GFF to use."
    )

    parser.add_argument(
        "annotations",
        type=argparse.FileType('r'),
        help="Annotation GFF to use (i.e. the GFF3 output of predector)"
    )

    parser.add_argument(
        "-o", "--outfile",
        type=str,
        default="stdout",
        help=(
            "Where to write the output to. If using the --split parameter "
            "this becomes the prefix. Default: stdout"
        )
    )

    parser.add_argument(
        "--split",
        action="store_true",
        default=False,
        help="Output distinct GFFs for source."
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

    return


def make_polypeptide(
    cds: GFFRecord,
    derives_from: bool = False
) -> GFFRecord:
    pep = cds.copy()
    pep.source = "source"
    pep.type = "polypeptide"
    pep.phase = Phase.NOT_CDS
    id_ = pep.attributes.id
    assert id_ is not None
    pep.attributes.id = "polypeptide-" + id_

    if derives_from:
        pep.attributes.derives_from = [cds.seqid]

    pep.attributes.alias = []
    pep.attributes.name = None
    pep.attributes.note = []
    pep.attributes.dbxref = []
    pep.attributes.custom = {}
    return pep


def cds_to_polypeptide(
    cdss: Iterable[GFFRecord],
    derives_from: bool = False
) -> list[GFFRecord]:
    peps = []
    for cds in cdss:
        pep = make_polypeptide(cds, derives_from=derives_from)
        peps.append(pep)

    return [pep]


def project_plus(
    cds: Sequence[GFFRecord],
    feature: GFFRecord
) -> Iterator[tuple[int, int]]:
    """ Projects a feature in protein space to genome coordinates
    on the forward strand.

    Note that at this point the coordinates should already
    be multiplied by 3 and the correct seqid and strand should be set.
    """

    from copy import copy
    feature = copy(feature)

    cds = sorted(cds, key=lambda c: c.start)
    for c in cds:
        if (feature.start + c.start) > c.end:
            feature.start -= len(c)
            feature.end -= len(c)
        elif (feature.end + c.start) <= c.end:
            start = feature.start + c.start
            end = feature.end + c.start
            yield (start, end)
            break
        else:
            start = feature.start + c.start
            end = c.end
            yield (start, end)
            feature.start = 0
            feature.end -= len(c)
    return


def project_minus(
    cds: Sequence[GFFRecord],
    feature: GFFRecord
) -> Iterator[tuple[int, int]]:
    """ Projects a feature in protein space to genome coordinates
    on the reverse strand.

    Note that at this point the coordinates should already
    be multiplied by 3 and the correct seqid and strand should be set.
    """

    from copy import copy
    feature = copy(feature)

    cds = sorted(cds, reverse=True, key=lambda c: c.start)
    for c in cds:
        if (c.end - feature.start) < c.start:
            feature.start -= len(c)
            feature.end -= len(c)
        elif (c.end - feature.end) >= c.start:
            start = c.end - feature.end
            end = c.end - feature.start
            yield (start, end)
            break
        else:
            start = c.start
            end = c.end - feature.start
            yield (start, end)
            feature.start = 0
            feature.end -= len(c)
    return


def preprocess_protein_feat(cdss, feature):
    """ Adds the new strand and seqid to the feature
    converts AA length to codon length. """

    from copy import copy
    strand = find_strand(cdss)
    seqid = find_seqid(cdss)

    feature = copy(feature)
    feature.seqid = seqid
    feature.strand = strand

    feature.start *= 3
    feature.end *= 3

    return feature


def find_source(cdss: Iterable[GFFRecord]) -> str:
    """ Checks that there is only one source in records. """
    sources = {c.source for c in cdss}
    if len(sources) > 1:
        raise ValueError(f"Got multiple strands {sources}")

    return sources.pop()


def find_strand(cdss: Iterable[GFFRecord]) -> Strand:
    """ Checks that there is only one strand in records. """
    strands = {c.strand for c in cdss}
    if len(strands) > 1:
        raise ValueError(f"Got multiple strands {strands}")

    return strands.pop()


def find_seqid(cdss: Iterable[GFFRecord]) -> str:
    """ Checks that there is only one seqid in records. """
    seqids = {c.seqid for c in cdss}
    if len(seqids) > 1:
        raise ValueError(f"Got multiple seqids {seqids}")

    return seqids.pop()


def split_feature_by_coords(
    coords: Iterable[tuple[int, int]],
    feature: GFFRecord
) -> Iterator[GFFRecord]:
    """ Just copies a feature for each set of new coordinates. """

    from copy import copy

    cs = []
    fs = []
    an = False
    for (start, end) in coords:
        if start == end:
            continue
        elif start > end:
            an = True
        new_feature = copy(feature)
        new_feature.start = start
        new_feature.end = end
        fs.append(new_feature)
        cs.append((start, end))
        yield new_feature

    # Neither of the below should happen, but It bit me before.
    if an:
        raise ValueError(
            "Splitting the feature yielded a feature with a higher start than end.\n"
            f"Initial feature: {repr(feature)}\n"
            f"Resulting split features: {fs}"
        )
    elif len(fs) == 0:
        raise ValueError(
            "Splitting the feature resulted in no outputs.\n"
            f"Initial feature: {repr(feature)}\n"
            f"Initial coordinates features: {cs}"
        )

    return


def project_to_cds(
    cdss: Sequence[GFFRecord],
    feature: GFFRecord
) -> list[GFFRecord]:
    """ This takes multiple CDS features from which the protein
    is derived and projects a single predicted protein feature
    back onto those CDSs, matching the intron structure.
    """

    from copy import copy
    cdss = list(copy(cdss))

    assert len(cdss) >= 1, cdss
    assert all(cds.end >= cds.start for cds in cdss)
    feature = preprocess_protein_feat(cdss, feature)

    if feature.strand == Strand.PLUS:
        coords = list(project_plus(cdss, feature))
    elif feature.strand == Strand.MINUS:
        coords = list(project_minus(cdss, feature))
    else:
        raise ValueError(f"Got an unexpected strand {feature.strand}.")

    feature.attributes.gap = None
    feats = split_feature_by_coords(coords, feature)
    return list(feats)


def create_match(
    parts: Sequence[GFFRecord],
    type_: str,
    prot_id: str,
    index: int | str
):
    """ Creates a parent feature that encapsulates all
    members of 'parts'.

    Intended to create "match" features from many "match_part"s.
    """

    from copy import deepcopy
    start = min(p.start for p in parts)
    end = max(p.end for p in parts)

    seqid = find_seqid(parts)
    strand = find_strand(parts)
    score = parts[0].score

    source_ = {p.source for p in parts}
    assert len(source_) == 1
    source = source_.pop()
    del source_

    attributes = deepcopy(parts[0].attributes)
    attributes.id = prot_id + "-" + type_ + "-" + str(index)

    parent = GFFRecord(
        seqid=seqid,
        source=source,
        type=type_,
        start=start,
        end=end,
        score=score,
        strand=strand,
        phase=Phase.NOT_CDS,
        attributes=attributes
    )

    return parent


def split_protein_features(  # noqa: W0611
    prots: Sequence[GFFRecord],
    cdss: Sequence[GFFRecord],
    derives_from: bool = False,
    separate_peps: bool = True,
    regions: bool = False
):
    from copy import deepcopy
    prots = deepcopy(prots)

    if len(cdss) == 0:
        return None

    if len(prots) == 0:
        return None

    prot_id = prots[0].seqid

    out_parents = []
    out_features = []

    index = 1
    peps = cds_to_polypeptide(cdss, derives_from=derives_from)
    if regions:
        region = peps[0]
        start = min(c.start for c in peps)
        end = max(c.end for c in peps)
        region.start = start
        region.end = end
        region.type = "region"
        peps = [region]

    for prot in prots:
        split_prot = project_to_cds(cdss, prot)

        type_ = prot.type
        source = prot.source

        if type_ in ("protein_hmm_match", "protein_match"):
            if len(cdss) > 1:
                for s in split_prot:
                    s.type = "match_part"

                parents = [create_match(split_prot, type_, prot_id, index)]
                index += 1
            else:
                # parent_type = type_
                parents = []

        else:
            # parent_type = "polypeptide"
            if separate_peps:
                parents = deepcopy(peps)
                for p in parents:
                    p.attributes.id = f"{p.attributes.id}-{index}"
                    p.source = source
                index += 1

            else:
                parents = peps

            for parent in parents:
                parent.attributes.custom["kind"] = type_

        for parent in parents:
            parent.add_children(split_prot)

        for feat in split_prot:
            feat.update_parents()

        if regions and (len(parents) > 0):
            assert len(parents) == 1
            parents[0].expand_to_children()
            parents[0].shrink_to_children()

        out_parents.extend([p for p in parents if p not in out_parents])
        out_features.extend(split_prot)

    return out_parents + out_features


def get_id(
    record: GFFRecord,
    id_field: str
) -> str | None:

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


def write_gff(
    gff: Iterable[GFFRecord],
    handle: TextIO
):
    seen: set[GFFRecord] = set()
    for feature in sorted(
        gff,
        key=lambda g: (g.seqid, g.start, g.end, g.type)
    ):
        if len(feature.parents) > 0:
            continue

        for child in feature.traverse_children(sort=True):
            if child in seen:
                continue

            seen.add(child)
            child.update_parents()
            print(child, file=handle)
        print("###", file=handle)


def split_on_source(
    records: Iterable[GFFRecord]
) -> defaultdict[str, list[GFFRecord]]:
    out: defaultdict[str, list[GFFRecord]] = defaultdict(list)
    for record in records:
        out[record.source].append(record)

    return out


def inner(  # noqa: C901
    genes: TextIO,
    annotations: TextIO,
    outfile: str,
    id_field: str,
    split: bool,
):
    genes_gff = list(GFFRecord.from_file(genes, unescape=True))

    cdss: defaultdict[str, list[GFFRecord]] = defaultdict(list)
    for g in genes_gff:
        if g.type != "CDS":
            continue

        id_ = get_id(g, id_field)
        if id_ is None:
            continue
        cdss[id_].append(g)

    prots: defaultdict[str, list[GFFRecord]] = defaultdict(list)
    for prot in GFFRecord.from_file(annotations):
        if prot.seqid not in cdss:
            continue
        prots[prot.seqid].append(prot)

    mapped: list[GFFRecord] = []
    for id_, these_prots in prots.items():
        these_cdss = cdss[id_]
        feats = split_protein_features(
            these_prots,
            these_cdss,
            separate_peps=True,
            regions=True
        )
        mapped.extend(feats)

    if (not split) and (outfile in ("stdout", "-")):
        write_gff(mapped, sys.stdout)
    elif not split:
        with open(outfile, "w") as handle:
            write_gff(mapped, handle)
    else:
        assert split
        for source, split_mapped in split_on_source(mapped).items():
            with open(f"{outfile}{source}.gff3", "w") as handle:
                write_gff(split_mapped, handle)
    return


def runner(args: argparse.Namespace) -> None:
    try:
        inner(
            genes=args.genes,
            annotations=args.annotations,
            outfile=args.outfile,
            id_field=args.id_field,
            split=args.split,
        )
    except Exception as e:
        raise e
    return
