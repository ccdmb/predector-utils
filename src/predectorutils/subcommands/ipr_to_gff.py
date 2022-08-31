#!/usr/bin/env python3
from __future__ import annotations

import sys
import argparse
from typing import Tuple, List
from typing import Iterable, Iterator
from typing import Optional
from typing import NamedTuple
from typing import TextIO
from typing import Callable

import xml.etree.ElementTree as ET

from predectorutils.gff import (
    GFFRecord,
    GFFAttributes,
    Strand,
    Target,
    Phase,
    attr_unescape,
)
from predectorutils.higher import fmap

NAMESPACE = "{http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5}"


def cli(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "xml",
        type=argparse.FileType('r'),
        help="Interproscan XML file."
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help=(
            "Where to write the GFF output to. Default: stdout"
        )
    )

    parser.add_argument(
        "--namespace",
        type=str,
        default=NAMESPACE,
        help="The XML dom namespace. This is mostly included for developers."
    )
    return


def get_source(
    ips_version: str,
    db: Optional[str],
    db_version: Optional[str]
) -> str:
    assert db is not None
    assert db_version is not None
    source = f"InterProScan:{ips_version}:{db}:{db_version}"
    return source


def get_tag(
    element: ET.Element,
    namespace: str
) -> str:
    assert element.tag.startswith(namespace), element
    return element.tag[len(namespace):]


class Signature(NamedTuple):

    accession: str
    name: Optional[str]
    desc: Optional[str]
    library: str
    library_version: str
    dbxrefs: List[str]
    goterms: List[str]
    kind: Optional[str]


def process_signature(
    start: ET.Element,
    itree: Iterator[Tuple[str, ET.Element]],
    namespace: str = NAMESPACE,
) -> Signature:
    accession = start.attrib["ac"]
    name: Optional[str] = start.attrib["name"]
    if name == "FAMILY NOT NAMED":
        name = None

    desc = start.attrib.get("desc")
    library = None
    library_version = None
    kind = None

    dbxrefs: List[str] = []
    goterms: List[str] = []

    event, element = next(itree)
    tag = get_tag(element, namespace)
    while not tag == "signature":
        if (event == "start") and (tag == "signature-library-release"):
            library = element.attrib["library"]
            library_version = element.attrib["version"]

        elif (event == "start") and (tag == "entry"):
            ac = element.attrib["ac"]
            if ac.startswith("IPR"):
                ac = f"InterPro:{ac}"
            dbxrefs.append(ac)

            if desc is None:
                desc = element.attrib.get("desc")
            if name is None:
                name = element.attrib.get("name")
            if kind is None:
                kind = element.attrib.get("type")

        elif (event == "start") and (tag == "pathway-xref"):
            db = element.attrib["db"]
            id_ = element.attrib["id"]
            dbxrefs.append(f"{db}:{id_}")
        elif (event == "start") and (tag == "go-xref"):
            goterms.append(element.attrib["id"])
            n = element.attrib.get("name")
            if (n is not None) and (n != ""):
                goterms.append(
                    "GO:" + element.attrib["name"].replace(" ", "_")
                )
        event, element = next(itree)
        tag = get_tag(element, namespace)

    assert library is not None
    assert library_version is not None
    return Signature(
        fmap(attr_unescape, accession),
        fmap(attr_unescape, name),
        fmap(attr_unescape, desc),
        attr_unescape(library),
        attr_unescape(library_version),
        dbxrefs,
        goterms,
        fmap(attr_unescape, kind),
    )


def process_coils(
    element: ET.Element,
    itree: Iterator[Tuple[str, ET.Element]],
    ips_version: Optional[str],
    query_id: Optional[str],
    namespace: str = NAMESPACE
) -> Iterator[GFFRecord]:
    assert ips_version is not None
    assert query_id is not None

    library = None
    library_version = None
    start = None
    end = None

    while True:
        event, element = next(itree)
        tag = get_tag(element, namespace)
        if (event == "end") and (tag == "coils-match"):
            break

        if (event == "start") and element.tag.endswith("signature"):
            signature = process_signature(element, itree)
        elif (
            (event == "start") and
            element.tag.endswith("coils-location-fragment")
        ):
            start = element.attrib["start"]
            end = element.attrib["end"]

    library = signature.library
    library_version = signature.library_version
    assert start is not None
    assert end is not None

    source = get_source(ips_version, library, library_version)
    type_ = "peptide_coil"  # SO:0100012

    rec = GFFRecord(
        seqid=query_id,
        source=source,
        type=type_,
        start=int(start) - 1,
        end=int(end),
        score=None,
        phase=Phase.NOT_CDS,
        strand=Strand.PLUS,
        attributes=GFFAttributes(note=["Coiled-coil"])
    )
    yield rec
    return


def process_hmmer3(
    element: ET.Element,
    itree: Iterator[Tuple[str, ET.Element]],
    ips_version: Optional[str],
    query_id: Optional[str],
    namespace: str = NAMESPACE
) -> Iterator[GFFRecord]:
    assert ips_version is not None
    assert query_id is not None

    library = None
    library_version = None

    seqid = query_id
    score = element.attrib["score"]
    evalue = element.attrib["evalue"]
    while True:
        event, element = next(itree)
        tag = get_tag(element, namespace)

        if (event == "end") and (tag == "hmmer3-match"):
            break
        if (event == "start") and (tag == "signature"):
            signature = process_signature(element, itree)
            accession = signature.accession
            name = signature.name
            desc = signature.desc
            library = signature.library
            library_version = signature.library_version
            dbxrefs = signature.dbxrefs
            goterms = signature.goterms
            kind = signature.kind

        elif (event == "start") and (tag == "hmmer3-location"):
            start = element.attrib["start"]
            end = element.attrib["end"]
            env_start = element.attrib["env-start"]
            env_end = element.attrib["env-end"]
            hmm_start = element.attrib["hmm-start"]
            hmm_end = element.attrib["hmm-end"]
            complete = element.attrib["hmm-bounds"] == "COMPLETE"

    source = get_source(
        ips_version,
        library,
        library_version
    )
    type_ = "protein_hmm_match"

    custom = {
        "evalue": str(evalue),
        "score": str(score),
        "accession": str(accession)
    }

    if complete:
        custom["complete"] = "true"

    if (desc is not None) and (desc.strip() != ""):
        custom["description"] = desc.strip()

    if (kind is not None) and (kind.strip() != ""):
        custom["kind"] = kind.strip()

    if (name is not None) and (name.strip() == ""):
        name = None

    if int(hmm_end) > 1:
        custom.update({
            "hmm_start": str(hmm_start),
            "hmm_end": str(hmm_end),
            "env_start": str(env_start),
            "env_end": str(env_end),
        })

        target = Target(accession, int(hmm_start) - 1, int(hmm_end))
        assert int(hmm_start) <= int(hmm_end), element
    else:
        target = None

    rec = GFFRecord(
        seqid=seqid,
        source=source,
        type=type_,
        start=int(start) - 1,
        end=int(end),
        score=float(evalue),
        phase=Phase.NOT_CDS,
        strand=Strand.PLUS,
        attributes=GFFAttributes(
            name=name,
            target=target,
            dbxref=dbxrefs,
            ontology_term=goterms,
            custom=custom
        )
    )
    assert int(start) <= int(end), str(rec)
    yield rec
    return


def process_mobidblite(
    element: ET.Element,
    itree: Iterator[Tuple[str, ET.Element]],
    ips_version: Optional[str],
    query_id: Optional[str],
    namespace: str = NAMESPACE
) -> Iterator[GFFRecord]:
    assert ips_version is not None
    assert query_id is not None

    library = None
    library_version = None
    start = None
    end = None

    while True:
        event, element = next(itree)
        tag = get_tag(element, namespace)
        if (event == "end") and (tag == "mobidblite-match"):
            break

        if (event == "start") and element.tag.endswith("signature"):
            signature = process_signature(element, itree)
            library = signature.library
            library_version = signature.library_version

        elif (
            (event == "start") and
            element.tag.endswith("mobidblite-location-fragment")
        ):
            start = element.attrib["start"]
            end = element.attrib["end"]

    assert library is not None
    assert library_version is not None
    assert start is not None
    assert end is not None

    source = get_source(ips_version, library, library_version)
    type_ = "intrinsically_unstructured_polypeptide_region"  # SO:0100012

    rec = GFFRecord(
        seqid=query_id,
        source=source,
        type=type_,
        start=int(start) - 1,
        end=int(end),
        score=None,
        phase=Phase.NOT_CDS,
        strand=Strand.PLUS,
        attributes=GFFAttributes(note=["consensus disorder prediction"])
    )
    yield rec
    return


def process_panther(
    element: ET.Element,
    itree: Iterator[Tuple[str, ET.Element]],
    ips_version: Optional[str],
    query_id: Optional[str],
    namespace: str = NAMESPACE
) -> Iterator[GFFRecord]:
    assert ips_version is not None
    assert query_id is not None

    library = None
    library_version = None

    seqid = query_id
    score = element.attrib["score"]
    evalue = element.attrib["evalue"]
    family_name = element.attrib.get("familyName")
    if family_name in ("Not available", "FAMILY NOT NAMED"):
        family_name = None

    while True:
        event, element = next(itree)
        tag = get_tag(element, namespace)

        if (event == "end") and (tag == "panther-match"):
            break

        if (event == "start") and (tag == "signature"):
            signature = process_signature(element, itree)
            accession = signature.accession
            name = signature.name
            desc = signature.desc
            library = signature.library
            library_version = signature.library_version
            dbxrefs = signature.dbxrefs
            goterms = signature.goterms
            kind = signature.kind
        elif (event == "start") and (tag == "panther-location"):
            start = element.attrib["start"]
            end = element.attrib["end"]
            env_start = element.attrib["env-start"]
            env_end = element.attrib["env-end"]
            hmm_start = element.attrib["hmm-start"]
            hmm_end = element.attrib["hmm-end"]
            complete = element.attrib["hmm-bounds"] == "COMPLETE"

    source = get_source(ips_version, library, library_version)
    type_ = "protein_hmm_match"

    custom = {
        "evalue": str(evalue),
        "score": str(score),
        "accession": str(accession)
    }

    if complete:
        custom["complete"] = "true"

    if (desc is not None) and (desc.strip() != ""):
        custom["description"] = desc.strip()

    if (kind is not None) and (kind.strip() != ""):
        custom["kind"] = kind.strip()
    else:
        custom["kind"] = "FAMILY"

    if (name is not None) and (name.strip() == ""):
        name = None

    if name is None:
        name = family_name

    if int(hmm_end) > 1:
        custom.update({
            "hmm_start": str(hmm_start),
            "hmm_end": str(hmm_end),
            "env_start": str(env_start),
            "env_end": str(env_end),
        })

        target = Target(accession, int(hmm_start) - 1, int(hmm_end))
        assert int(hmm_start) <= int(hmm_end), element
    else:
        target = None

    rec = GFFRecord(
        seqid=seqid,
        source=source,
        type=type_,
        start=int(start) - 1,
        end=int(end),
        score=float(evalue),
        phase=Phase.NOT_CDS,
        strand=Strand.PLUS,
        attributes=GFFAttributes(
            name=name,
            target=target,
            dbxref=dbxrefs,
            ontology_term=goterms,
            custom=custom
        )
    )
    yield rec
    return


def process_hmmer2(
    element: ET.Element,
    itree: Iterator[Tuple[str, ET.Element]],
    ips_version: Optional[str],
    query_id: Optional[str],
    namespace: str = NAMESPACE
) -> Iterator[GFFRecord]:
    assert ips_version is not None
    assert query_id is not None

    library = None
    library_version = None

    seqid = query_id
    score = element.attrib["score"]
    evalue = element.attrib["evalue"]
    while True:
        event, element = next(itree)
        tag = get_tag(element, namespace)

        if (event == "end") and (tag == "hmmer2-match"):
            break

        if (event == "start") and (tag == "signature"):
            signature = process_signature(element, itree)
            accession = signature.accession
            name = signature.name
            desc = signature.desc
            library = signature.library
            library_version = signature.library_version
            dbxrefs = signature.dbxrefs
            goterms = signature.goterms
            kind = signature.kind
        elif (event == "start") and (tag == "hmmer2-location"):
            start = element.attrib["start"]
            end = element.attrib["end"]
            hmm_start = element.attrib["hmm-start"]
            hmm_end = element.attrib["hmm-end"]
            complete = element.attrib["hmm-bounds"] == "COMPLETE"

    source = get_source(ips_version, library, library_version)
    type_ = "protein_hmm_match"

    custom = {
        "evalue": str(evalue),
        "score": str(score),
        "accession": str(accession)
    }

    if complete:
        custom["complete"] = "true"

    if (desc is not None) and (desc.strip() != ""):
        custom["description"] = desc.strip()

    if (kind is not None) and (kind.strip() != ""):
        custom["kind"] = kind.strip()

    if (name is not None) and (name.strip() == ""):
        name = None

    if int(hmm_end) > 1:
        custom.update({
            "hmm_start": str(hmm_start),
            "hmm_end": str(hmm_end),
        })

        target = Target(accession, int(hmm_start) - 1, int(hmm_end))
        assert int(hmm_start) <= int(hmm_end), element
    else:
        target = None

    rec = GFFRecord(
        seqid=seqid,
        source=source,
        type=type_,
        start=int(start) - 1,
        end=int(end),
        score=float(evalue),
        phase=Phase.NOT_CDS,
        strand=Strand.PLUS,
        attributes=GFFAttributes(
            name=name,
            target=target,
            dbxref=dbxrefs,
            ontology_term=goterms,
            custom=custom
        )
    )
    assert int(start) <= int(end), str(rec)
    yield rec
    return


def process_profilescan(
    element: ET.Element,
    itree: Iterator[Tuple[str, ET.Element]],
    ips_version: Optional[str],
    query_id: Optional[str],
    namespace: str = NAMESPACE
) -> Iterator[GFFRecord]:
    assert ips_version is not None
    assert query_id is not None

    library = None
    library_version = None

    seqid = query_id
    score = element.attrib.get("score")
    while True:
        event, element = next(itree)
        tag = get_tag(element, namespace)

        if (event == "end") and (tag == "profilescan-match"):
            break

        if (event == "start") and (tag == "signature"):
            signature = process_signature(element, itree)
            accession = signature.accession
            name = signature.name
            desc = signature.desc
            library = signature.library
            library_version = signature.library_version
            dbxrefs = signature.dbxrefs
            goterms = signature.goterms
            kind = signature.kind
        elif (event == "start") and (tag == "profilescan-location"):
            start = element.attrib["start"]
            end = element.attrib["end"]
            if score is None:
                score = element.attrib.get("score")

    assert library is not None
    assert library_version is not None
    source = get_source(ips_version, library, library_version)
    type_ = "protein_match"

    custom = {
        "accession": str(accession)
    }

    if score is not None:
        custom["score"] = str(score)
        score_ = float(score)
    else:
        score_ = None

    if (desc is not None) and (desc.strip() != ""):
        custom["description"] = desc.strip()

    if (kind is not None) and (kind.strip() != ""):
        custom["kind"] = kind.strip()

    if (name is not None) and (name.strip() == ""):
        name = None

    target = None

    rec = GFFRecord(
        seqid=seqid,
        source=source,
        type=type_,
        start=int(start) - 1,
        end=int(end),
        score=score_,
        phase=Phase.NOT_CDS,
        strand=Strand.PLUS,
        attributes=GFFAttributes(
            name=name,
            target=target,
            dbxref=dbxrefs,
            ontology_term=goterms,
            custom=custom
        )
    )
    yield rec
    return


class RPSBlastSite(NamedTuple):

    desc: Optional[str]
    start: Optional[str]
    end: Optional[str]
    residue: Optional[str]


def process_rpsblast_sites(
    element: ET.Element,
    itree: Iterator[Tuple[str, ET.Element]],
    namespace: str = NAMESPACE
) -> Iterator[RPSBlastSite]:
    desc = fmap(attr_unescape, element.attrib["description"])

    while True:
        event, element = next(itree)
        tag = get_tag(element, namespace)

        if (event == "end") and (tag == "rpsblast-site"):
            break

        if (event == "start") and (tag == "site-location"):
            start = element.attrib["start"]
            end = element.attrib["end"]
            residue = element.attrib["residue"]
            yield RPSBlastSite(desc, start, end, residue)
    return


def process_rpsblast(  # noqa: C901
    element: ET.Element,
    itree: Iterator[Tuple[str, ET.Element]],
    ips_version: Optional[str],
    query_id: Optional[str],
    namespace: str = NAMESPACE
) -> Iterator[GFFRecord]:
    assert ips_version is not None
    assert query_id is not None

    library = None
    library_version = None

    seqid = query_id
    sites: List[RPSBlastSite] = []

    while True:
        event, element = next(itree)
        tag = get_tag(element, namespace)

        if (event == "end") and (tag == "rpsblast-match"):
            break

        if (event == "start") and (tag == "signature"):
            signature = process_signature(element, itree)
            accession = signature.accession
            name = signature.name
            desc = signature.desc
            library = signature.library
            library_version = signature.library_version
            dbxrefs = signature.dbxrefs
            goterms = signature.goterms
            kind = signature.kind
        elif (event == "start") and (tag == "rpsblast-location"):
            start = element.attrib["start"]
            end = element.attrib["end"]
            score = element.attrib.get("score")
            evalue = element.attrib.get("evalue")
        elif (event == "start") and (tag == "rpsblast-site"):
            sites.extend(process_rpsblast_sites(element, itree, namespace))

    source = get_source(ips_version, library, library_version)
    type_ = "protein_match"

    custom = {
        "accession": str(accession),
        "name": str(name),
    }

    if score is not None:
        custom["score"] = str(score)

    if evalue is not None:
        custom["evalue"] = str(evalue)
        evalue_ = float(evalue)
    else:
        evalue_ = None

    if (desc is not None) and (desc.strip() != ""):
        custom["description"] = desc.strip()

    if (kind is not None) and (kind.strip() != ""):
        custom["kind"] = kind.strip()

    if (name is not None) and (name.strip() == ""):
        name = None

    target = None

    rec = GFFRecord(
        seqid=seqid,
        source=source,
        type=type_,
        start=int(start) - 1,
        end=int(end),
        score=evalue_,
        phase=Phase.NOT_CDS,
        strand=Strand.PLUS,
        attributes=GFFAttributes(
            name=name,
            target=target,
            dbxref=dbxrefs,
            ontology_term=goterms,
            custom=custom
        )
    )
    yield rec

    for site in sites:
        custom = {"kind": "ACT_SITE"}
        if site.desc is not None:
            custom["description"] = site.desc

        if site.residue is not None:
            custom["residue"] = site.residue

        assert site.start is not None
        assert site.end is not None

        type_ = "biochemical_region_of_peptide"  # SO:0100001
        rec = GFFRecord(
            seqid=seqid,
            source=source,
            type=type_,
            start=int(site.start) - 1,
            end=int(site.end),
            phase=Phase.NOT_CDS,
            strand=Strand.PLUS,
            attributes=GFFAttributes(
                note=[f"Part of {name}"],
                custom=custom,
            )
        )
        yield rec
    return


def process_superfamilyhmmer3(
    element: ET.Element,
    itree: Iterator[Tuple[str, ET.Element]],
    ips_version: Optional[str],
    query_id: Optional[str],
    namespace: str = NAMESPACE
) -> Iterator[GFFRecord]:
    assert ips_version is not None
    assert query_id is not None

    library = None
    library_version = None

    seqid = query_id
    evalue = element.attrib["evalue"]
    while True:
        event, element = next(itree)
        tag = get_tag(element, namespace)

        if (event == "end") and (tag == "superfamilyhmmer3-match"):
            break

        if (event == "start") and (tag == "signature"):
            signature = process_signature(element, itree)
            accession = signature.accession
            name = signature.name
            desc = signature.desc
            library = signature.library
            library_version = signature.library_version
            dbxrefs = signature.dbxrefs
            goterms = signature.goterms
            kind = signature.kind
        elif (event == "start") and (tag == "superfamilyhmmer3-location"):
            start = element.attrib["start"]
            end = element.attrib["end"]

    source = get_source(ips_version, library, library_version)
    type_ = "protein_hmm_match"

    custom = {
        "evalue": str(evalue),
        "accession": str(accession)
    }

    if (desc is not None) and (desc.strip() != ""):
        custom["description"] = desc.strip()

    if (kind is not None) and (kind.strip() != ""):
        custom["kind"] = kind.strip()

    if (name is not None) and (name.strip() == ""):
        name = None

    target = None

    rec = GFFRecord(
        seqid=seqid,
        source=source,
        type=type_,
        start=int(start) - 1,
        end=int(end),
        score=float(evalue),
        phase=Phase.NOT_CDS,
        strand=Strand.PLUS,
        attributes=GFFAttributes(
            name=name,
            target=target,
            dbxref=dbxrefs,
            ontology_term=goterms,
            custom=custom
        )
    )
    yield rec
    return


def process_fingerprints(
    element: ET.Element,
    itree: Iterator[Tuple[str, ET.Element]],
    ips_version: Optional[str],
    query_id: Optional[str],
    namespace: str = NAMESPACE
) -> Iterator[GFFRecord]:
    from copy import copy
    assert ips_version is not None
    assert query_id is not None

    library = None
    library_version = None

    seqid = query_id
    locations = []

    while True:
        event, element = next(itree)
        tag = get_tag(element, namespace)

        if (event == "end") and (tag == "fingerprints-match"):
            break

        if (event == "start") and (tag == "signature"):
            signature = process_signature(element, itree)
            accession = signature.accession
            name = signature.name
            desc = signature.desc
            library = signature.library
            library_version = signature.library_version
            dbxrefs = signature.dbxrefs
            goterms = signature.goterms
            kind = signature.kind
        elif (event == "start") and (tag == "fingerprints-location"):
            start = element.attrib["start"]
            end = element.attrib["end"]
            pvalue = element.attrib["pvalue"]
            score = element.attrib["score"]

            locations.append((start, end, pvalue, score))

    source = get_source(ips_version, library, library_version)
    type_ = "protein_match"

    custom = {
        "accession": str(accession)
    }

    if (desc is not None) and (desc.strip() != ""):
        custom["description"] = desc.strip()

    if (kind is not None) and (kind.strip() != ""):
        custom["kind"] = kind.strip()

    if (name is not None) and (name.strip() == ""):
        name = None

    target = None

    for start, end, pvalue, score in locations:
        this_custom = copy(custom)
        this_custom["pvalue"] = pvalue
        this_custom["score"] = score

        rec = GFFRecord(
            seqid=seqid,
            source=source,
            type=type_,
            start=int(start) - 1,
            end=int(end),
            score=float(pvalue),
            phase=Phase.NOT_CDS,
            strand=Strand.PLUS,
            attributes=GFFAttributes(
                name=name,
                target=target,
                dbxref=dbxrefs,
                ontology_term=goterms,
                custom=this_custom
            )
        )
        yield rec
    return


class HMMER3Site(NamedTuple):

    desc: Optional[str]
    start: Optional[str]
    end: Optional[str]
    hmm_start: Optional[str]
    hmm_end: Optional[str]
    residue: Optional[str]


def process_hmmer3_sites(
    element: ET.Element,
    itree: Iterator[Tuple[str, ET.Element]],
    namespace: str = NAMESPACE
) -> Iterator[HMMER3Site]:
    desc = fmap(attr_unescape, element.attrib.get("description"))

    while True:
        event, element = next(itree)
        tag = get_tag(element, namespace)

        if (event == "end") and (tag == "hmmer3-site"):
            break

        if (event == "start") and (tag == "site-location"):
            start = element.attrib["start"]
            end = element.attrib["end"]
            residue = element.attrib["residue"]
        elif (event == "start") and (tag == "hmmEnd"):
            hmm_end = element.text
        elif (event == "start") and (tag == "hmmStart"):
            hmm_start = element.text

    yield HMMER3Site(desc, start, end, hmm_start, hmm_end, residue)
    return


def process_hmmer3_with_sites(  # noqa: C901
    element: ET.Element,
    itree: Iterator[Tuple[str, ET.Element]],
    ips_version: Optional[str],
    query_id: Optional[str],
    namespace: str = NAMESPACE
) -> Iterator[GFFRecord]:
    assert ips_version is not None
    assert query_id is not None

    library = None
    library_version = None

    seqid = query_id
    score = element.attrib["score"]
    evalue = element.attrib["evalue"]

    sites: List[HMMER3Site] = []

    while True:
        event, element = next(itree)
        tag = get_tag(element, namespace)

        if (event == "end") and (tag == "hmmer3-match-with-sites"):
            break

        if (event == "start") and (tag == "signature"):
            signature = process_signature(element, itree)
            accession = signature.accession
            name = signature.name
            desc = signature.desc
            library = signature.library
            library_version = signature.library_version
            dbxrefs = signature.dbxrefs
            goterms = signature.goterms
            kind = signature.kind
        elif (event == "start") and (tag == "hmmer3-location-with-sites"):
            start = element.attrib["start"]
            end = element.attrib["end"]
            env_start = element.attrib["env-start"]
            env_end = element.attrib["env-end"]
            hmm_start = element.attrib["hmm-start"]
            hmm_end = element.attrib["hmm-end"]
            complete = element.attrib["hmm-bounds"] == "COMPLETE"

        elif (event == "start") and (tag == "hmmer3-site"):
            sites.extend(process_hmmer3_sites(element, itree, namespace))

    source = get_source(ips_version, library, library_version)
    type_ = "protein_hmm_match"

    custom = {
        "evalue": str(evalue),
        "score": str(score),
        "accession": str(accession)
    }

    if complete:
        custom["complete"] = "true"

    if (desc is not None) and (desc.strip() != ""):
        custom["description"] = desc.strip()

    if (kind is not None) and (kind.strip() != ""):
        custom["kind"] = kind.strip()

    if (name is not None) and (name.strip() == ""):
        name = None

    if int(hmm_end) > 1:
        custom.update({
            "hmm_start": str(hmm_start),
            "hmm_end": str(hmm_end),
            "env_start": str(env_start),
            "env_end": str(env_end),
        })

        target = Target(accession, int(hmm_start) - 1, int(hmm_end))
        assert int(hmm_start) <= int(hmm_end), element
    else:
        target = None

    assert evalue is not None

    rec = GFFRecord(
        seqid=seqid,
        source=source,
        type=type_,
        start=int(start) - 1,
        end=int(end),
        score=float(evalue),
        phase=Phase.NOT_CDS,
        strand=Strand.PLUS,
        attributes=GFFAttributes(
            name=name,
            target=target,
            dbxref=dbxrefs,
            ontology_term=goterms,
            custom=custom
        )
    )
    assert int(start) <= int(end), str(rec)
    yield rec

    for site in sites:
        type_ = "biochemical_region_of_peptide"  # SO:0100001
        custom = {"kind": "ACT_SITE"}

        if site.residue is not None:
            custom["residue"] = site.residue

        if site.hmm_start is not None:
            custom["hmm_start"] = site.hmm_start

        if site.hmm_end is not None:
            custom["hmm_end"] = site.hmm_end

        if name is not None:
            note = [f"Part of {name}"]

        if site.desc is not None:
            custom["description"] = site.desc

        assert site.start is not None
        assert site.end is not None

        rec = GFFRecord(
            seqid=seqid,
            source=source,
            type=type_,
            start=int(site.start) - 1,
            end=int(site.end),
            phase=Phase.NOT_CDS,
            strand=Strand.PLUS,
            attributes=GFFAttributes(
                note=note,
                custom=custom,
            )
        )
        assert int(site.start) <= int(site.end), str(rec)
        yield rec
    return


def process_patternscan(
    element: ET.Element,
    itree: Iterator[Tuple[str, ET.Element]],
    ips_version: Optional[str],
    query_id: Optional[str],
    namespace: str = NAMESPACE
) -> Iterator[GFFRecord]:
    from copy import copy
    assert ips_version is not None
    assert query_id is not None

    library = None
    library_version = None

    seqid = query_id
    locations = []

    while True:
        event, element = next(itree)
        tag = get_tag(element, namespace)

        if (event == "end") and (tag == "patternscan-match"):
            break

        if (event == "start") and (tag == "signature"):
            signature = process_signature(element, itree)
            accession = signature.accession
            name = signature.name
            desc = signature.desc
            library = signature.library
            library_version = signature.library_version
            dbxrefs = signature.dbxrefs
            goterms = signature.goterms
            kind = signature.kind
        elif (event == "start") and (tag == "patternscan-location"):
            start = element.attrib["start"]
            end = element.attrib["end"]
            level = element.attrib["level"]

            locations.append((start, end, level))

    source = get_source(ips_version, library, library_version)
    type_ = "protein_match"

    custom = {
        "accession": str(accession)
    }

    if (desc is not None) and (desc.strip() != ""):
        custom["description"] = desc.strip()

    if (kind is not None) and (kind.strip() != ""):
        custom["kind"] = kind.strip()

    if (name is not None) and (name.strip() == ""):
        name = None

    target = None

    for start, end, level in locations:
        this_custom = copy(custom)
        this_custom["level"] = level

        rec = GFFRecord(
            seqid=seqid,
            source=source,
            type=type_,
            start=int(start) - 1,
            end=int(end),
            score=None,
            phase=Phase.NOT_CDS,
            strand=Strand.PLUS,
            attributes=GFFAttributes(
                name=name,
                target=target,
                dbxref=dbxrefs,
                ontology_term=goterms,
                custom=this_custom
            )
        )
        yield rec
    return


def printit(
    outfile: TextIO
) -> Callable[[Iterable[GFFRecord]], None]:
    def inner(records: Iterable[GFFRecord]):
        for record in records:
            print(record.as_str(), file=outfile)
    return inner


def inner(  # noqa: C901
    infile: TextIO,
    outfile: TextIO,
    namespace: str,
):
    itree = ET.iterparse(
        infile,
        events=["start", "end"]
    )
    ips_version = None
    query_id = None
    printer = printit(outfile)

    while True:
        try:
            event, element = next(itree)
        except StopIteration:
            break

        tag = get_tag(element, namespace)

        if (event == "start") and (tag == "protein-matches"):
            ips_version = element.attrib["interproscan-version"]
        elif (tag == "protein"):
            query_id = None
        elif (event == "start") and (tag == "xref"):
            query_id = element.attrib["id"]
        elif (event == "start") and (tag == "coils-match"):
            coils = process_coils(element, itree, ips_version, query_id)
            printer(coils)
        elif (event == "start") and (tag == "mobidblite-match"):
            mdb = process_mobidblite(element, itree, ips_version, query_id)
            printer(mdb)
        elif (event == "start") and (tag == "hmmer3-match"):
            hmm = process_hmmer3(element, itree, ips_version, query_id)
            printer(hmm)
        elif (event == "start") and (tag == "panther-match"):
            hmm = process_panther(element, itree, ips_version, query_id)
            printer(hmm)
        elif (event == "start") and (tag == "hmmer2-match"):
            hmm = process_hmmer2(element, itree, ips_version, query_id)
            printer(hmm)
        elif (event == "start") and (tag == "profilescan-match"):
            hmm = process_profilescan(element, itree, ips_version, query_id)
            printer(hmm)
        elif (event == "start") and (tag == "rpsblast-match"):
            hmm = process_rpsblast(element, itree, ips_version, query_id)
            printer(hmm)
        elif (event == "start") and (tag == "superfamilyhmmer3-match"):
            hmm = process_superfamilyhmmer3(
                element,
                itree,
                ips_version,
                query_id
            )
            printer(hmm)
        elif (event == "start") and (tag == "fingerprints-match"):
            hmm = process_fingerprints(element, itree, ips_version, query_id)
            printer(hmm)
        elif (event == "start") and (tag == "hmmer3-match-with-sites"):
            hmm = process_hmmer3_with_sites(
                element,
                itree,
                ips_version,
                query_id
            )
            printer(hmm)
        elif (event == "start") and (tag == "patternscan-match"):
            hmm = process_patternscan(element, itree, ips_version, query_id)
            printer(hmm)
        elif tag not in {'matches', 'protein-matches', 'sequence', 'xref'}:
            raise ValueError(f"unexpected tag {tag} {element}")

    return


def runner(args: argparse.Namespace) -> None:
    try:
        inner(
            infile=args.xml,
            outfile=args.outfile,
            namespace=args.namespace,
        )
    except Exception as e:
        raise e
    return
