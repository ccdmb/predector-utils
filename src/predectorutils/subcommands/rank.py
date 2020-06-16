#!/usr/bin/env python3

import sys
import argparse
import json
from collections import defaultdict

from typing import (Set, Dict, List, Union, Tuple)

from predectorutils.analyses import (
    Analyses,
    Analysis,
    ApoplastP,
    DeepSig,
    EffectorP1,
    EffectorP2,
    Phobius,
    SignalP3HMM,
    SignalP3NN,
    SignalP4,
    SignalP5,
    TargetPNonPlant,
    TMHMM,
    LOCALIZER,
    DeepLoc,
    DBCAN,
    PfamScan,
    PepStats,
    PHIBase,
    EffectorSearch
)

COLUMNS = [
    "name",
    "phibase_effector",
    "phibase_virulence",
    "phibase_lethal",
    "phibase_phenotypes",
    "phibase_matches",
    "effector_match",
    "effector_matches",
    "pfam_match",
    "pfam_matches",
    "dbcan_match",
    "dbcan_matches",
    "molecular_weight",
    "residues",
    "charge",
    "isoelectric_point",
    "aa_c_number",
    "aa_tiny_number",
    "aa_small_number",
    "aa_aliphatic_number",
    "aa_aromatic_number",
    "aa_nonpolar_number",
    "aa_charged_number",
    "aa_basic_number",
    "aa_acidic_number",
    "effectorp1",
    "effectorp2",
    "apoplastp",
    "localizer_nuclear",
    "localizer_chloro",
    "localizer_mito",
    "is_secreted",
    "any_signal_peptide",
    "signalp3_nn",
    "signalp3_hmm",
    "signalp4",
    "signalp5",
    "deepsig",
    "phobius_sp",
    "phobius_tmcount",
    "tmhmm_tmcount",
    "tmhmm_first60",
    "targetp_secreted",
    "targetp_mitochondrial",
    "deeploc_membrane",
    "deeploc_nucleus",
    "deeploc_cytoplasm",
    "deeploc_extracellular",
    "deeploc_mitochondrion",
    "deeploc_cell_membrane",
    "deeploc_endoplasmic_reticulum",
    "deeploc_plastid",
    "deeploc_golgi",
    "deeploc_lysosome",
    "deeploc_peroxisome",
]


def cli(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "dbcan",
        type=argparse.FileType('r'),
        help="The dbcan matches to parse as input. Use '-' for stdin."
    )

    parser.add_argument(
        "pfam",
        type=argparse.FileType('r'),
        help="The pfam domains to parse as input. Use '-' for stdin."
    )

    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help="The ldjson file to parse as input. Use '-' for stdin."
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Where to write the output to. Default: stdout"
    )

    return


def rank_it(
    record: Dict[str, Union[None, int, float, str]],
    secreted: float = 3,
    less_trustworthy_signal_prediction: float = 0.25,
    trustworthy_signal_prediction: float = 0.5,
    transmembrane: float = -10,
    deeploc_extracellular: float = 1,
    deeploc_intracellular: float = -2,
    targetp_secreted: float = 1,
    targetp_mitochondrial: float = -2,
    effectorp1: float = 3,
    effectorp2: float = 3,
    effector: float = 5,
    virulence: float = 2,
    lethal: float = -5,
):
    """ """

    score = 0

    score += int(record.get("is_secreted", 0)) * secreted
    for k in ["signalp3_hmm", "signalp3_nn", "phobius_sp"]:
        score += record.get(k, 0) * less_trustworthy_signal_prediction

    for k in ["signalp4", "signalp5", "deepsig"]:
        score += record.get(k, 0) * trustworthy_signal_prediction

    score += int(record.get("is_transmembrane", 0)) * transmembrane

    score += float(record.get("deeploc_extracellular", 0)) * deeploc_extracellular  # noqa
    for k in [
        'deeploc_membrane', 'deeploc_nucleus', 'deeploc_cytoplasm',
        'deeploc_mitochondrion', 'deeploc_cell_membrane', 'deeploc_plastid',
        'deeploc_lysosome', 'deeploc_peroxisome'
    ]:
        score += float(record.get(k, 0)) * deeploc_intracellular

    score += record["targetp_secreted"] * targetp_secreted
    score += record["targetp_mitochondrial"] * targetp_mitochondrial

    score += 2 * (float(record["effectorp1"]) - 0.5) * effectorp1
    score += 2 * (float(record["effectorp2"]) - 0.5) * effectorp2

    return


def get_analysis(dline):
    cls = Analyses.from_string(dline["analysis"]).get_analysis()
    analysis = cls.from_dict(dline["data"])
    return analysis


def parse_phibase_header(i):
    si = i.split("#")
    assert len(si) == 6
    return set(si[5].lower().split("__"))


def get_phibase_phis(i):
    si = i.split("#")
    assert len(si) == 6
    return set(si[1].split("__"))


def decide_any_signal(
    an: DeepLoc,
    record: Dict[str, Union[None, int, float, str]]
) -> None:
    record["any_signal_peptide"] = int(any([
        record.get(k, None) == 1
        for k
        in [
            'signalp3_nn', 'signalp3_hmm', 'signalp4',
            'signalp5', 'deepsig', 'phobius_sp'
        ]
    ]))
    return


def decide_is_transmembrane(
    an: DeepLoc,
    record: Dict[str, Union[None, int, float, str]]
) -> None:

    record["is_transmembrane"] = int(
        (record["tmhmm_tmcount"] > 1) or
        (record["phobius_tmcount"] > 0) or
        ((record["tmhmm_first60"] > 10) and
         (record["tmhmm_tmcount"] == 1) and
         bool(record["any_signal_peptide"]))
    )
    return


def decide_is_secreted(
    an: DeepLoc,
    record: Dict[str, Union[None, int, float, str]]
) -> None:
    record["is_secreted"] = int(
        bool(record["any_signal_peptide"]) and not
        bool(record["is_transmembrane"])
    )
    return


def get_deeploc_cols(
    an: DeepLoc,
    record: Dict[str, Union[None, int, float, str]]
) -> None:
    record["deeploc_membrane"] = an.membrane
    record["deeploc_nucleus"] = an.nucleus
    record["deeploc_cytoplasm"] = an.cytoplasm
    record["deeploc_extracellular"] = an.extracellular
    record["deeploc_mitochondrion"] = an.mitochondrion
    record["deeploc_cell_membrane"] = an.cell_membrane
    record["deeploc_endoplasmic_reticulum"] = an.endoplasmic_reticulum
    record["deeploc_plastid"] = an.plastid
    record["deeploc_golgi"] = an.golgi_apparatus
    record["deeploc_lysosome"] = an.lysosome_vacuole
    record["deeploc_peroxisome"] = an.peroxisome
    return


def get_pepstats_cols(
    an: DeepLoc,
    record: Dict[str, Union[None, int, float, str]]
) -> None:
    record["molecular_weight"] = an.molecular_weight
    record["residues"] = an.residues
    record["charge"] = an.charge
    record["isoelectric_point"] = an.isoelectric_point
    record["aa_c_number"] = an.residue_c_number
    record["aa_tiny_number"] = an.property_tiny_number
    record["aa_small_number"] = an.property_small_number
    record["aa_aliphatic_number"] = an.property_aliphatic_number
    record["aa_aromatic_number"] = an.property_aromatic_number
    record["aa_nonpolar_number"] = an.property_nonpolar_number
    record["aa_charged_number"] = an.property_charged_number
    record["aa_basic_number"] = an.property_basic_number
    record["aa_acidic_number"] = an.property_acidic_number
    return


def get_phibase_cols(
    matches: Set[str],
    phenotypes: Set[str],
    record: Dict[str, Union[None, int, float, str]]
) -> None:
    record["phibase_effector"] = int(len(
        phenotypes.intersection([
            "loss_of_pathogenicity",
            "increased_virulence_(hypervirulence)",
            "effector_(plant_avirulence_determinant)"
        ])
    ) > 0)

    record["phibase_virulence"] = int("reduced_virulence" in phenotypes)
    record["phibase_lethal"] = int("lethal" in phenotypes)

    if len(phenotypes) > 0:
        record["phibase_phenotypes"] = ",".join(phenotypes)

    if len(matches) > 0:
        record["phibase_matches"] = ",".join(matches)
    return


def get_sper_prob_col(
    an: DeepLoc,
    positive: Tuple[str]
) -> None:
    if an.prediction in positive:
        return an.prob
    else:
        return 1 - an.prob

def construct_row(  # noqa
    name,
    analyses: List[Analysis],
    pfam_domains: Set[str],
    dbcan_domains: Set[str]
):

    phibase_matches: Set[str] = set()
    phibase_phenotypes: Set[str] = set()
    effector_matches: Set[str] = set()
    pfam_matches: Set[str] = set()
    dbcan_matches: Set[str] = set()

    record: Dict[str, Union[None, int, float, str]] = {"name": name}
    record["effector_match"] = 0

    for an in analyses:
        if isinstance(an, ApoplastP):
            record["apoplastp"] = get_sper_prob_col(an, ("Apoplastic",))
        elif isinstance(an, DeepSig):
            record["deepsig"] = int(an.prediction == "SignalPeptide")
        elif isinstance(an, EffectorP1):
            record["effectorp1"] = get_sper_prob_col(an, ("Effector",))
        elif isinstance(an, EffectorP2):
            record["effectorp2"] = get_sper_prob_col(
                an,
                ("Effector", "Unlikely effector")
            )
        elif isinstance(an, Phobius):
            record["phobius_sp"] = int(an.sp)
            record["phobius_tmcount"] = an.tm
        elif isinstance(an, SignalP3HMM):
            record["signalp3_hmm"] = int(an.is_secreted)
            record["signalp3_hmm_s"] = an.sprob
        elif isinstance(an, SignalP3NN):
            record["signalp3_nn"] = int(an.d_decision)
            record["signalp3_nn_d"] = an.d
        elif isinstance(an, SignalP4):
            record["signalp4"] = int(an.decision)
            record["signalp4_d"] = an.d
            record["signalp4_dmax_cut"] = an.dmax_cut

        elif isinstance(an, SignalP5):
            record["signalp5"] = int(an.prediction == "SP(Sec/SPI)")
            record["signalp5_prob"] = an.prob_signal
        elif isinstance(an, TargetPNonPlant):
            record["targetp_secreted"] = an.sp
            record["targetp_mitochondrial"] = an.mtp
        elif isinstance(an, TMHMM):
            record["tmhmm_tmcount"] = an.pred_hel
            record["tmhmm_first60"] = an.first_60
        elif isinstance(an, LOCALIZER):
            record["localizer_nuclear"] = int(an.nucleus_decision)
            record["localizer_chloro"] = int(an.chloroplast_decision)
            record["localizer_mito"] = int(an.mitochondria_decision)
        elif isinstance(an, DeepLoc):
            get_deeploc_cols(an, record)
        elif isinstance(an, DBCAN):
            if an.decide_significant():
                dbcan_matches.add(an.hmm)
        elif isinstance(an, PfamScan):
            hmm = an.hmm.split('.', maxsplit=1)[0]
            if hmm in pfam_domains:
                pfam_matches.add(hmm)

        elif isinstance(an, PepStats):
            get_pepstats_cols(an, record)

        elif isinstance(an, PHIBase):
            if an.decide_significant():
                phibase_phenotypes.update(parse_phibase_header(an.target))
                phibase_matches.update(get_phibase_phis(an.target))

        elif isinstance(an, EffectorSearch):
            if an.decide_significant():
                record["effector_match"] = 1
                effector_matches.add(an.target)

    decide_any_signal(an, record)
    decide_is_transmembrane(an, record)
    decide_is_secreted(an, record)

    get_phibase_cols(phibase_matches, phibase_phenotypes, record)

    if len(effector_matches) > 0:
        record["effector_matches"] = ",".join(effector_matches)

    record["pfam_match"] = int(len(
        pfam_matches.intersection(pfam_domains)
    ) > 0)

    record["dbcan_match"] = int(len(
        dbcan_matches.intersection(dbcan_domains)
    ) > 0)

    if len(pfam_matches) > 0:
        record["pfam_matches"] = ",".join(pfam_matches)

    if len(dbcan_matches) > 0:
        record["dbcan_matches"] = ",".join(dbcan_matches)

    line = "\t".join(str(record.get(c, ".")) for c in COLUMNS)
    return line


def runner(args: argparse.Namespace) -> None:
    records = defaultdict(list)

    dbcan = {l.strip() for l in args.dbcan.readlines()}
    pfam = {l.strip() for l in args.pfam.readlines()}

    for line in args.infile:
        sline = line.strip()
        if sline == "":
            continue

        dline = json.loads(sline)
        record = get_analysis(dline)
        records[dline["protein_name"]].append(record)

    for name, protein_records in records.items():
        line = construct_row(name, protein_records, pfam, dbcan)
        print(line, file=args.outfile)
    return
