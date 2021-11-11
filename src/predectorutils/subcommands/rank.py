#!/usr/bin/env python3

import sys
import argparse
import json

from collections import defaultdict
from statistics import median

from typing import (Set, Dict, List, Sequence, Union)

import numpy as np
import pandas as pd
import xgboost as xgb

from predectorutils.data import (
    get_interesting_dbcan_ids,
    get_interesting_pfam_ids,
    get_ltr_model,
)
from predectorutils.gff import GFFRecord
from predectorutils.analyses import (
    Analyses,
    Analysis,
    ApoplastP,
    DeepSig,
    EffectorP1,
    EffectorP2,
    EffectorP3,
    Phobius,
    SignalP3HMM,
    SignalP3NN,
    SignalP4,
    SignalP5,
    SignalP6,
    TargetPNonPlant,
    TMHMM,
    LOCALIZER,
    DeepLoc,
    DBCAN,
    PfamScan,
    PepStats,
    PHIBase,
    EffectorDB,
    DeepredeffFungi,
    Kex2SiteAnalysis,
    RXLRLikeAnalysis
)

#  RegexAnalysis,

COLUMNS = [
    "name",
    "effector_score",
    "manual_effector_score",
    "manual_secretion_score",
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
    "effectorp1",
    "effectorp2",
    "effectorp3_cytoplasmic",
    "effectorp3_apoplastic",
    "effectorp3_noneffector",
    "deepredeff_fungi",
    "is_secreted",
    "any_signal_peptide",
    "apoplastp",
    "single_transmembrane",
    "multiple_transmembrane",
    "molecular_weight",
    "residue_number",
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
    "fykin_gap",
    "kex2_cutsites",
    "rxlr_like_motifs",
    "localizer_nuclear",
    "localizer_chloro",
    "localizer_mito",
    "signalp3_nn",
    "signalp3_hmm",
    "signalp4",
    "signalp5",
    "signalp6",
    "deepsig",
    "phobius_sp",
    "phobius_tmcount",
    "tmhmm_tmcount",
    "tmhmm_first_60",
    "tmhmm_exp_aa",
    "tmhmm_first_tm_sp_coverage",
    "targetp_secreted",
    "targetp_secreted_prob",
    "targetp_mitochondrial_prob",
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
    "signalp3_nn_d",
    "signalp3_hmm_s",
    "signalp4_d",
    "signalp5_prob",
    "signalp6_prob",
]


def cli(parser: argparse.ArgumentParser) -> None:

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

    parser.add_argument(
        "--dbcan",
        type=argparse.FileType('r'),
        default=None,
        help="The dbcan matches to parse as input. Use '-' for stdin."
    )

    parser.add_argument(
        "--pfam",
        type=argparse.FileType('r'),
        default=None,
        help="The pfam domains to parse as input. Use '-' for stdin."
    )

    parser.add_argument(
        "--secreted-weight",
        type=float,
        default=3,
        help=(
            "The weight to give a protein if it is predicted to be secreted "
            "by any signal peptide method."
        )
    )

    parser.add_argument(
        "--sigpep-good-weight",
        type=float,
        default=0.5,
        help=(
            "The weight to give a protein if it is predicted to have a signal "
            "peptide by one of the more reliable methods."
        )
    )

    parser.add_argument(
        "--sigpep-ok-weight",
        type=float,
        default=0.25,
        help=(
            "The weight to give a protein if it is predicted to have a signal "
            "peptide by one of the reasonably reliable methods."
        )
    )

    parser.add_argument(
        "--single-transmembrane-weight",
        type=float,
        default=-2,
        help=(
            "The weight to give a protein if it is predicted to have "
            "> 0 TM domains by either method and not both > 1 (mutually "
            "exclusive with multiple-transmembrane-score). "
            "This is not applied if TMHMM first 60 is > 10. "
            "Use negative numbers to penalise."
        )
    )

    parser.add_argument(
        "--multiple-transmembrane-weight",
        type=float,
        default=-6,
        help=(
            "The weight to give a protein if it is predicted to have"
            "transmembrane have > 1 TM domains by both TMHMM and Phobius."
            "Use negative numbers to penalise."
        )
    )

    parser.add_argument(
        "--deeploc-extracellular-weight",
        type=float,
        default=1,
        help=(
            "The weight to give a protein if it is predicted to be "
            "extracellular by deeploc."
        )
    )

    parser.add_argument(
        "--deeploc-intracellular-weight",
        type=float,
        default=-0.5,
        help=(
            "The weigt to give a protein if it is predicted to be "
            "intracellular by deeploc. Use negative numbers to penalise."
        )
    )

    parser.add_argument(
        "--deeploc-membrane-weight",
        type=float,
        default=-1,
        help=(
            "The weight to give a protein if it is predicted to be "
            "membrane associated by deeploc. Use negative numbers to penalise."
        )
    )

    parser.add_argument(
        "--targetp-mitochondrial-weight",
        type=float,
        default=-0.5,
        help=(
            "The weight to give a protein if it is predicted to be "
            "mitochondrial by targetp. Use negative numbers to penalise."
        )
    )

    parser.add_argument(
        "--effectorp1-weight",
        type=float,
        default=3,
        help=(
            "The weight to give a protein if it is predicted to be "
            "an effector by effectorp1."
        )
    )

    parser.add_argument(
        "--effectorp2-weight",
        type=float,
        default=3,
        help=(
            "The weight to give a protein if it is predicted to be "
            "an effector by effectorp2."
        )
    )

    parser.add_argument(
        "--effector-homology-weight",
        type=float,
        default=5,
        help=(
            "The weight to give a protein if it is similar to a known "
            "effector or effector domain."
        )
    )

    parser.add_argument(
        "--virulence-homology-weight",
        type=float,
        default=1,
        help=(
            "The weight to give a protein if it is similar to a known "
            "protein that may be involved in virulence."
        )
    )

    parser.add_argument(
        "--lethal-homology-weight",
        type=float,
        default=-5,
        help=(
            "The weight to give a protein if it is similar to a known "
            "protein in phibase which caused a lethal phenotype."
        )
    )

    parser.add_argument(
        "--tmhmm-first-60-threshold",
        type=float,
        default=10,
        help=(
            "The minimum number of AAs predicted to be transmembrane in the "
            "first 60 AAs to consider a protein with a single TM domain "
            "a false positive (caused by hydrophobic region in sp)."
        )
    )

    return


def effector_score_it(
    record: Dict[str, Union[None, int, float, str]],
    effectorp1: float = 3,
    effectorp2: float = 3,
    effector: float = 5,
    virulence: float = 2,
    lethal: float = -5,
) -> float:
    """ """
    score: float = 0

    secretion_score = record["manual_secretion_score"]
    assert isinstance(secretion_score, float)
    score += secretion_score

    effectorp1_prob = record["effectorp1"]
    effectorp2_prob = record["effectorp2"]

    assert isinstance(effectorp1_prob, float)
    assert isinstance(effectorp2_prob, float)
    score += 2 * (effectorp1_prob - 0.5) * effectorp1
    score += 2 * (effectorp2_prob - 0.5) * effectorp2

    assert isinstance(record["phibase_effector_match"], int)
    assert isinstance(record["effector_match"], int)
    assert isinstance(record["dbcan_match"], int)
    assert isinstance(record["pfam_match"], int)
    has_effector_match = sum([
        record["phibase_effector_match"],
        record["effector_match"],
        record["dbcan_match"],
        record["pfam_match"]
    ]) > 0

    score += int(has_effector_match) * effector

    assert isinstance(record["phibase_virulence_match"], int)
    if not has_effector_match:
        score += record["phibase_virulence_match"] * virulence

    assert isinstance(record["phibase_lethal_match"], int)
    score += record["phibase_lethal_match"] * lethal
    return score


def secretion_score_it(
    record: Dict[str, Union[None, int, float, str]],
    secreted: float = 3,
    less_trustworthy_signal_prediction: float = 0.25,
    trustworthy_signal_prediction: float = 0.5,
    single_transmembrane: float = -1,
    multiple_transmembrane: float = -10,
    deeploc_extracellular: float = 1,
    deeploc_intracellular: float = -2,
    deeploc_membrane: float = -2,
    targetp_mitochondrial: float = -2,
) -> float:
    score: float = 0

    assert isinstance(record["is_secreted"], int)
    score += int(record["is_secreted"]) * secreted

    for k in ["signalp3_hmm", "signalp3_nn", "phobius_sp", "deepsig"]:
        v = record.get(k, 0)
        assert isinstance(v, int)
        score += v * less_trustworthy_signal_prediction

    for k in ["signalp4", "signalp5", "targetp_secreted"]:
        v = record.get(k, 0)
        assert isinstance(v, int)
        score += v * trustworthy_signal_prediction

    assert isinstance(record["multiple_transmembrane"], int)
    score += record["multiple_transmembrane"] * multiple_transmembrane

    assert isinstance(record["single_transmembrane"], int)
    score += record["single_transmembrane"] * single_transmembrane

    assert isinstance(record["deeploc_extracellular"], float)
    score += record["deeploc_extracellular"] * deeploc_extracellular  # noqa
    for k in [
        'deeploc_nucleus',
        'deeploc_cytoplasm',
        'deeploc_mitochondrion',
        'deeploc_cell_membrane',
        'deeploc_endoplasmic_reticulum',
        'deeploc_plastid',
        'deeploc_golgi',
        'deeploc_lysosome',
        'deeploc_peroxisome'
    ]:
        v = record.get(k, 0.0)
        assert isinstance(v, float)
        score += v * deeploc_intracellular

    assert isinstance(record["deeploc_membrane"], float)
    score += record["deeploc_membrane"] * deeploc_membrane

    assert isinstance(record["targetp_mitochondrial_prob"], float)
    score += record["targetp_mitochondrial_prob"] * targetp_mitochondrial
    return score


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
    record: Dict[str, Union[None, int, float, str]]
) -> None:
    record["any_signal_peptide"] = int(any([
        record.get(k, None) == 1
        for k
        in [
            'signalp3_nn', 'signalp3_hmm', 'signalp4',
            'signalp5', 'deepsig', 'phobius_sp', 'targetp_secreted'
        ]
    ]))
    return


def gff_intersection(left: GFFRecord, right: GFFRecord) -> int:
    lstart = min([left.start, left.end])
    lend = max([left.start, left.end])

    rstart = min([right.start, right.end])
    rend = max([right.start, right.end])

    start = max([lstart, rstart])
    end = min([lend, rend])

    # This will be < 0 if they don't overlap
    if start < end:
        return end - start
    else:
        return 0


def gff_coverage(left: GFFRecord, right: GFFRecord) -> float:
    noverlap = gff_intersection(left, right)
    return noverlap / (right.end - right.start)


def get_tm_sp_coverage(
    sp_gff: Sequence[GFFRecord],
    tm_gff: Sequence[GFFRecord],
) -> float:
    if len(sp_gff) == 0:
        return 0.0

    if len(tm_gff) == 0:
        return 0.0

    tm = sorted(tm_gff, key=lambda x: x.start)[0]
    covs = [gff_coverage(sp, tm) for sp in sp_gff]
    return median(covs)


def decide_is_transmembrane(
    record: Dict[str, Union[None, int, float, str]],
    sp_gff: Sequence[GFFRecord],
    tm_gff: Sequence[GFFRecord],
    tmhmm_first_60_threshold: float = 10,
) -> None:

    assert isinstance(record["tmhmm_tmcount"], int)
    assert isinstance(record["phobius_tmcount"], int)
    assert isinstance(record["tmhmm_first_60"], float)

    record["tmhmm_first_tm_sp_coverage"] = get_tm_sp_coverage(sp_gff, tm_gff)

    record["multiple_transmembrane"] = int(
        (record["tmhmm_tmcount"] > 1) or
        (record["phobius_tmcount"] > 1)
    )

    record["single_transmembrane"] = int(
        not bool(record["multiple_transmembrane"])
        and (
            record["phobius_tmcount"] == 1
            or (
                record["tmhmm_tmcount"] == 1
                and not bool(record["any_signal_peptide"])
            )
            or (
                record["tmhmm_tmcount"] == 1
                and bool(record["any_signal_peptide"])
                and (record["tmhmm_first_60"] < tmhmm_first_60_threshold)
            )
        )
    )

    return


def decide_is_secreted(
    record: Dict[str, Union[None, int, float, str]]
) -> None:
    record["is_secreted"] = int(
        bool(record["any_signal_peptide"]) and not
        bool(record["multiple_transmembrane"])
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
    an: PepStats,
    record: Dict[str, Union[None, int, float, str]]
) -> None:
    record["molecular_weight"] = an.molecular_weight
    record["residue_number"] = an.residues
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

    fykin = (
        an.residue_f_number +
        an.residue_k_number +
        an.residue_y_number +
        an.residue_i_number +
        an.residue_n_number
    )

    gap = (
        an.residue_g_number +
        an.residue_a_number +
        an.residue_p_number
    )

    record["fykin_gap"] = (float(fykin) + 1) / (float(gap) + 1)
    return


def get_phibase_cols(
    matches: Set[str],
    phenotypes: Set[str],
    record: Dict[str, Union[None, int, float, str]]
) -> None:
    record["phibase_effector_match"] = int(len(
        phenotypes.intersection([
            "loss_of_pathogenicity",
            "increased_virulence_(hypervirulence)",
            "effector_(plant_avirulence_determinant)"
        ])
    ) > 0)

    record["phibase_virulence_match"] = int("reduced_virulence" in phenotypes)
    record["phibase_lethal_match"] = int("lethal" in phenotypes)

    if len(phenotypes) > 0:
        record["phibase_phenotypes"] = ",".join(phenotypes)

    if len(matches) > 0:
        record["phibase_matches"] = ",".join(matches)
    return


def get_sper_prob_col(
    an: Union[EffectorP1, EffectorP2, ApoplastP],
    positive: List[str]
) -> float:
    if an.prediction in positive:
        return an.prob
    else:
        return 1 - an.prob


def construct_row(  # noqa
    name,
    analyses: List[Analysis],
    pfam_domains: Set[str],
    dbcan_domains: Set[str],
    tmhmm_first_60_threshold: float,
) -> Dict[str, Union[None, int, float, str]]:

    phibase_matches: Set[str] = set()
    phibase_phenotypes: Set[str] = set()
    effector_matches: Set[str] = set()
    pfam_matches: Set[str] = set()
    dbcan_matches: Set[str] = set()
    kex2_cutsites: Set[str] = set()
    rxlr_like_motifs: Set[str] = set()

    record: Dict[str, Union[None, int, float, str]] = {"name": name}
    record["effector_match"] = 0

    tm_gff: List[GFFRecord] = []
    sp_gff: List[GFFRecord] = []

    for an in analyses:
        if isinstance(an, ApoplastP):
            record["apoplastp"] = get_sper_prob_col(an, ["Apoplastic"])

        elif isinstance(an, DeepSig):
            record["deepsig"] = int(an.prediction == "SignalPeptide")
            sp_gff.extend(an.as_gff())

        elif isinstance(an, EffectorP1):
            record["effectorp1"] = get_sper_prob_col(an, ["Effector"])

        elif isinstance(an, EffectorP2):
            record["effectorp2"] = get_sper_prob_col(
                an,
                ["Effector", "Unlikely effector"]
            )

        elif isinstance(an, EffectorP3):
            record["effectorp3_cytoplasmic"] = an.cytoplasmic_prob
            record["effectorp3_apoplastic"] = an.apoplastic_prob
            record["effectorp3_noneffector"] = an.noneffector_prob

        elif isinstance(an, Phobius):
            record["phobius_sp"] = int(an.sp)
            record["phobius_tmcount"] = an.tm
            sp_gff.extend(
                r for r in an.as_gff()
                if r.type == "signal_peptide"
            )

        elif isinstance(an, SignalP3HMM):
            record["signalp3_hmm"] = int(an.is_secreted)
            record["signalp3_hmm_s"] = an.sprob
            sp_gff.extend(an.as_gff())

        elif isinstance(an, SignalP3NN):
            record["signalp3_nn"] = int(an.d_decision)
            record["signalp3_nn_d"] = an.d
            sp_gff.extend(an.as_gff())

        elif isinstance(an, SignalP4):
            record["signalp4"] = int(an.decision)
            record["signalp4_d"] = an.d
            record["signalp4_dmax_cut"] = an.dmax_cut
            sp_gff.extend(an.as_gff())

        elif isinstance(an, SignalP5):
            record["signalp5"] = int(an.prediction == "SP(Sec/SPI)")
            # For some proteins, this outputs a very high number
            # so we constrain it here.
            record["signalp5_prob"] = min([an.prob_signal, 1.0])
            sp_gff.extend(an.as_gff())

        elif isinstance(an, SignalP6):
            record["signalp6"] = int(an.prediction == "SP")
            record["signalp6_prob"] = min([an.prob_signal, 1.0])
            sp_gff.extend(an.as_gff())

        elif isinstance(an, TargetPNonPlant):
            record["targetp_secreted"] = int(an.prediction == "SP")
            record["targetp_secreted_prob"] = an.sp
            record["targetp_mitochondrial_prob"] = an.mtp

        elif isinstance(an, TMHMM):
            record["tmhmm_tmcount"] = an.pred_hel
            record["tmhmm_first_60"] = an.first_60
            record["tmhmm_exp_aa"] = an.exp_aa
            tm_gff.extend(an.as_gff())

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

        elif isinstance(an, EffectorDB):
            if an.decide_significant():
                record["effector_match"] = 1
                effector_matches.add(an.hmm)

        elif isinstance(an, DeepredeffFungi):
            record["deepredeff_fungi"] = float(an.s_score)

        elif isinstance(an, Kex2SiteAnalysis):
            kex2_cutsites.add(f"{an.pattern}:{an.start + 1}-{an.end}")

        elif isinstance(an, RXLRLikeAnalysis):
            rxlr_like_motifs.add(f"{an.pattern}:{an.start + 1}-{an.end}")

    decide_any_signal(record)
    decide_is_transmembrane(record, sp_gff, tm_gff, tmhmm_first_60_threshold)
    decide_is_secreted(record)

    get_phibase_cols(phibase_matches, phibase_phenotypes, record)

    if len(kex2_cutsites) > 0:
        record["kex2_cutsites"] = ",".join(kex2_cutsites)

    if len(rxlr_like_motifs) > 0:
        record["rxlr_like_motifs"] = ",".join(rxlr_like_motifs)

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

    return record


def run_ltr(df: pd.DataFrame) -> np.ndarray:
    df = df.copy()

    df["aa_c_prop"] = df["aa_c_number"] / df["residue_number"]
    df["aa_tiny_prop"] = df["aa_tiny_number"] / df["residue_number"]
    df["aa_small_prop"] = df["aa_small_number"] / df["residue_number"]
    df["aa_aliphatic_prop"] = df["aa_aliphatic_number"] / df["residue_number"]
    df["aa_aromatic_prop"] = df["aa_aromatic_number"] / df["residue_number"]
    df["aa_nonpolar_prop"] = df["aa_nonpolar_number"] / df["residue_number"]
    df["aa_charged_prop"] = df["aa_charged_number"] / df["residue_number"]
    df["aa_basic_prop"] = df["aa_basic_number"] / df["residue_number"]
    df["aa_acidic_prop"] = df["aa_acidic_number"] / df["residue_number"]

    df_features = df[[
        'molecular_weight',
        'aa_c_prop',
        'aa_tiny_prop',
        'aa_small_prop',
        'aa_nonpolar_prop',
        'aa_basic_prop',
        'effectorp1',
        'effectorp2',
        'apoplastp',
        'phobius_tmcount',
        'tmhmm_tmcount',
        'tmhmm_first_60',
        'deeploc_membrane',
        'deeploc_extracellular',
        'deepsig',
        'phobius_sp',
        'signalp3_nn_d',
        'signalp3_hmm_s',
        'signalp4_d',
        'signalp5_prob',
        'targetp_secreted_prob',
    ]]

    dmat = xgb.DMatrix(df_features)
    model = get_ltr_model()
    return model.predict(dmat)


def record_to_series(record: Dict[str, Union[None, int, float, str]]) -> str:
    return pd.Series([record.get(c, None) for c in COLUMNS], index=COLUMNS)


def write_line(record: Dict[str, Union[None, int, float, str]]) -> str:
    line = "\t".join(str(record.get(c, ".")) for c in COLUMNS)
    return line


def runner(args: argparse.Namespace) -> None:
    records = defaultdict(list)

    if args.dbcan is not None:
        dbcan: Set[str] = {d.strip() for d in args.dbcan.readlines()}
    else:
        dbcan = set(get_interesting_dbcan_ids())

    if args.pfam is not None:
        pfam: Set[str] = {d.strip() for d in args.pfam.readlines()}
    else:
        pfam = set(get_interesting_pfam_ids())

    for line in args.infile:
        sline = line.strip()
        if sline == "":
            continue

        dline = json.loads(sline)
        record = get_analysis(dline)
        records[getattr(record, record.name_column)].append(record)

    out_records = []

    for name, protein_records in records.items():
        record = construct_row(
            name,
            protein_records,
            pfam,
            dbcan,
            args.tmhmm_first_60_threshold
        )

        record["manual_secretion_score"] = secretion_score_it(
            record,
            args.secreted_weight,
            args.sigpep_ok_weight,
            args.sigpep_good_weight,
            args.single_transmembrane_weight,
            args.multiple_transmembrane_weight,
            args.deeploc_extracellular_weight,
            args.deeploc_intracellular_weight,
            args.deeploc_membrane_weight,
            args.targetp_mitochondrial_weight,
        )
        record["manual_effector_score"] = effector_score_it(
            record,
            args.effectorp1_weight,
            args.effectorp2_weight,
            args.effector_homology_weight,
            args.virulence_homology_weight,
            args.lethal_homology_weight,
        )
        out_records.append(record)

    df = pd.DataFrame([record_to_series(r) for r in out_records])
    df["effector_score"] = run_ltr(df)
    df.sort_values("effector_score", ascending=False, inplace=True)

    df.to_csv(args.outfile, sep="\t", index=False, na_rep=".")
    return
