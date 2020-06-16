#!/usr/bin/env python3

import sys
import argparse
import json
from collections import defaultdict

from typing import (Set, Dict, List, Union)

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
    "score",
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
    "signalp3_nn_d",
    "signalp3_hmm_s",
    "signalp4_d",
    "signalp4_dmax_cut",
    "signalp5_prob",
]

DBCAN_DEFAULT = {
    "CBM50",
    "AA9",
    "GH10",
    "GH11"
}


PFAM_DEFAULT = {
    "PF00734", "PF03443", "PF01083", "PF03211", "PF00331",
    "PF00457", "PF00160", "PF04126", "PF10231", "PF18050",
    "PF00188", "PF00293", "PF00964", "PF03377", "PF05730",
    "PF07249", "PF16541", "PF03330", "PF11529", "PF18661",
    "PF18247", "PF18241", "PF14521", "PF11327", "PF01764",
    "PF14856", "PF08995", "PF05630", "PF09461", "PF11584",
    "PF18224", "PF00187", "PF01476", "PF01607", "PF02128",
    "PF07504", "PF02221", "PF01363", "PF06602", "PF16810",
    "PF18488", "PF18634", "PF00161", "PF00545", "PF00652",
    "PF11663", "PF14200", "PF16536", "PF16542", "PF00014",
    "PF00030", "PF00059", "PF00068", "PF00087", "PF00200",
    "PF00451", "PF00537", "PF00555", "PF00706", "PF00819",
    "PF01117", "PF01123", "PF01324", "PF01335", "PF01338",
    "PF01375", "PF01376", "PF01421", "PF01742", "PF01821",
    "PF02048", "PF02258", "PF02452", "PF02691", "PF02763",
    "PF02764", "PF02794", "PF02819", "PF02876", "PF02917",
    "PF02918", "PF02950", "PF03077", "PF03318", "PF03440",
    "PF03495", "PF03497", "PF03498", "PF03505", "PF03507",
    "PF03768", "PF03769", "PF03784", "PF03944", "PF03945",
    "PF04221", "PF04365", "PF04829", "PF05015", "PF05016",
    "PF05294", "PF05353", "PF05374", "PF05394", "PF05453",
    "PF05627", "PF05638", "PF05707", "PF05791", "PF05819",
    "PF05932", "PF06254", "PF06255", "PF06261", "PF06286",
    "PF06296", "PF06340", "PF06357", "PF06414", "PF06416",
    "PF06451", "PF06453", "PF06457", "PF06755", "PF06769",
    "PF07254", "PF07365", "PF07442", "PF07473", "PF07591",
    "PF07737", "PF07740", "PF07771", "PF07822", "PF07829",
    "PF07906", "PF07927", "PF07936", "PF07945", "PF07951",
    "PF07952", "PF07953", "PF07968", "PF08024", "PF08025",
    "PF08086", "PF08087", "PF08088", "PF08089", "PF08090",
    "PF08091", "PF08092", "PF08093", "PF08094", "PF08095",
    "PF08096", "PF08097", "PF08098", "PF08099", "PF08104",
    "PF08115", "PF08116", "PF08119", "PF08120", "PF08121",
    "PF08131", "PF08178", "PF08249", "PF08251", "PF08396",
    "PF08493", "PF08562", "PF08843", "PF08845", "PF08888",
    "PF09009", "PF09075", "PF09101", "PF09102", "PF09119",
    "PF09131", "PF09156", "PF09207", "PF09275", "PF09276",
    "PF09407", "PF09857", "PF09907", "PF10279", "PF10530",
    "PF10550", "PF11047", "PF11410", "PF11592", "PF12207",
    "PF12255", "PF12256", "PF12563", "PF12703", "PF12918",
    "PF13304", "PF13939", "PF13940", "PF13955", "PF13956",
    "PF13957", "PF13958", "PF13979", "PF13981", "PF14021",
    "PF14113", "PF14171", "PF14449", "PF14866", "PF15500",
    "PF15520", "PF15521", "PF15522", "PF15523", "PF15524",
    "PF15526", "PF15527", "PF15528", "PF15529", "PF15530",
    "PF15531", "PF15532", "PF15533", "PF15534", "PF15535",
    "PF15536", "PF15537", "PF15538", "PF15540", "PF15541",
    "PF15542", "PF15543", "PF15544", "PF15545", "PF15604",
    "PF15605", "PF15606", "PF15607", "PF15635", "PF15636",
    "PF15637", "PF15638", "PF15639", "PF15640", "PF15641",
    "PF15643", "PF15644", "PF15645", "PF15646", "PF15647",
    "PF15648", "PF15649", "PF15650", "PF15651", "PF15652",
    "PF15653", "PF15657", "PF15658", "PF15659", "PF15723",
    "PF15738", "PF15781", "PF15935", "PF16754", "PF16847",
    "PF16873", "PF16981", "PF17056", "PF17454", "PF17475",
    "PF17476", "PF17486", "PF17491", "PF17492", "PF17499",
    "PF17521", "PF17556", "PF17557", "PF17563", "PF17997",
    "PF18022", "PF18078", "PF18276", "PF18449", "PF18648",
    "PF18807", "PF06766", "PF16073", "PF11807", "PF00314",
    "PF16361"
}


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
        "--secreted-score",
        type=float,
        default=3,
        help="The score to give a protein if it is predicted to be secreted."
    )

    parser.add_argument(
        "--sigpep-good-score",
        type=float,
        default=0.5,
        help=(
            "The score to give a protein if it is predicted to have a signal "
            "peptide by one of the more reliable methods."
        )
    )

    parser.add_argument(
        "--sigpep-ok-score",
        type=float,
        default=0.25,
        help=(
            "The score to give a protein if it is predicted to have a signal "
            "peptide by one of the reasonably reliable methods."
        )
    )

    parser.add_argument(
        "--transmembrane-score",
        type=float,
        default=-10,
        help=(
            "The score to give a protein if it is predicted to be "
            "transmembrane. Use negative numbers to penalise."
        )
    )

    parser.add_argument(
        "--deeploc-extracellular-score",
        type=float,
        default=1,
        help=(
            "The score to give a protein if it is predicted to be "
            "extracellular by deeploc."
        )
    )

    parser.add_argument(
        "--deeploc-intracellular-score",
        type=float,
        default=-2,
        help=(
            "The score to give a protein if it is predicted to be "
            "intracellular by deeploc. Use negative numbers to penalise."
        )
    )

    parser.add_argument(
        "--targetp-secreted-score",
        type=float,
        default=1,
        help=(
            "The score to give a protein if it is predicted to be "
            "secreted by targetp."
        )
    )

    parser.add_argument(
        "--targetp-mitochondrial-score",
        type=float,
        default=-2,
        help=(
            "The score to give a protein if it is predicted to be "
            "mitochondrial by targetp. Use negative numbers to penalise."
        )
    )

    parser.add_argument(
        "--effectorp1-score",
        type=float,
        default=3,
        help=(
            "The score to give a protein if it is predicted to be "
            "an effector by effectorp1."
        )
    )

    parser.add_argument(
        "--effectorp2-score",
        type=float,
        default=3,
        help=(
            "The score to give a protein if it is predicted to be "
            "an effector by effectorp2."
        )
    )

    parser.add_argument(
        "--effector-homology-score",
        type=float,
        default=5,
        help=(
            "The score to give a protein if it is similar to a known "
            "effector or effector domain."
        )
    )

    parser.add_argument(
        "--virulence-homology-score",
        type=float,
        default=1,
        help=(
            "The score to give a protein if it is similar to a known "
            "protein that may be involved in virulence."
        )
    )

    parser.add_argument(
        "--lethal-homology-score",
        type=float,
        default=-5,
        help=(
            "The score to give a protein if it is similar to a known "
            "protein in phibase which caused a lethal phenotype."
        )
    )

    return


def score_it(
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
) -> float:
    """ """

    score: float = 0

    is_secreted = record.get("is_secreted", 0)
    assert isinstance(is_secreted, int)
    score += int(is_secreted) * secreted

    for k in ["signalp3_hmm", "signalp3_nn", "phobius_sp"]:
        v = record.get(k, 0)
        assert isinstance(v, int)
        score += v * less_trustworthy_signal_prediction

    for k in ["signalp4", "signalp5", "deepsig"]:
        v = record.get(k, 0)
        assert isinstance(v, int)
        score += v * trustworthy_signal_prediction

    is_transmembrane = record.get("is_transmembrane", 0)
    assert isinstance(is_transmembrane, int)
    score += is_transmembrane * transmembrane

    deeploc_extracellular_prob = record.get("deeploc_extracellular", 0.0)
    assert isinstance(deeploc_extracellular_prob, float)
    score += deeploc_extracellular_prob * deeploc_extracellular  # noqa
    for k in [
        'deeploc_membrane', 'deeploc_nucleus', 'deeploc_cytoplasm',
        'deeploc_mitochondrion', 'deeploc_cell_membrane', 'deeploc_plastid',
        'deeploc_lysosome', 'deeploc_peroxisome'
    ]:
        v = record.get(k, 0.0)
        assert isinstance(v, float)
        score += v * deeploc_intracellular

    targetp_secreted_prob = record.get("targetp_secreted", 0.0)
    targetp_mitochondrial_prob = record.get("targetp_mitochondrial", 0.0)

    assert isinstance(targetp_secreted_prob, float)
    assert isinstance(targetp_mitochondrial_prob, float)
    score += targetp_secreted_prob * targetp_secreted
    score += targetp_mitochondrial_prob * targetp_mitochondrial

    effectorp1_prob = record.get("effectorp1", 0.0)
    effectorp2_prob = record.get("effectorp2", 0.0)

    assert isinstance(effectorp1_prob, float)
    assert isinstance(effectorp2_prob, float)
    score += 2 * (effectorp1_prob - 0.5) * effectorp1
    score += 2 * (effectorp2_prob - 0.5) * effectorp2

    assert isinstance(record["phibase_effector"], int)
    assert isinstance(record["effector_match"], int)
    assert isinstance(record["dbcan_match"], int)
    assert isinstance(record["pfam_match"], int)
    has_effector_match = sum([
        record["phibase_effector"],
        record["effector_match"],
        record["dbcan_match"],
        record["pfam_match"]
    ]) > 0

    score += int(has_effector_match) * effector

    assert isinstance(record["phibase_virulence"], int)
    if not has_effector_match:
        score += record["phibase_virulence"] * virulence

    assert isinstance(record["phibase_lethal"], int)
    score += record["phibase_lethal"] * lethal

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
            'signalp5', 'deepsig', 'phobius_sp'
        ]
    ]))
    return


def decide_is_transmembrane(
    record: Dict[str, Union[None, int, float, str]]
) -> None:

    assert isinstance(record["tmhmm_tmcount"], int)
    assert isinstance(record["phobius_tmcount"], int)
    assert isinstance(record["tmhmm_first60"], float)
    assert isinstance(record["any_signal_peptide"], int)

    record["is_transmembrane"] = int(
        (record["tmhmm_tmcount"] > 1) or
        (record["phobius_tmcount"] > 0) or
        ((record["tmhmm_first60"] > 10) and
         (record["tmhmm_tmcount"] == 1) and
         bool(record["any_signal_peptide"]))
    )
    return


def decide_is_secreted(
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
    an: PepStats,
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
    dbcan_domains: Set[str]
) -> Dict[str, Union[None, int, float, str]]:

    phibase_matches: Set[str] = set()
    phibase_phenotypes: Set[str] = set()
    effector_matches: Set[str] = set()
    pfam_matches: Set[str] = set()
    dbcan_matches: Set[str] = set()

    record: Dict[str, Union[None, int, float, str]] = {"name": name}
    record["effector_match"] = 0

    for an in analyses:
        if isinstance(an, ApoplastP):
            record["apoplastp"] = get_sper_prob_col(an, ["Apoplastic"])
        elif isinstance(an, DeepSig):
            record["deepsig"] = int(an.prediction == "SignalPeptide")
        elif isinstance(an, EffectorP1):
            record["effectorp1"] = get_sper_prob_col(an, ["Effector"])
        elif isinstance(an, EffectorP2):
            record["effectorp2"] = get_sper_prob_col(
                an,
                ["Effector", "Unlikely effector"]
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

    decide_any_signal(record)
    decide_is_transmembrane(record)
    decide_is_secreted(record)

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

    return record


def write_line(record: Dict[str, Union[None, int, float, str]]) -> str:
    line = "\t".join(str(record.get(c, ".")) for c in COLUMNS)
    return line


def runner(args: argparse.Namespace) -> None:
    records = defaultdict(list)

    if args.dbcan is not None:
        dbcan = {l.strip() for l in args.dbcan.readlines()}

    if args.pfam is not None:
        pfam = {l.strip() for l in args.pfam.readlines()}

    for line in args.infile:
        sline = line.strip()
        if sline == "":
            continue

        dline = json.loads(sline)
        record = get_analysis(dline)
        records[dline["protein_name"]].append(record)

    print("\t".join(COLUMNS), file=args.outfile)

    for name, protein_records in records.items():
        record = construct_row(name, protein_records, pfam, dbcan)
        record["score"] = score_it(
            record,
            args.secreted_score,
            args.sigpep_ok_score,
            args.sigpep_good_score,
            args.transmembrane_score,
            args.deeploc_extracellular_score,
            args.deeploc_intracellular_score,
            args.targetp_secreted_score,
            args.targetp_mitochondrial_score,
            args.effectorp1_score,
            args.effectorp2_score,
            args.effector_homology_score,
            args.virulence_homology_score,
            args.lethal_homology_score,
        )

        line = write_line(record)
        print(line, file=args.outfile)
    return
