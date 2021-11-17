#!/usr/bin/env python3

import sys
import argparse

from collections import defaultdict
from statistics import median

from typing import (Tuple, Dict, Set, List, Sequence)

import sqlite3

import numpy as np
import pandas as pd
import xgboost as xgb

from predectorutils.data import (
    get_interesting_dbcan_ids,
    get_interesting_pfam_ids,
    get_ltr_model,
)

from predectorutils.indexedresults import ResultsTable, ResultRow

from predectorutils.gff import GFFRecord
from predectorutils.analyses import (
    GFFAble,
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
    DomTbl,
    PfamScan,
    PepStats,
    MMSeqs,
    DeepredeffFungi,
    DeepredeffOomycete,
    RegexAnalysis,
)


COLUMNS = [
    "name",
    "effector_score",
    "manual_effector_score",
    "manual_secretion_score",
    "effector_matches",
    "phibase_genes",
    "phibase_phenotypes",
    "phibase_ids",
    "has_phibase_effector_match",
    "has_phibase_virulence_match",
    "has_phibase_lethal_match",
    "pfam_ids",
    "pfam_names",
    "has_pfam_virulence_match",
    "dbcan_matches",
    "has_dbcan_virulence_match",
    "effectorp1",
    "effectorp2",
    "effectorp3_cytoplasmic",
    "effectorp3_apoplastic",
    "effectorp3_noneffector",
    "deepredeff_fungi",
    "deepredeff_oomycete",
    "apoplastp",
    "is_secreted",
    "any_signal_peptide",
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
    "localizer_nucleus",
    "localizer_chloro",
    "localizer_mito",
    "signal_peptide_cutsites",
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
        "--db",
        type=str,
        default=":memory:",
        help="Where to store the sqlite database."
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
        "--effectorp3-apoplastic-weight",
        type=float,
        default=3,
        help=(
            "The weight to give a protein if it is predicted to be "
            "an apoplastic effector by effectorp3."
        )
    )

    parser.add_argument(
        "--effectorp3-cytoplasmic-weight",
        type=float,
        default=3,
        help=(
            "The weight to give a protein if it is predicted to be "
            "a cytoplasmic effector by effectorp3."
        )
    )

    parser.add_argument(
        "--effectorp3-noneffector-weight",
        type=float,
        default=-3,
        help=(
            "The weight to give a protein if it is predicted to be "
            "a non-effector by effectorp3."
        )
    )

    parser.add_argument(
        "--deepredeff-fungi-weight",
        type=float,
        default=2,
        help=(
            "The weight to give a protein if it is predicted to be "
            "a fungal effector by deepredeff."
        )
    )

    parser.add_argument(
        "--deepredeff-oomycete-weight",
        type=float,
        default=2,
        help=(
            "The weight to give a protein if it is predicted to be "
            "an oomycete effector by deepredeff."
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


def create_tables(
    tab: ResultsTable,
    con: sqlite3.Connection,
    cur: sqlite3.Cursor,
    pfam_targets: Set[str],
    dbcan_targets: Set[str],
    tmhmm_first_60_threshold: float = 10,
) -> pd.DataFrame:
    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS apoplastp
    AS
    SELECT
        json_extract(data, '$.{ApoplastP.name_column}') as name,
        checksum,
        json_extract(data, '$.prob') as apoplastp
    FROM results
    WHERE (analysis = 'apoplastp'
        AND json_extract(data, '$.prediction') = 'Apoplastic')
    UNION
    SELECT
        json_extract(data, '$.{ApoplastP.name_column}') as name,
        checksum,
        (1 - json_extract(data, '$.prob')) as apoplastp
    FROM results
    WHERE (analysis = 'apoplastp'
        AND NOT json_extract(data, '$.prediction') = 'Apoplastic')
    """)
    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS effectorp1
    AS
    SELECT
        json_extract(data, '$.{EffectorP1.name_column}') as name,
        checksum,
        json_extract(data, '$.prob') as effectorp1
    FROM results
    WHERE (analysis = 'effectorp1'
        AND json_extract(data, '$.prediction') = 'Effector')
    UNION
    SELECT
        json_extract(data, '$.{EffectorP1.name_column}') as name,
        checksum,
        (1 - json_extract(data, '$.prob')) as effectorp1
    FROM results
    WHERE (analysis = 'effectorp1'
        AND NOT json_extract(data, '$.prediction') = 'Effector')
    """)

    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS effectorp2
    AS
    SELECT
        json_extract(data, '$.{EffectorP2.name_column}') as name,
        checksum,
        json_extract(data, '$.prob') as effectorp2
    FROM results
    WHERE (
        analysis = 'effectorp2'
        AND (json_extract(data, '$.prediction')
                IN ('Effector', 'Unlikely effector'))
    )
    UNION
    SELECT
        json_extract(data, '$.{EffectorP2.name_column}') as name,
        checksum,
        (1 - json_extract(data, '$.prob')) as effectorp2
    FROM results
    WHERE (
        analysis = 'effectorp2'
        AND (json_extract(data, '$.prediction')
                NOT IN ('Effector', 'Unlikely effector'))
    )
    """)

    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS effectorp3
    AS
    SELECT
        json_extract(data, '$.{EffectorP3.name_column}') as name,
        checksum,
        json_extract(data, '$.cytoplasmic_prob') as effectorp3_cytoplasmic,
        json_extract(data, '$.apoplastic_prob') as effectorp3_apoplastic,
        json_extract(data, '$.noneffector_prob') as effectorp3_noneffector
    FROM results
    WHERE analysis = 'effectorp3'
    """)

    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS localizer
    AS
    SELECT
        json_extract(data, '$.{LOCALIZER.name_column}') as name,
        checksum,
        json_extract(data, '$.nucleus_decision') as localizer_nucleus,
        json_extract(data, '$.chloroplast_decision') as localizer_chloro,
        json_extract(data, '$.mitochondria_decision') as localizer_mito
    FROM results
    WHERE analysis = 'localizer'
    """)

    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS deepredeff_fungi
    AS
    SELECT
        json_extract(data, '$.{DeepredeffFungi.name_column}') as name,
        checksum,
        json_extract(data, '$.s_score') as deepredeff_fungi
    FROM results
    WHERE analysis = 'deepredeff_fungi'
    """)

    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS deepredeff_oomycete
    AS
    SELECT
        json_extract(data, '$.{DeepredeffOomycete.name_column}') as name,
        checksum,
        json_extract(data, '$.s_score') as deepredeff_oomycete
    FROM results
    WHERE analysis = 'deepredeff_oomycete'
    """)

    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS signalp3_nn
    AS
    SELECT
        json_extract(data, '$.{SignalP3NN.name_column}') as name,
        checksum,
        json_extract(data, '$.d_decision') as signalp3_nn,
        json_extract(data, '$.d') as signalp3_nn_d
    FROM results
    WHERE analysis = 'signalp3_nn'
    """)

    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS signalp3_hmm
    AS
    SELECT
        json_extract(data, '$.{SignalP3HMM.name_column}') as name,
        checksum,
        json_extract(data, '$.is_secreted') as signalp3_hmm,
        json_extract(data, '$.sprob') as signalp3_hmm_s
    FROM results
    WHERE analysis = 'signalp3_hmm'
    """)

    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS signalp4
    AS
    SELECT
        json_extract(data, '$.{SignalP4.name_column}') as name,
        checksum,
        json_extract(data, '$.decision') as signalp4,
        json_extract(data, '$.d') as signalp4_d,
        json_extract(data, '$.dmax_cut') as signalp4_dmax_cut
    FROM results
    WHERE analysis = 'signalp4'
    """)

    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS signalp5
    AS
    SELECT
        json_extract(data, '$.{SignalP5.name_column}') as name,
        checksum,
        json_extract(data, '$.prediction') == 'SP(Sec/SPI)' as signalp5,
        MIN(json_extract(data, '$.prob_signal'), 1.0) as signalp5_prob
    FROM results
    WHERE analysis = 'signalp5'
    """)

    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS signalp6
    AS
    SELECT
        json_extract(data, '$.{SignalP6.name_column}') as name,
        checksum,
        json_extract(data, '$.prediction') == 'SP' as signalp6,
        MIN(json_extract(data, '$.prob_signal'), 1.0) as signalp6_prob
    FROM results
    WHERE analysis = 'signalp6'
    """)

    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS deepsig
    AS
    SELECT
        json_extract(data, '$.{DeepSig.name_column}') as name,
        checksum,
        json_extract(data, '$.prediction') == 'SignalPeptide' as deepsig,
        json_extract(data, '$.prob') as deepsig_signal_prob,
        NULL as deepsig_transmembrane_prob,
        NULL as deepsig_other_prob
    FROM results
    WHERE
        analysis = 'deepsig'
        AND json_extract(data, '$.prediction') == 'SignalPeptide'
    UNION
    SELECT
        json_extract(data, '$.{DeepSig.name_column}') as name,
        checksum,
        json_extract(data, '$.prediction') == 'SignalPeptide' as deepsig,
        NULL as deepsig_signal_prob,
        json_extract(data, '$.prob') as deepsig_transmembrane_prob,
        NULL as deepsig_other_prob
    FROM results
    WHERE
        analysis = 'deepsig'
        AND json_extract(data, '$.prediction') == 'Transmembrane'
    UNION
    SELECT
        json_extract(data, '$.{DeepSig.name_column}') as name,
        checksum,
        json_extract(data, '$.prediction') == 'SignalPeptide' as deepsig,
        NULL as deepsig_signal_prob,
        NULL as deepsig_transmembrane_prob,
        json_extract(data, '$.prob') as deepsig_other_prob
    FROM results
    WHERE
        analysis = 'deepsig'
        AND json_extract(data, '$.prediction') == 'Other'
    """)

    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS phobius
    AS
    SELECT
        json_extract(data, '$.{Phobius.name_column}') as name,
        checksum,
        json_extract(data, '$.sp') as phobius_sp,
        json_extract(data, '$.tm') as phobius_tmcount
    FROM results
    WHERE analysis = 'phobius'
    """)

    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS tmhmm
    AS
    SELECT
        json_extract(data, '$.{TMHMM.name_column}') as name,
        checksum,
        json_extract(data, '$.pred_hel') as tmhmm_tmcount,
        json_extract(data, '$.first_60') as tmhmm_first_60,
        json_extract(data, '$.exp_aa') as tmhmm_exp_aa
    FROM results
    WHERE analysis = 'tmhmm'
    """)

    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS targetp
    AS
    SELECT
        json_extract(data, '$.{TargetPNonPlant.name_column}') as name,
        checksum,
        json_extract(data, '$.prediction') == 'SP' as targetp_secreted,
        json_extract(data, '$.sp') as targetp_secreted_prob,
        MIN(json_extract(data, '$.mtp'), 1.0) as targetp_mitochondrial_prob
    FROM results
    WHERE analysis = 'targetp_non_plant'
    """)

    cur.execute(f"""  -- # noqa
    CREATE TEMP TABLE IF NOT EXISTS deeploc
    AS
    SELECT
        json_extract(data, '$.{DeepLoc.name_column}') as name,
        checksum,
        json_extract(data, '$.membrane') as deeploc_membrane,
        json_extract(data, '$.nucleus') as deeploc_nucleus,
        json_extract(data, '$.cytoplasm') as deeploc_cytoplasm,
        json_extract(data, '$.extracellular') as deeploc_extracellular,
        json_extract(data, '$.mitochondrion') as deeploc_mitochondrion,
        json_extract(data, '$.cell_membrane') as deeploc_cell_membrane,
        json_extract(data, '$.endoplasmic_reticulum') as deeploc_endoplasmic_reticulum,
        json_extract(data, '$.plastid') as deeploc_plastid,
        json_extract(data, '$.golgi_apparatus') as deeploc_golgi,
        json_extract(data, '$.lysosome_vacuole') as deeploc_lysosome,
        json_extract(data, '$.peroxisome') as deeploc_peroxisome
    FROM results
    WHERE analysis = 'deeploc'
    """)

    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS hmmer_significant
    AS
    SELECT
        checksum,
        analysis,
        json_extract(data, '$.{DomTbl.name_column}') as name,
        ((json_extract(data, '$.hmm_to') - json_extract(data, '$.hmm_from')) /
                (1.0 * json_extract(data, '$.hmm_len'))) as coverage,
        (json_extract(data, '$.query_to') -
                json_extract(data, '$.query_from')) as length,
        json_extract(data, '$.domain_i_evalue') as evalue,
        json_extract(data, '$.hmm') as hmm
    FROM
        results
    WHERE analysis IN ("dbcan", "effectordb")
    """)

    cur.execute("""
    CREATE TEMP TABLE IF NOT EXISTS hmmer
    AS
    SELECT DISTINCT
        name, checksum, analysis, hmm
    FROM (
        SELECT *
        FROM hmmer_significant
        WHERE (coverage > 0.3) AND (length > 80) AND (evalue < 1e-5)
        UNION
        SELECT *
        FROM hmmer_significant
        WHERE (coverage > 0.3) AND (length <= 80) AND (evalue < 1e-3)
    )
    ORDER BY evalue ASC
    """)

    cur.execute("""
    CREATE TEMP TABLE IF NOT EXISTS effectordb
    AS
    SELECT
        name as name,
        checksum as checksum,
        analysis as analysis,
        GROUP_CONCAT(hmm, ',') as effector_matches
    FROM
        hmmer
    WHERE
        analysis = "effectordb"
    GROUP BY
        name, checksum, analysis
    """)

    cur.execute("DROP TABLE IF EXISTS dbcan_targets")
    cur.execute("CREATE TEMP TABLE dbcan_targets (dbcan_ids text NOT NULL)")
    cur.executemany(
        "INSERT INTO dbcan_targets VALUES (?)",
        [(c,) for c in dbcan_targets]
    )

    cur.execute("""
    CREATE TEMP TABLE IF NOT EXISTS dbcan
    AS
    SELECT
        name as name,
        checksum as checksum,
        analysis as analysis,
        MAX(dbcan_virulence_related) as has_dbcan_virulence_match,
        GROUP_CONCAT(hmm, ',') as dbcan_matches
    FROM (
        SELECT
            name,
            checksum,
            analysis,
            (hmm IN dbcan_targets) as dbcan_virulence_related,
            hmm
        FROM
            hmmer
        WHERE
            analysis = "dbcan"
    )
    GROUP BY
        name, checksum, analysis
    """)

    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS mmseqs
    AS
    SELECT DISTINCT
        json_extract(data, '$.{MMSeqs.name_column}') as name,
        checksum,
        analysis,
        json_extract(data, '$.target') as target
    FROM results
    WHERE
        analysis = 'phibase' AND
        json_extract(data, '$.evalue') <= 1e-5 AND
        json_extract(data, '$.tcov') >= 0.5
    ORDER BY json_extract(data, '$.evalue') ASC
    """)

    cur.execute("""
    CREATE TEMP TABLE IF NOT EXISTS phibase
    AS
    SELECT
        name,
        checksum,
        GROUP_CONCAT(phibase_ids, ',') as phibase_ids,
        GROUP_CONCAT(phibase_gene, ',') as phibase_genes,
        GROUP_CONCAT(phibase_phenotypes, ',') as phibase_phenotypes
    FROM (
        SELECT DISTINCT
            name,
            checksum,
            REPLACE(phibase_ids, '__', ',') as phibase_ids,
            REPLACE(phibase_gene, '__', ',') as phibase_gene,
            REPLACE(
                SUBSTR(remaining, INSTR(remaining, '#') + 1),
                '__',
                ','
            ) as phibase_phenotypes
        FROM (
            SELECT
                name,
                checksum,
                phibase_ids,
                phibase_gene,
                SUBSTR(remaining, INSTR(remaining, '#') + 1) as remaining
            FROM (
                SELECT
                    name,
                    checksum,
                    phibase_ids,
                    SUBSTR(
                        remaining,
                        1,
                        INSTR(remaining, '#') - 1
                    ) as phibase_gene,
                    SUBSTR(remaining, INSTR(remaining, '#') + 1) as remaining
                FROM (
                    SELECT
                        name,
                        checksum,
                        uniprot_id,
                        SUBSTR(
                            remaining,
                            1,
                            INSTR(remaining, '#') - 1
                        ) as phibase_ids,
                        SUBSTR(
                            remaining,
                            INSTR(remaining, '#') + 1
                        ) as remaining
                    FROM (
                        SELECT
                            name,
                            checksum,
                            SUBSTR(
                                target,
                                1,
                                INSTR(target, '#') - 1
                            ) as uniprot_id,
                            SUBSTR(target, INSTR(target, '#') + 1) as remaining
                        FROM mmseqs
                        where analysis = 'phibase'
                    )
                )
            )
        )
    ) GROUP BY name, checksum
    """)

    cur.execute("DROP TABLE IF EXISTS pfam_targets")
    cur.execute("CREATE TEMP TABLE pfam_targets (pfam_ids text NOT NULL)")
    cur.executemany(
        "INSERT INTO pfam_targets VALUES (?)",
        [(c,) for c in pfam_targets]
    )

    cur.execute(f"""
    CREATE TEMP TABLE pfamscan
    AS
    SELECT
        name,
        checksum,
        MAX(virulence_related) as has_pfam_virulence_match,
        GROUP_CONCAT(hmm, ',') as pfam_ids,
        GROUP_CONCAT(hmm_name, ',') as pfam_names
    FROM (
        SELECT DISTINCT
            json_extract(data, '$.{PfamScan.name_column}') as name,
            checksum,
            (
                SUBSTR(
                    json_extract(data, '$.hmm'),
                    1,
                    INSTR(json_extract(data, '$.hmm') || '.', '.') - 1
                )
                IN pfam_targets
            ) as virulence_related,
            SUBSTR(
                json_extract(data, '$.hmm'),
                1,
                INSTR(json_extract(data, '$.hmm') || '.', '.') - 1
            ) as hmm,
            (json_extract(data, '$.hmm_type') ||
                ':' ||
                json_extract(data, '$.hmm_name')) as hmm_name
        FROM
            results
        WHERE
            analysis = 'pfamscan'
        ORDER BY
            json_extract(data, '$.{PfamScan.name_column}'),
            checksum,
            json_extract(data, '$.evalue') ASC
    )
    GROUP BY
        name,
        checksum
    """)

    cur.execute(f"""  -- # noqa
    CREATE TEMP TABLE IF NOT EXISTS pepstats
    AS
    SELECT
        name,
        checksum,
        molecular_weight,
        residue_number,
        charge,
        isoelectric_point,
        aa_c_number,
        aa_tiny_number,
        aa_small_number,
        aa_aliphatic_number,
        aa_aromatic_number,
        aa_nonpolar_number,
        aa_charged_number,
        aa_basic_number,
        aa_acidic_number,
        ((f_number + k_number + y_number + i_number + n_number + 1.0) /
        (1.0 * (g_number + a_number + p_number + 1.0))) as fykin_gap
    FROM (
        SELECT
            json_extract(data, '$.{PepStats.name_column}') as name,
            checksum,
            json_extract(data, '$.molecular_weight') as molecular_weight,
            json_extract(data, '$.residues') as residue_number,
            json_extract(data, '$.charge') as charge,
            json_extract(data, '$.isoelectric_point') as isoelectric_point,
            json_extract(data, '$.residue_c_number') as aa_c_number,
            json_extract(data, '$.property_tiny_number') as aa_tiny_number,
            json_extract(data, '$.property_small_number') as aa_small_number,
            json_extract(data, '$.property_aliphatic_number') as aa_aliphatic_number,
            json_extract(data, '$.property_aromatic_number') as aa_aromatic_number,
            json_extract(data, '$.property_nonpolar_number') as aa_nonpolar_number,
            json_extract(data, '$.property_charged_number') as aa_charged_number,
            json_extract(data, '$.property_basic_number') as aa_basic_number,
            json_extract(data, '$.property_acidic_number') as aa_acidic_number,
            json_extract(data, '$.residue_f_number') as f_number,
            json_extract(data, '$.residue_k_number') as k_number,
            json_extract(data, '$.residue_y_number') as y_number,
            json_extract(data, '$.residue_i_number') as i_number,
            json_extract(data, '$.residue_n_number') as n_number,
            json_extract(data, '$.residue_g_number') as g_number,
            json_extract(data, '$.residue_a_number') as a_number,
            json_extract(data, '$.residue_p_number') as p_number
        FROM results
        WHERE
            analysis = 'pepstats'
    )
    """)

    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS kex2_cutsite
    AS
    SELECT
        name,
        checksum,
        GROUP_CONCAT(pos, ',') as kex2_cutsites
    FROM (
        SELECT DISTINCT
            json_extract(data, '$.{RegexAnalysis.name_column}') as name,
            checksum,
            (json_extract(data, '$.pattern') || ':' ||
            (json_extract(data, '$.start') + 1) || '-' ||
            json_extract(data, '$.end')) as pos
        FROM results
        WHERE analysis = 'kex2_cutsite'
        ORDER BY json_extract(data, '$.start')
    )
    GROUP BY name, checksum
    """)

    cur.execute(f"""
    CREATE TEMP TABLE IF NOT EXISTS rxlr_like_motif
    AS
    SELECT
        name,
        checksum,
        GROUP_CONCAT(pos, ',') as rxlr_like_motifs
    FROM (
        SELECT DISTINCT
            json_extract(data, '$.{RegexAnalysis.name_column}') as name,
            checksum,
            ((json_extract(data, '$.start') + 1) || '-' ||
            json_extract(data, '$.end')) as pos
        FROM results
        WHERE analysis = 'rxlr_like_motif'
        ORDER BY json_extract(data, '$.start')
    )
    GROUP BY name, checksum
    """)

    print("created subtables")

    cur.execute("""
    CREATE TEMP VIEW all_records
    AS
    SELECT DISTINCT name, checksum
    FROM (
        SELECT name, checksum FROM deeploc
        UNION
        SELECT name, checksum FROM apoplastp
        UNION
        SELECT name, checksum FROM effectorp1
        UNION
        SELECT name, checksum FROM effectorp2
        UNION
        SELECT name, checksum FROM effectorp3
        UNION
        SELECT name, checksum FROM localizer
        UNION
        SELECT name, checksum FROM signalp3_nn
        UNION
        SELECT name, checksum FROM signalp3_hmm
        UNION
        SELECT name, checksum FROM signalp4
        UNION
        SELECT name, checksum FROM signalp5
        UNION
        SELECT name, checksum FROM signalp6
        UNION
        SELECT name, checksum FROM targetp
        UNION
        SELECT name, checksum FROM deepsig
        UNION
        SELECT name, checksum FROM tmhmm
        UNION
        SELECT name, checksum FROM deepredeff_fungi
        UNION
        SELECT name, checksum FROM deepredeff_oomycete
        UNION
        SELECT name, checksum FROM phobius
        UNION
        SELECT name, checksum FROM effectordb
        UNION
        SELECT name, checksum FROM dbcan
        UNION
        SELECT name, checksum FROM phibase
        UNION
        SELECT name, checksum FROM pfamscan
        UNION
        SELECT name, checksum FROM pepstats
        UNION
        SELECT name, checksum FROM kex2_cutsite
        UNION
        SELECT name, checksum FROM rxlr_like_motif
    )
    ORDER BY name, checksum
    """)

    print("created views")

    table = pd.read_sql_query(""" --  # noqa
    SELECT
        all_records.name,
        all_records.checksum,
        phibase.phibase_ids,
        phibase.phibase_genes,
        phibase.phibase_phenotypes,
        effectordb.effector_matches,
        pfamscan.pfam_ids,
        pfamscan.pfam_names,
        pfamscan.has_pfam_virulence_match,
        dbcan.dbcan_matches,
        dbcan.has_dbcan_virulence_match,
        effectorp1.effectorp1,
        effectorp2.effectorp2,
        effectorp3.effectorp3_cytoplasmic,
        effectorp3.effectorp3_apoplastic,
        effectorp3.effectorp3_noneffector,
        deepredeff_fungi.deepredeff_fungi,
        deepredeff_oomycete.deepredeff_oomycete,
        apoplastp.apoplastp,
        pepstats.molecular_weight,
        pepstats.residue_number,
        pepstats.charge,
        pepstats.isoelectric_point,
        pepstats.aa_c_number,
        pepstats.aa_tiny_number,
        pepstats.aa_small_number,
        pepstats.aa_aliphatic_number,
        pepstats.aa_aromatic_number,
        pepstats.aa_nonpolar_number,
        pepstats.aa_charged_number,
        pepstats.aa_basic_number,
        pepstats.aa_acidic_number,
        pepstats.fykin_gap,
        kex2_cutsite.kex2_cutsites,
        rxlr_like_motif.rxlr_like_motifs,
        localizer.localizer_nucleus,
        localizer.localizer_chloro,
        localizer.localizer_mito,
        signalp3_nn.signalp3_nn,
        signalp3_hmm.signalp3_hmm,
        signalp4.signalp4,
        signalp5.signalp5,
        signalp6.signalp6,
        deepsig.deepsig,
        phobius.phobius_sp,
        phobius.phobius_tmcount,
        tmhmm.tmhmm_tmcount,
        tmhmm.tmhmm_first_60,
        tmhmm.tmhmm_exp_aa,
        targetp.targetp_secreted,
        targetp.targetp_secreted_prob,
        targetp.targetp_mitochondrial_prob,
        deeploc.deeploc_membrane,
        deeploc.deeploc_nucleus,
        deeploc.deeploc_cytoplasm,
        deeploc.deeploc_extracellular,
        deeploc.deeploc_mitochondrion,
        deeploc.deeploc_cell_membrane,
        deeploc.deeploc_endoplasmic_reticulum,
        deeploc.deeploc_plastid,
        deeploc.deeploc_golgi,
        deeploc.deeploc_lysosome,
        deeploc.deeploc_peroxisome,
        signalp3_nn.signalp3_nn_d,
        signalp3_hmm.signalp3_hmm_s,
        signalp4.signalp4_d,
        signalp5.signalp5_prob,
        signalp6.signalp6_prob
    FROM all_records
    LEFT JOIN phibase
        ON (all_records.checksum = phibase.checksum) AND (all_records.name = phibase.name)
    LEFT JOIN effectordb
        ON (all_records.checksum = effectordb.checksum) AND (all_records.name = effectordb.name)
    LEFT JOIN pfamscan
        ON (all_records.checksum = pfamscan.checksum) AND (all_records.name = pfamscan.name)
    LEFT JOIN dbcan
        ON (all_records.checksum = dbcan.checksum) AND (all_records.name = dbcan.name)
    LEFT JOIN effectorp1
        ON (all_records.checksum = effectorp1.checksum) AND (all_records.name = effectorp1.name)
    LEFT JOIN effectorp2
        ON (all_records.checksum = effectorp2.checksum) AND (all_records.name = effectorp2.name)
    LEFT JOIN effectorp3
        ON (all_records.checksum = effectorp3.checksum) AND (all_records.name = effectorp3.name)
    LEFT JOIN deepredeff_fungi
        ON (all_records.checksum = deepredeff_fungi.checksum) AND (all_records.name = deepredeff_fungi.name)
    LEFT JOIN deepredeff_oomycete
        ON (all_records.checksum = deepredeff_oomycete.checksum) AND (all_records.name = deepredeff_oomycete.name)
    LEFT JOIN apoplastp
        ON (all_records.checksum = apoplastp.checksum) AND (all_records.name = apoplastp.name)
    LEFT JOIN pepstats
        ON (all_records.checksum = pepstats.checksum) AND (all_records.name = pepstats.name)
    LEFT JOIN kex2_cutsite
        ON (all_records.checksum = kex2_cutsite.checksum) AND (all_records.name = kex2_cutsite.name)
    LEFT JOIN rxlr_like_motif
        ON (all_records.checksum = rxlr_like_motif.checksum) AND (all_records.name = rxlr_like_motif.name)
    LEFT JOIN localizer
        ON (all_records.checksum = localizer.checksum) AND (all_records.name = localizer.name)
    LEFT JOIN signalp3_nn
        ON (all_records.checksum = signalp3_nn.checksum) AND (all_records.name = signalp3_nn.name)
    LEFT JOIN signalp3_hmm
        ON (all_records.checksum = signalp3_hmm.checksum) AND (all_records.name = signalp3_hmm.name)
    LEFT JOIN signalp4
        ON (all_records.checksum = signalp4.checksum) AND (all_records.name = signalp4.name)
    LEFT JOIN signalp5
        ON (all_records.checksum = signalp5.checksum) AND (all_records.name = signalp5.name)
    LEFT JOIN signalp6
        ON (all_records.checksum = signalp6.checksum) AND (all_records.name = signalp6.name)
    LEFT JOIN deepsig
        ON (all_records.checksum = deepsig.checksum) AND (all_records.name = deepsig.name)
    LEFT JOIN phobius
        ON (all_records.checksum = phobius.checksum) AND (all_records.name = phobius.name)
    LEFT JOIN tmhmm
        ON (all_records.checksum = tmhmm.checksum) AND (all_records.name = tmhmm.name)
    LEFT JOIN targetp
        ON (all_records.checksum = targetp.checksum) AND (all_records.name = targetp.name)
    LEFT JOIN deeploc
        ON (all_records.checksum = deeploc.checksum) AND (all_records.name = deeploc.name)
    """, con)

    print("loaded table")

    table["has_pfam_virulence_match"] = (
        table["has_pfam_virulence_match"].fillna(0)
    )

    table["has_dbcan_virulence_match"] = (
        table["has_dbcan_virulence_match"].fillna(0)
    )

    sp_gffs = get_gff_records("""
    SELECT * FROM results
    WHERE analysis IN ('signalp3_nn', 'signalp3_hmm', 'signalp4',
                       'signalp5', 'signalp6', 'phobius', 'deepsig')
    """, cur)

    tm_gffs = get_gff_records(
        "SELECT * FROM results WHERE analysis = 'tmhmm'",
        cur
    )

    tm_sp_coverage = get_tm_sp_coverage(sp_gffs, tm_gffs)
    table["tmhmm_first_tm_sp_coverage"] = table[["name", "checksum"]].apply(
        lambda r: tm_sp_coverage.get((r["name"], r["checksum"]), 0.0),
        axis=1
    )

    sp_sites = fetch_sp_sites(sp_gffs)
    table["signal_peptide_cutsites"] = table[["name", "checksum"]].apply(
        lambda r: sp_sites.get((r["name"], r["checksum"]), None),
        axis=1
    )

    # mutates
    decide_any_signal(table)

    # mutates
    decide_is_transmembrane(
        table,
        tmhmm_first_60_threshold=tmhmm_first_60_threshold
    )

    # mutates
    decide_is_secreted(table)

    # mutates
    get_phibase_cols(table)

    return table


def get_gff_records(
    query: str,
    cur: sqlite3.Cursor,
) -> List[Tuple[str, str, GFFRecord]]:
    out = []
    for r in cur.execute(query):
        row = ResultRow(*r)
        an = row.as_analysis()
        assert isinstance(an, GFFAble), f"{an.__class__.__name__} not gffable"
        for gffrow in an.as_gff():
            gffrow.source = an.__class__.__name__
            out.append((getattr(an, an.name_column), row.checksum, gffrow))
    return out


def fetch_sp_sites(
    gff: Sequence[Tuple[str, str, GFFRecord]]
) -> Dict[Tuple[str, str], str]:
    d = defaultdict(list)

    for name, chk, gffrow in gff:
        if gffrow.type != "signal_peptide":
            continue
        d[(name, chk)].append(gffrow)

    out = dict()
    for (name, chk), gffrows in d.items():
        str_gffrows = ",".join([
            f"{r.source}:{r.end}"
            for r
            in sorted(gffrows, key=lambda k: k.end)
        ])
        out[(name, chk)] = str_gffrows

    return out


def decide_any_signal(
    table: pd.DataFrame
) -> None:
    table["any_signal_peptide"] = (
        table.loc[:, [
            'signalp3_nn', 'signalp3_hmm', 'signalp4', 'signalp5',
            'signalp6', 'deepsig', 'phobius_sp'
        ]]
        .copy()
        .fillna(0)
        .astype(bool)
        .any(axis=1)
        .astype(int)
    )
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
    sp_gff: Sequence[Tuple[str, str, GFFRecord]],
    tm_gff: Sequence[Tuple[str, str, GFFRecord]],
) -> Dict[Tuple[str, str], float]:
    sp_d = defaultdict(list)
    tm_d = defaultdict(list)

    for name, chk, gffrow in sp_gff:
        if gffrow.type != "signal_peptide":
            continue
        sp_d[(name, chk)].append(gffrow)

    for name, chk, gffrow in tm_gff:
        tm_d[(name, chk)].append(gffrow)

    all_keys = set(sp_d.keys())
    all_keys.update(tm_d.keys())

    out = dict()
    for name, chk in all_keys:
        sp_gffrows = sp_d.get((name, chk), [])
        tm_gffrows = tm_d.get((name, chk), [])

        if len(sp_gffrows) == 0:
            out[(name, chk)] = 0.0
            continue

        elif len(tm_gffrows) == 0:
            out[(name, chk)] = 0.0
            continue

        tm = sorted(tm_gffrows, key=lambda x: x.start)[0]
        covs = [gff_coverage(sp, tm) for sp in sp_gffrows]
        out[(name, chk)] = median(covs)

    return out


def decide_is_transmembrane(
    table: pd.DataFrame,
    tmhmm_first_60_threshold: float = 10,
) -> None:

    table["multiple_transmembrane"] = (
        (table["tmhmm_tmcount"] > 1) |
        (table["phobius_tmcount"] > 1)
    ).astype(int)

    table["single_transmembrane"] = (
        ~table["multiple_transmembrane"].astype(bool)
        & (
            table["phobius_tmcount"] == 1
            | (
                (table["tmhmm_tmcount"] == 1)
                & ~table["any_signal_peptide"].astype(bool)
            )
            | (
                (table["tmhmm_tmcount"] == 1)
                & table["any_signal_peptide"].astype(bool)
                & (table["tmhmm_first_60"] < tmhmm_first_60_threshold)
            )
        )
    ).astype(int)

    return


def decide_is_secreted(
    table: pd.DataFrame
) -> None:
    table["is_secreted"] = (
        table["any_signal_peptide"].astype(bool) & ~
        table["multiple_transmembrane"].astype(bool)
    ).astype(int)
    return


def get_phibase_cols(table: pd.DataFrame) -> None:
    effector_phenotypes = {
            "loss_of_pathogenicity",
            "increased_virulence_(hypervirulence)",
            "effector_(plant_avirulence_determinant)"
    }

    phenotypes = (
        table["phibase_phenotypes"].fillna('').str.split(',')
        .apply(set)
    )

    table["has_phibase_effector_match"] = (
        phenotypes
        .apply(lambda s: len(s.intersection(effector_phenotypes)) > 0)
        .astype(int)
    )

    table["has_phibase_virulence_match"] = (
        phenotypes
        .apply(lambda s: "reduced_virulence" in s)
        .astype(int)
    )

    table["has_phibase_lethal_match"] = (
        phenotypes
        .apply(lambda s: "lethal" in s)
        .astype(int)
    )

    return


def effector_score_it(
    table: pd.DataFrame,
    effectorp1: float = 3,
    effectorp2: float = 3,
    effectorp3_apoplastic: float = 3,
    effectorp3_cytoplasmic: float = 3,
    effectorp3_noneffector: float = -3,
    deepredeff_fungi: float = 2,
    deepredeff_oomycete: float = 2,
    effector: float = 5,
    virulence: float = 2,
    lethal: float = -5,
) -> pd.Series:
    """ """
    score = table["manual_secretion_score"].copy()

    score += (
        table["effectorp1"]
        .fillna(0.5)
        .astype(float)
        - 0.5
    ) * 2 * effectorp1

    score += (
        table["effectorp2"]
        .fillna(0.5)
        .astype(float)
        - 0.5
    ) * 2 * effectorp2

    score += (
        table["effectorp3_apoplastic"]
        .fillna(0.0)
        .astype(float)
        * effectorp3_apoplastic
    )

    score += (
        table["effectorp3_cytoplasmic"]
        .fillna(0.0)
        .astype(float)
        * effectorp3_cytoplasmic
    )

    score += (
        table["effectorp3_noneffector"]
        .fillna(0.0)
        .astype(float)
        * effectorp3_noneffector
    )

    score += (
        table["deepredeff_fungi"]
        .fillna(0.5)
        .astype(float)
        - 0.5
    ) * 2 * deepredeff_fungi

    score += (
        table["deepredeff_oomycete"]
        .fillna(0.5)
        .astype(float)
        - 0.5
    ) * 2 * deepredeff_oomycete

    has_effector_match = (
        table["has_phibase_effector_match"].fillna(0).astype(bool) |
        table["effector_matches"].notnull().astype(bool) |
        table["has_dbcan_virulence_match"].fillna(0).astype(bool) |
        table["has_pfam_virulence_match"].fillna(0).astype(bool)
    ).astype(int)

    score += has_effector_match * effector

    score += (
        (~has_effector_match) *
        table["has_phibase_virulence_match"].fillna(0).astype(int) *
        virulence
    )

    score += table["has_phibase_lethal_match"].fillna(0).astype(int) * lethal
    return score


def secretion_score_it(
    table,
    secreted: float = 3.,
    less_trustworthy_signal_prediction: float = 0.25,
    trustworthy_signal_prediction: float = 0.5,
    single_transmembrane: float = -1,
    multiple_transmembrane: float = -10,
    deeploc_extracellular: float = 1,
    deeploc_intracellular: float = -2,
    deeploc_membrane: float = -2,
    targetp_mitochondrial: float = -2,
) -> pd.Series:
    score = 1.0 * table["is_secreted"].astype(int).fillna(0) * secreted

    for k in ["signalp3_hmm", "signalp3_nn", "phobius_sp", "deepsig"]:
        col = table[k].astype(int).fillna(0)
        score += 1.0 * col * less_trustworthy_signal_prediction

    for k in ["signalp4", "signalp5", "signalp6", "targetp_secreted"]:
        col = table[k].astype(int).fillna(0)
        score += 1.0 * col * trustworthy_signal_prediction

    score += (
        1.0 *
        table["multiple_transmembrane"].astype(int).fillna(0) *
        multiple_transmembrane
    )

    score += (
        1.0 *
        table["single_transmembrane"].astype(int).fillna(0) *
        single_transmembrane
    )

    score += (
        table["deeploc_extracellular"].astype(float).fillna(0.0) *
        deeploc_extracellular
    )

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
        col = table[k].astype(float).fillna(0.0)
        score += col * deeploc_intracellular

    score += (
        table["deeploc_membrane"].astype(float).fillna(0.0) *
        deeploc_membrane
    )

    score += (
        table["targetp_mitochondrial_prob"].astype(float).fillna(0.0) *
        targetp_mitochondrial
    )
    return score.astype(float)


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


def inner(con: sqlite3.Connection, args: argparse.Namespace) -> None:
    cur = con.cursor()

    tab = ResultsTable(con, cur)
    tab.create_tables()

    tab.insert_results(
        ResultRow.from_file(args.infile, replace_name=False)
    )
    print("loaded results")

    if args.dbcan is not None:
        dbcan: Set[str] = {d.strip() for d in args.dbcan.readlines()}
    else:
        dbcan = set(get_interesting_dbcan_ids())

    if args.pfam is not None:
        pfam: Set[str] = {d.strip() for d in args.pfam.readlines()}
    else:
        pfam = set(get_interesting_pfam_ids())

    df = create_tables(
        tab,
        con,
        cur,
        pfam,
        dbcan,
        args.tmhmm_first_60_threshold
    )

    df["manual_secretion_score"] = secretion_score_it(
        df,
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

    df["manual_effector_score"] = effector_score_it(
        df,
        args.effectorp1_weight,
        args.effectorp2_weight,
        args.effectorp3_apoplastic_weight,
        args.effectorp3_cytoplasmic_weight,
        args.effectorp3_noneffector_weight,
        args.deepredeff_fungi_weight,
        args.deepredeff_oomycete_weight,
        args.effector_homology_weight,
        args.virulence_homology_weight,
        args.lethal_homology_weight,
    )

    df["effector_score"] = run_ltr(df)
    df.sort_values("effector_score", ascending=False, inplace=True)
    df[COLUMNS].round(3).to_csv(
        args.outfile,
        sep="\t",
        index=False,
        na_rep=".",
    )
    return


def runner(args: argparse.Namespace) -> None:
    con = sqlite3.connect(args.db)
    try:
        inner(con, args)
    except Exception as e:
        raise e
    finally:
        con.commit()
        con.close()
    return
