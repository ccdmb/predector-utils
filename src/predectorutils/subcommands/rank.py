#!/usr/bin/env python3

import sys
import argparse

from math import floor
from statistics import median

from typing import (Optional, Tuple, Dict, Set, List)

import sqlite3

import numpy as np
import pandas as pd
import xgboost as xgb

from predectorutils.analyses import Analyses
from predectorutils.database import ResultsTable
from predectorutils.data import (
    get_interesting_dbcan_ids,
    get_interesting_pfam_ids,
    get_ltr_model,
)


def cli(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "db",
        type=str,
        help="Where the sqlite database is."
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Where to write the output to. Default: stdout"
    )

    parser.add_argument(
        "--mem",
        type=float,
        default=1.0,
        help=(
            "The amount of RAM in gibibytes to let "
            "SQLite use for cache."
        )
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


def create_pfam_targets(
    con: sqlite3.Connection,
    cur: sqlite3.Cursor,
    targets: Set[str]
) -> None:
    cur.execute("DROP TABLE IF EXISTS pfam_targets")
    cur.execute(
        "CREATE TEMP TABLE pfam_targets (pfam_ids text NOT NULL UNIQUE)"
    )
    cur.executemany(
        "INSERT INTO pfam_targets VALUES (?)",
        ((c,) for c in targets)
    )
    cur.execute(
        "CREATE UNIQUE INDEX pfam_targets_index ON pfam_targets (pfam_ids)"
    )
    con.commit()
    return


def create_dbcan_targets(
    con: sqlite3.Connection,
    cur: sqlite3.Cursor,
    targets: Set[str]
) -> None:
    cur.execute("DROP TABLE IF EXISTS dbcan_targets")
    cur.execute(
        "CREATE TEMP TABLE dbcan_targets (dbcan_ids text NOT NULL UNIQUE)"
    )
    cur.executemany(
        "INSERT INTO dbcan_targets VALUES (?)",
        ((c,) for c in targets)
    )
    cur.execute(
        "CREATE UNIQUE INDEX dbcan_targets_index ON dbcan_targets (dbcan_ids)"
    )
    con.commit()
    return


def agg_dbcan_vir() -> str:
    s = (  # noqa
        "IFNULL("
            "MAX(json_extract(data, '$.hmm') IN pfam_targets"
            f") FILTER (WHERE analysis == {int(Analyses.dbcan)}), "
            "0"
        ")"
    )
    return s


def agg_pfam_vir() -> str:
    s = (  # noqa
        "IFNULL("
            "MAX("
                "SUBSTR("
                    "json_extract(data, '$.hmm'), "
                    "1, "
                    "INSTR(json_extract(data, '$.hmm') || '.', '.') - 1"
                ") IN pfam_targets"
            f") FILTER (WHERE analysis == {int(Analyses.pfamscan)}), "
            "0"
        ")"
    )
    return s


class AggBase(object):

    nargs: int
    analysis: Analyses
    sql_fname: str

    def __init__(self) -> None:
        raise NotImplementedError("this is a baseclass")


class AggPhobiusTMDomains(AggBase):

    nargs = 1
    analysis = Analyses.phobius
    sql_fname = "agg_phobius_domains"

    def __init__(self):
        self.matches: List[Tuple[int, int]] = []
        return

    def step(self, data: str) -> None:
        from ..analyses import GFFAble
        an = (
            Analyses(self.analysis)
            .get_analysis()
            .from_json_str(data)
        )

        assert isinstance(an, GFFAble)

        for gffrow in an.as_gff():
            if gffrow.type != "transmembrane_polypeptide_region":
                continue

            self.matches.append((gffrow.start + 1, gffrow.end))
        return

    def finalize(self) -> Optional[str]:
        if len(self.matches) == 0:
            return None

        return "tm:" + ",".join(
            f"{s}-{e}"
            for s, e
            in sorted(self.matches, key=lambda t: t[0])
        )

    @classmethod
    def sql_query(cls) -> str:
        return (
            f"{cls.sql_fname}(data)"
            f"FILTER (WHERE analysis = {int(cls.analysis)})"
        )


class AggTMHMMDomains(AggBase):

    nargs = 1
    analysis = Analyses.tmhmm
    sql_fname = "agg_tmhmm_domains"

    def __init__(self):
        self.matches: List[Tuple[int, int]] = []
        return

    def step(self, data: str) -> None:
        from ..analyses import GFFAble
        an = (
            Analyses(self.analysis)
            .get_analysis()
            .from_json_str(data)
        )

        assert isinstance(an, GFFAble)

        for gffrow in an.as_gff():
            if gffrow.type != "transmembrane_polypeptide_region":
                continue

            self.matches.append((gffrow.start + 1, gffrow.end))
        return

    def finalize(self) -> Optional[str]:
        if len(self.matches) == 0:
            return None

        return "tm:" + ",".join(
            f"{s}-{e}"
            for s, e
            in sorted(self.matches, key=lambda t: t[0])
        )

    @classmethod
    def sql_query(cls) -> str:
        return (
            f"{cls.sql_fname}(data)"
            f"FILTER (WHERE analysis = {int(cls.analysis)})"
        )


class AggTMSPCoverage(AggBase):

    nargs = 2
    analyses = [
        Analyses.signalp3_hmm,
        Analyses.signalp3_nn,
        Analyses.signalp4,
        Analyses.signalp5,
        Analyses.signalp6,
        Analyses.deepsig,
        Analyses.phobius,
        Analyses.tmhmm
    ]
    sql_fname = "agg_tmsp_coverage"

    def __init__(self):
        self.signals: List[int] = []
        self.membranes: List[Tuple[int, int]] = []
        return

    def step(self, analysis: int, data: str) -> None:
        from ..analyses import GFFAble
        an = (
            Analyses(analysis)
            .get_analysis()
            .from_json_str(data)
        )

        assert isinstance(an, GFFAble)

        for gffrow in an.as_gff():
            if gffrow.type == "signal_peptide":
                self.signals.append(gffrow.end)
            elif gffrow.type == "transmembrane_polypeptide_region":
                self.membranes.append((gffrow.start, gffrow.end))
        return

    def finalize(self) -> float:
        if len(self.signals) == 0:
            return 0.0

        if len(self.membranes) == 0:
            return 0.0

        tm = sorted(self.membranes, key=lambda t: t[0])[0]
        covs = [self.gff_coverage(sp, tm) for sp in self.signals]
        return median(covs)

    @classmethod
    def sql_query(cls) -> str:
        sources = ",".join(str(int(a)) for a in cls.analyses)

        return (
            f"{cls.sql_fname}(analysis, data)"
            f"FILTER (WHERE analysis IN ({sources}))"
        )

    def gff_intersection(self, sp: int, tm: Tuple[int, int]) -> int:
        lstart = 0
        lend = sp

        rstart = min(tm)
        rend = max(tm)

        start = max([lstart, rstart])
        end = min([lend, rend])

        # This will be < 0 if they don't overlap
        if start < end:
            return end - start
        else:
            return 0

    def gff_coverage(self, sp: int, tm: Tuple[int, int]) -> float:
        noverlap = self.gff_intersection(sp, tm)
        rstart, rend = tm
        return noverlap / (rend - rstart)


class AggSPCutsite(AggBase):

    nargs = 2
    analyses = [
        Analyses.signalp3_hmm,
        Analyses.signalp3_nn,
        Analyses.signalp4,
        Analyses.signalp5,
        Analyses.signalp6,
        Analyses.deepsig,
        Analyses.phobius
    ]
    sql_fname = "agg_sp_cutsite"

    def __init__(self):
        self.matches: List[Tuple[str, int]] = []
        return

    def step(self, analysis: int, data: str) -> None:
        from ..analyses import GFFAble
        an = (
            Analyses(analysis)
            .get_analysis()
            .from_json_str(data)
        )

        assert isinstance(an, GFFAble)

        for gffrow in an.as_gff():
            if gffrow.type != "signal_peptide":
                continue

            self.matches.append((an.__class__.__name__, gffrow.end))
        return

    def finalize(self) -> Optional[str]:
        if len(self.matches) == 0:
            return None

        return ",".join(
            f"{n}:{e}"
            for n, e
            in sorted(self.matches, key=lambda t: t[1])
        )

    @classmethod
    def sql_query(cls) -> str:
        sources = ",".join(str(int(a)) for a in cls.analyses)

        return (
            f"{cls.sql_fname}(analysis, data)"
            f"FILTER (WHERE analysis IN ({sources}))"
        )


class AggKex2(AggBase):
    nargs = 4
    analysis = Analyses.kex2_cutsite
    sql_fname = "agg_kex2"

    def __init__(self):
        self.matches: Dict[Tuple[str, int, int], Set[str]] = dict()
        return

    def step(
        self,
        pattern: str,
        match: str,
        start: int,
        end: int
    ) -> None:
        pattern = pattern.replace("[A-Z]", "X")
        tup = (match, start + 1, end)
        if tup in self.matches:
            self.matches[tup].add(pattern)
        else:
            self.matches[tup] = {pattern}
        return

    def finalize(self) -> Optional[str]:
        if len(self.matches) == 0:
            return None
        else:
            return ",".join(
                f"{m}:{'&'.join(ps)}:{s}-{e}"
                for (m, s, e), ps
                in sorted(self.matches.items(), key=lambda t: t[0][1])
            )

    @classmethod
    def sql_query(cls) -> str:
        return (
            f"{cls.sql_fname}("
            "json_extract(data, '$.pattern'), "
            "json_extract(data, '$.match'), "
            "json_extract(data, '$.start'), "
            "json_extract(data, '$.end')"
            ") "
            f"FILTER (WHERE analysis = {int(cls.analysis)})"
        )


class AggRxLR(AggBase):
    nargs = 3
    analysis = Analyses.rxlr_like_motif
    sql_fname = "agg_rxlr"

    def __init__(self):
        self.matches: List[Tuple[str, int, int]] = []
        return

    def step(
        self,
        match: str,
        start: int,
        end: int
    ) -> None:
        self.matches.append((match, start + 1, end))
        return

    def finalize(self) -> Optional[str]:
        if len(self.matches) == 0:
            return None
        else:
            return ",".join(
                f"{m}:{s}-{e}"
                for (m, s, e)
                in sorted(self.matches, key=lambda t: t[1])
            )

    @classmethod
    def sql_query(cls) -> str:
        return (
            f"{cls.sql_fname}("
            "json_extract(data, '$.match'), "
            "json_extract(data, '$.start'), "
            "json_extract(data, '$.end')"
            ") "
            f"FILTER (WHERE analysis = {int(cls.analysis)})"
        )


class AggPHIMatches(AggBase):
    nargs = 1
    index: int

    def __init__(self):
        self.matches: Dict[str, float] = dict()
        return

    def step(self, data: str) -> None:
        from ..analyses import MMSeqs
        an = (
            self.analysis
            .get_analysis()
            .from_json_str(data)
        )

        assert isinstance(an, MMSeqs)

        if an.decide_significant():
            si = an.target.strip().split("#")
            assert len(si) == 6
            matches = set(si[self.index].split("__"))
            for match in matches:
                if match in self.matches:
                    if an.evalue < self.matches[match]:
                        self.matches[match] = an.evalue
                else:
                    self.matches[match] = an.evalue
        return

    def finalize(self) -> Optional[str]:
        if len(self.matches) == 0:
            return None
        else:
            return ",".join(
                k
                for k, _
                in sorted(self.matches.items(), key=lambda t: t[1])
            )

    @classmethod
    def sql_query(cls) -> str:
        return (
            f"{cls.sql_fname}(data) "
            f"FILTER (WHERE analysis = {int(cls.analysis)})"
        )


class AggPHIPhenos(AggPHIMatches):
    nargs = 1
    sql_fname = "agg_phi_phenos"
    analysis = Analyses.phibase
    index = 5


class AggPHIIDs(AggPHIMatches):
    nargs = 1
    sql_fname = "agg_phi_ids"
    analysis = Analyses.phibase
    index = 1


class AggPHIGenes(AggPHIMatches):
    nargs = 1
    sql_fname = "agg_phi_genes"
    analysis = Analyses.phibase
    index = 2


class AggPHIHasMatch(AggBase):
    nargs = 1
    analysis = Analyses.phibase
    index = 5
    targets: Set[str]

    def __init__(self):
        self.phenos: Set[str] = set()
        return

    def step(self, data: str) -> None:
        from ..analyses import MMSeqs
        an = (
            self.analysis
            .get_analysis()
            .from_json_str(data)
        )

        assert isinstance(an, MMSeqs)
        if an.decide_significant():
            si = an.target.strip().split("#")
            assert len(si) == 6
            phenos = set(si[self.index].lower().split("__"))
            self.phenos.update(phenos)
        return

    def finalize(self) -> bool:
        if len(self.phenos) == 0:
            return False
        return len(
            self.phenos.intersection(self.targets)
        ) > 0

    @classmethod
    def sql_query(cls) -> str:
        return (
            "IFNULL("
            f"{cls.sql_fname}(data) "
            f"FILTER (WHERE analysis = {int(cls.analysis)})"
            ", 0)"
        )


class AggPHIEffectorMatch(AggPHIHasMatch):
    sql_fname = "agg_phi_ematch"
    targets = {
        "loss_of_pathogenicity",
        "increased_virulence_(hypervirulence)",
        "effector_(plant_avirulence_determinant)"
    }


class AggPHIVirulenceMatch(AggPHIHasMatch):
    sql_fname = "agg_phi_vmatch"
    targets = {"reduced_virulence"}


class AggPHILethalMatch(AggPHIHasMatch):
    sql_fname = "agg_phi_lmatch"
    targets = {"lethal"}


class AggHMMER(AggBase):
    nargs = 1

    def __init__(self):
        self.matches: Dict[str, float] = {}
        return

    def step(self, data: str) -> None:
        from ..analyses import DomTbl
        an = (
            self.analysis
            .get_analysis()
            .from_json_str(data)
        )

        assert isinstance(an, DomTbl)

        if an.decide_significant():
            if an.hmm in self.matches:
                if self.matches[an.hmm] > an.domain_i_evalue:
                    self.matches[an.hmm] = an.domain_i_evalue
            else:
                self.matches[an.hmm] = an.domain_i_evalue

        return

    def finalize(self) -> Optional[str]:
        if len(self.matches) == 0:
            return None

        return ",".join(
            k
            for k, _
            in sorted(self.matches.items(), key=lambda t: t[1])
        )

    @classmethod
    def sql_query(cls) -> str:
        return (
            f"{cls.sql_fname}(data) "
            f"FILTER (WHERE analysis = {int(cls.analysis)})"
        )


class AggEffectorDB(AggHMMER):
    analysis = Analyses.effectordb
    sql_fname = "agg_effdb"


class AggDBCAN(AggHMMER):
    analysis = Analyses.dbcan
    sql_fname = "agg_dbcan"


class AggPfamscanIDS(AggBase):
    nargs = 1
    analysis = Analyses.pfamscan
    sql_fname = "agg_pfamscan_ids"

    def __init__(self):
        self.matches: Dict[str, float] = {}
        return

    def step(self, data: str) -> None:
        from ..analyses import PfamScan
        an = (
            self.analysis
            .get_analysis()
            .from_json_str(data)
        )

        assert isinstance(an, PfamScan)
        hmm = self.split_hmm(an.hmm)

        if hmm in self.matches:
            if self.matches[hmm] > an.evalue:
                self.matches[hmm] = an.evalue
        else:
            self.matches[hmm] = an.evalue
        return

    def finalize(self) -> Optional[str]:
        if len(self.matches) == 0:
            return None

        return ",".join(
            k
            for k, _
            in sorted(self.matches.items(), key=lambda t: t[1])
        )

    @classmethod
    def sql_query(cls) -> str:
        return (
            f"{cls.sql_fname}(data) "
            f"FILTER (WHERE analysis = {int(cls.analysis)})"
        )

    @staticmethod
    def split_hmm(hmm: str) -> str:
        return hmm.strip().split(".")[0]


class AggPfamscanNames(AggBase):
    nargs = 1
    analysis = Analyses.pfamscan
    sql_fname = "agg_pfamscan_names"

    def __init__(self):
        self.matches: Dict[str, float] = {}
        return

    def step(self, data: str) -> None:
        from ..analyses import PfamScan
        an = (
            self.analysis
            .get_analysis()
            .from_json_str(data)
        )

        assert isinstance(an, PfamScan)
        name = an.hmm_type + ":" + an.hmm_name

        if name in self.matches:
            if self.matches[name] > an.evalue:
                self.matches[name] = an.evalue
        else:
            self.matches[name] = an.evalue
        return

    def finalize(self) -> Optional[str]:
        if len(self.matches) == 0:
            return None

        return ",".join(
            k
            for k, _
            in sorted(self.matches.items(), key=lambda t: t[1])
        )

    @classmethod
    def sql_query(cls) -> str:
        return (
            f"{cls.sql_fname}(data) "
            f"FILTER (WHERE analysis = {int(cls.analysis)})"
        )


class AggSperProb(AggBase):
    nargs = 2
    pred_col: str
    prob_col: str
    pos_values: List[str]

    def __init__(self):
        self.prob: Optional[float] = None
        return

    def step(self, prob: float, pred: str) -> None:
        if (self.prob is not None) and (self.prob != prob):
            raise ValueError("This shouldn't happen!")

        if pred in self.pos_values:
            self.prob = prob
        else:
            self.prob = 1 - prob
        return

    def finalize(self) -> Optional[float]:
        return self.prob

    @classmethod
    def sql_query(cls) -> str:
        return (  # noqa
            f"{cls.sql_fname}("
                f"json_extract(data, '$.{cls.prob_col}'), "
                f"json_extract(data, '$.{cls.pred_col}')"
            f') FILTER (WHERE analysis = {int(cls.analysis)})'
        )


class AggEP1(AggSperProb):

    sql_fname = "agg_ep1"
    prob_col = "prob"
    pred_col = "prediction"
    analysis = Analyses.effectorp1

    pos_values = ["Effector"]


class AggEP2(AggSperProb):

    sql_fname = "agg_ep2"
    prob_col = "prob"
    pred_col = "prediction"
    analysis = Analyses.effectorp2

    pos_values = ["Effector", "Unlikely effector"]


class AggApoplastP(AggSperProb):

    sql_fname = "agg_apoplastp"
    prob_col = "prob"
    pred_col = "prediction"
    analysis = Analyses.apoplastp

    pos_values = ["Apoplastic"]


def agg_json(field: str, analysis: str) -> str:
    an = Analyses.from_string(analysis)
    return (
        f"MAX(json_extract(data, '$.{field}')) "
        f"FILTER (WHERE analysis = {int(an)})"
    )


# This variant makes sure it's between 0 and 1
def agg_json_prob(field: str, analysis: str) -> str:
    an = Analyses.from_string(analysis)

    filtered = (
        f"(MAX(json_extract(data, '$.{field}')) "
        f"FILTER (WHERE analysis = {int(an)}))"
    )

    return f"MAX(MIN(({filtered}), 1.0), 0.0)"


# Checks if it is a value
def agg_json_eq_str(field: str, value: str, analysis: str) -> str:
    an = Analyses.from_string(analysis)

    return (
        f"MAX(json_extract(data, '$.{field}') == '{value}') "
        f"FILTER (WHERE analysis = {int(an)})"
    )


def agg_deepsig(pred: str) -> str:
    an = Analyses.deepsig

    filtered = (  # noqa
        "MAX(json_extract(data, '$.prob'))"
        "FILTER ("
            f"WHERE analysis = {int(an)} "
            f"AND json_extract(data, '$.prediction') == '{pred}'"
        ")"
    )

    return f"MAX(MIN(({filtered}), 1.0), 0.0)"


def agg_fkyin_gap() -> str:
    an = Analyses.pepstats

    template = (
        "(MAX(json_extract(data, '$.{field}')) "
        f"FILTER (WHERE analysis = {int(an)}))"
    )

    col_template = "residue_{}_number"

    numerator = "+".join([
        template.format(field=col_template.format(b))
        for b
        in "fkyin"
    ])

    denominator = " + ".join([
        template.format(field=col_template.format(b))
        for b
        in "gap"
    ])

    return f"((1.0 * {numerator} + 1) / (1.0 * {denominator} + 1))"


def load_db(
    path: str,
    mem: float
) -> Tuple[sqlite3.Connection, sqlite3.Cursor]:
    sqlite3.register_converter("analyses", Analyses.from_bytes_)
    con = sqlite3.connect(path, detect_types=sqlite3.PARSE_DECLTYPES)
    con.row_factory = sqlite3.Row

    for c in [
        AggEP1, AggEP2, AggApoplastP, AggEffectorDB, AggDBCAN,
        AggPHIPhenos, AggPHIIDs, AggPHIGenes,
        AggPHIEffectorMatch, AggPHIVirulenceMatch, AggPHILethalMatch,
        AggPfamscanIDS, AggPfamscanNames, AggKex2, AggRxLR,
        AggSPCutsite, AggTMSPCoverage, AggTMHMMDomains, AggPhobiusTMDomains
    ]:
        con.create_aggregate(c.sql_fname, c.nargs, c)  # type: ignore

    cur = con.cursor()
    # Allow it to use 1GB RAM for cache
    mem_gb: int = floor(1000000 * mem)
    cur.execute(f"PRAGMA cache_size = -{mem_gb}")
    cur.execute("PRAGMA journal_mode = WAL")
    cur.execute("PRAGMA locking_mode = EXCLUSIVE")
    cur.execute("PRAGMA synchronous = NORMAL")

    con.commit()
    return con, cur


def create_select_all_table(tab: ResultsTable) -> None:
    assert tab.exists_table("decoder"), "no decoder table"

    query = []
    for an in Analyses:
        ans = str(an)
        anv = an.value

        if not tab.exists_table(f"results_{ans}"):
            continue

        if an.needs_database():
            db_cols = (
                """
                database,
                database_version,
                """
            )
        else:
            db_cols = (
                """
                '' as database,
                '' as database_version,
                """
            )

        q = f"""
        SELECT
            CAST({anv} AS analyses) as analysis,
            software,
            software_version,
            {db_cols}
            pipeline_version,
            r.checksum as checksum,
            d.name as name,
            data
        FROM results_{ans} r
        INNER JOIN decoder d ON r.checksum = d.checksum
        """
        query.append(q)

    qstring = (
        "CREATE TEMP VIEW rank_table AS \n" +
        "\nUNION\n".join(query)
    )
    tab.cur.execute(qstring)
    return


def create_tables(
    con: sqlite3.Connection,
    cur: sqlite3.Cursor,
    pfam_targets: Set[str],
    dbcan_targets: Set[str],
    tmhmm_first_60_threshold: float = 10,
) -> pd.DataFrame:
    tab = ResultsTable(con, cur)
    create_select_all_table(tab)

    table = pd.read_sql_query(
        f""" --  # noqa
        SELECT
            name,
            NULL as effector_score,
            NULL as manual_effector_score,
            NULL as manual_secretion_score,
            {AggEffectorDB.sql_query()} as effector_matches,
            {AggPHIGenes.sql_query()} as phibase_genes,
            {AggPHIPhenos.sql_query()} as phibase_phenotypes,
            {AggPHIIDs.sql_query()} as phibase_ids,
            {AggPHIEffectorMatch.sql_query()} as has_phibase_effector_match,
            {AggPHIVirulenceMatch.sql_query()} as has_phibase_virulence_match,
            {AggPHILethalMatch.sql_query()} as has_phibase_lethal_match,
            {AggPfamscanIDS.sql_query()} as pfam_ids,
            {AggPfamscanNames.sql_query()} as pfam_names,
            {agg_pfam_vir()} as has_pfam_virulence_match,
            {AggDBCAN.sql_query()} as dbcan_matches,
            {agg_dbcan_vir()} as has_dbcan_virulence_match,
            {AggEP1.sql_query()} as effectorp1,
            {AggEP2.sql_query()} as effectorp2,
            {agg_json('cytoplasmic_prob', 'effectorp3')} as effectorp3_cytoplasmic,
            {agg_json('apoplastic_prob', 'effectorp3')} as effectorp3_apoplastic,
            {agg_json('noneffector_prob', 'effectorp3')} as effectorp3_noneffector,
            {agg_json('s_score', 'deepredeff_fungi')} as deepredeff_fungi,
            {agg_json('s_score', 'deepredeff_oomycete')} as deepredeff_oomycete,
            {AggApoplastP.sql_query()} as apoplastp,
            NULL as is_secreted,
            NULL as any_signal_peptide,
            NULL as single_transmembrane,
            NULL as multiple_transmembrane,
            {agg_json('molecular_weight', 'pepstats')} as molecular_weight,
            {agg_json('residues', 'pepstats')} as residue_number,
            {agg_json('charge', 'pepstats')} as charge,
            {agg_json('isoelectric_point', 'pepstats')} as isoelectric_point,
            {agg_json('residue_c_number', 'pepstats')} as aa_c_number,
            {agg_json('property_tiny_number', 'pepstats')} as aa_tiny_number,
            {agg_json('property_small_number', 'pepstats')} as aa_small_number,
            {agg_json('property_aliphatic_number', 'pepstats')} as aa_aliphatic_number,
            {agg_json('property_aromatic_number', 'pepstats')} as aa_aromatic_number,
            {agg_json('property_nonpolar_number', 'pepstats')} as aa_nonpolar_number,
            {agg_json('property_charged_number', 'pepstats')} as aa_charged_number,
            {agg_json('property_basic_number', 'pepstats')} as aa_basic_number,
            {agg_json('property_acidic_number', 'pepstats')} as aa_acidic_number,
            {agg_fkyin_gap()} as fykin_gap,
            {AggKex2.sql_query()} as kex2_cutsites,
            {AggRxLR.sql_query()} as rxlr_like_motifs,
            {agg_json('nucleus_decision', 'localizer')} as localizer_nucleus,
            {agg_json('chloroplast_decision', 'localizer')} as localizer_chloro,
            {agg_json('mitochondria_decision', 'localizer')} as localizer_mito,
            {AggSPCutsite.sql_query()} as signal_peptide_cutsites,
            {agg_json('d_decision', 'signalp3_nn')} as signalp3_nn,
            {agg_json('is_secreted', 'signalp3_hmm')} as signalp3_hmm,
            {agg_json('decision', 'signalp4')} as signalp4,
            {agg_json_eq_str("prediction", "SP(Sec/SPI)", "signalp5")} as signalp5,
            {agg_json_eq_str("prediction", "SP", "signalp6")} as signalp6,
            {agg_json_eq_str("prediction", "SignalPeptide", "deepsig")} as deepsig,
            {agg_json('sp', 'phobius')} as phobius_sp,
            {agg_json('tm', 'phobius')} as phobius_tmcount,
            {AggPhobiusTMDomains.sql_query()} as phobius_tm_domains,
            {agg_json('pred_hel', 'tmhmm')} as tmhmm_tmcount,
            {agg_json('first_60', 'tmhmm')} as tmhmm_first_60,
            {agg_json('exp_aa', 'tmhmm')} as tmhmm_exp_aa,
            {AggTMSPCoverage.sql_query()} as tmhmm_first_tm_sp_coverage,
            {AggTMHMMDomains.sql_query()} as tmhmm_domains,
            {agg_json_eq_str('prediction', 'SP', 'targetp_non_plant')} as targetp_secreted,
            {agg_json_prob('sp', 'targetp_non_plant')} as targetp_secreted_prob,
            {agg_json_prob('mtp', 'targetp_non_plant')} as targetp_mitochondrial_prob,
            {agg_json_prob('membrane', 'deeploc')} as deeploc_membrane,
            {agg_json_prob('nucleus', 'deeploc')} as deeploc_nucleus,
            {agg_json_prob('cytoplasm', 'deeploc')} as deeploc_cytoplasm,
            {agg_json_prob('extracellular', 'deeploc')} as deeploc_extracellular,
            {agg_json_prob('mitochondrion', 'deeploc')} as deeploc_mitochondrion,
            {agg_json_prob('cell_membrane', 'deeploc')} as deeploc_cell_membrane,
            {agg_json_prob('endoplasmic_reticulum', 'deeploc')} as deeploc_endoplasmic_reticulum,
            {agg_json_prob('plastid', 'deeploc')} as deeploc_plastid,
            {agg_json_prob('golgi_apparatus', 'deeploc')} as deeploc_golgi,
            {agg_json_prob('lysosome_vacuole', 'deeploc')} as deeploc_lysosome,
            {agg_json_prob('peroxisome', 'deeploc')} as deeploc_peroxisome,
            {agg_json('d', 'signalp3_nn')} as signalp3_nn_d,
            {agg_json('sprob', 'signalp3_hmm')} as signalp3_hmm_s,
            {agg_json('d', 'signalp4')} as signalp4_d,
            {agg_json_prob('prob_signal', 'signalp5')} as signalp5_prob,
            {agg_json_prob('prob_signal', 'signalp6')} as signalp6_prob,
            {agg_deepsig("SignalPeptide")} as deepsig_signal_prob,
            {agg_deepsig("Transmembrane")} as deepsig_transmembrane_prob,
            {agg_deepsig("Other")} as deepsig_other_prob
        FROM rank_table
        GROUP BY name, checksum
        """,
        con
    )

    table["any_signal_peptide"] = decide_any_signal(table)
    table["multiple_transmembrane"] = decide_is_multi_tm(table)
    table["single_transmembrane"] = decide_is_single_tm(
        table,
        tmhmm_first_60_threshold=tmhmm_first_60_threshold
    )

    table["is_secreted"] = decide_is_secreted(table)
    return table


def decide_any_signal(
    table: pd.DataFrame
) -> pd.Series:
    return (
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


def decide_is_multi_tm(table: pd.DataFrame) -> pd.Series:
    return (
        (table["tmhmm_tmcount"] > 1) |
        (table["phobius_tmcount"] > 1)
    ).astype(int)


def decide_is_single_tm(
    table: pd.DataFrame,
    tmhmm_first_60_threshold: float = 10,
) -> pd.Series:
    return (
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


def decide_is_secreted(
    table: pd.DataFrame
) -> pd.Series:
    return (
        table["any_signal_peptide"].astype(bool) & ~
        table["multiple_transmembrane"].astype(bool)
    ).astype(int)


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
    score = 1.0 * table["is_secreted"].fillna(0).astype(int) * secreted

    for k in ["signalp3_hmm", "signalp3_nn", "phobius_sp", "deepsig"]:
        col = table[k].fillna(0).astype(int)
        score += 1.0 * col * less_trustworthy_signal_prediction

    for k in ["signalp4", "signalp5", "signalp6", "targetp_secreted"]:
        col = table[k].fillna(0).astype(int)
        score += 1.0 * col * trustworthy_signal_prediction

    score += (
        1.0 *
        table["multiple_transmembrane"].fillna(0).astype(int) *
        multiple_transmembrane
    )

    score += (
        1.0 *
        table["single_transmembrane"].fillna(0).astype(int) *
        single_transmembrane
    )

    score += (
        table["deeploc_extracellular"].fillna(0.0).astype(float) *
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


def inner(
    con: sqlite3.Connection,
    cur: sqlite3.Cursor,
    args: argparse.Namespace
) -> None:
    if args.dbcan is not None:
        dbcan: Set[str] = {d.strip() for d in args.dbcan.readlines()}
    else:
        dbcan = set(get_interesting_dbcan_ids())

    if args.pfam is not None:
        pfam: Set[str] = {d.strip() for d in args.pfam.readlines()}
    else:
        pfam = set(get_interesting_pfam_ids())

    create_pfam_targets(con, cur, pfam)
    create_dbcan_targets(con, cur, dbcan)

    df = create_tables(
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
    df.round(3).to_csv(
        args.outfile,
        sep="\t",
        index=False,
        na_rep=".",
    )
    return


def runner(args: argparse.Namespace) -> None:
    try:
        con, cur = load_db(args.db, args.mem)
        inner(con, cur, args)
    except Exception as e:
        raise e
    finally:
        con.commit()
        con.close()
    return
