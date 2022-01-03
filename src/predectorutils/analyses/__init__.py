#!/usr/bin/env python3

from typing import Dict, Type
import enum
import sqlite3


from predectorutils.analyses.base import Analysis, GFFAble
from predectorutils.analyses.apoplastp import ApoplastP
from predectorutils.analyses.deepsig import DeepSig
from predectorutils.analyses.effectorp import (
    EffectorP1, EffectorP2, EffectorP3, EffectorP3Fungal
)
from predectorutils.analyses.phobius import Phobius
from predectorutils.analyses.signalp import (
    SignalP3NN,
    SignalP3HMM,
    SignalP4,
    SignalP5,
    SignalP6,
)
from predectorutils.analyses.targetp import TargetPNonPlant, TargetPPlant
from predectorutils.analyses.tmhmm import TMHMM
from predectorutils.analyses.localizer import LOCALIZER
from predectorutils.analyses.deeploc import DeepLoc
from predectorutils.analyses.hmmer import DomTbl, DBCAN, EffectorDB
from predectorutils.analyses.pfamscan import PfamScan
from predectorutils.analyses.pepstats import PepStats
from predectorutils.analyses.mmseqs import MMSeqs, PHIBase
from predectorutils.analyses.hhr import HHRAlignment  # noqa
from predectorutils.analyses.regex import (
        RegexAnalysis,
        Kex2SiteAnalysis,
        RXLRLikeAnalysis
)

from predectorutils.analyses.deepredeff import (
    DeepredeffFungi,
    DeepredeffOomycete,
    DeepredeffBacteria
)


__all__ = [
    "Analysis", "ApoplastP", "DeepSig",
    "EffectorP1", "EffectorP2", "EffectorP3", "EffectorP3Fungal"
    "Phobius",
    "SignalP3NN", "SignalP3HMM", "SignalP4", "SignalP5", "SignalP6",
    "TargetPNonPlant", "TargetPPlant", "TMHMM", "LOCALIZER", "DeepLoc",
    "DomTbl", "DBCAN", "GFFAble", "PfamScan", "PepStats", "MMSeqs",
    "PHIBase", "HHRAligmment", "EffectorDB", "RegexAnalysis",
    "DeepredeffFungi", "DeepredeffOomycete", "DeepredeffBacteria",
    "Kex2SiteAnalysis", "RXLRLikeAnalysis"
]


class Analyses(enum.IntEnum):

    signalp3_nn = 1
    signalp3_hmm = 2
    signalp4 = 3
    signalp5 = 4
    deepsig = 5
    phobius = 6
    tmhmm = 7
    deeploc = 8
    targetp_plant = 9
    targetp_non_plant = 10
    effectorp1 = 11
    effectorp2 = 12
    apoplastp = 13
    localizer = 14
    pfamscan = 15
    dbcan = 16
    phibase = 17
    pepstats = 18
    effectordb = 19
    signalp6 = 20
    effectorp3 = 21
    effectorp3_fungal = 22
    regexanalysis = 23
    deepredeff_fungi = 24
    deepredeff_oomycete = 25
    deepredeff_bacteria = 26
    kex2_cutsite = 27
    rxlr_like_motif = 28

    def __str__(self) -> str:
        return self.name

    @classmethod
    def from_string(cls, s: str) -> "Analyses":
        try:
            return cls[s]
        except KeyError:
            raise ValueError(f"{s} is not a valid result type to parse.")

    def get_analysis(self) -> Type[Analysis]:
        return NAME_TO_ANALYSIS[self]

    def multiple_ok(self) -> bool:
        return MULTIPLE_ACCEPTABLE[self]

    def needs_database(self) -> bool:
        return NAME_TO_ANALYSIS[self].database is not None

    def name_column(self) -> str:
        return NAME_TO_ANALYSIS[self].name_column

    def __conform__(self, protocol):
        if protocol is sqlite3.PrepareProtocol:
            return self.value

    # int defines a from bytes method, hence underscore
    @classmethod
    def from_bytes_(cls, b: bytes) -> "Analyses":
        return cls(int(b))


NAME_TO_ANALYSIS: Dict[Analyses, Type[Analysis]] = {
    Analyses.signalp3_nn: SignalP3NN,
    Analyses.signalp3_hmm: SignalP3HMM,
    Analyses.signalp4: SignalP4,
    Analyses.signalp5: SignalP5,
    Analyses.signalp6: SignalP6,
    Analyses.deepsig: DeepSig,
    Analyses.phobius: Phobius,
    Analyses.tmhmm: TMHMM,
    Analyses.targetp_plant: TargetPPlant,
    Analyses.targetp_non_plant: TargetPNonPlant,
    Analyses.effectorp1: EffectorP1,
    Analyses.effectorp2: EffectorP2,
    Analyses.effectorp3: EffectorP3,
    Analyses.effectorp3_fungal: EffectorP3Fungal,
    Analyses.apoplastp: ApoplastP,
    Analyses.localizer: LOCALIZER,
    Analyses.deeploc: DeepLoc,
    Analyses.dbcan: DBCAN,
    Analyses.pfamscan: PfamScan,
    Analyses.pepstats: PepStats,
    Analyses.phibase: PHIBase,
    Analyses.effectordb: EffectorDB,
    Analyses.regexanalysis: RegexAnalysis,
    Analyses.deepredeff_fungi: DeepredeffFungi,
    Analyses.deepredeff_oomycete: DeepredeffOomycete,
    Analyses.deepredeff_bacteria: DeepredeffBacteria,
    Analyses.kex2_cutsite: Kex2SiteAnalysis,
    Analyses.rxlr_like_motif: RXLRLikeAnalysis,
}


MULTIPLE_ACCEPTABLE = {
    Analyses.signalp3_nn: False,
    Analyses.signalp3_hmm: False,
    Analyses.signalp4: False,
    Analyses.signalp5: False,
    Analyses.signalp6: False,
    Analyses.deepsig: False,
    Analyses.phobius: False,
    Analyses.tmhmm: False,
    Analyses.targetp_plant: False,
    Analyses.targetp_non_plant: False,
    Analyses.effectorp1: False,
    Analyses.effectorp2: False,
    Analyses.effectorp3: False,
    Analyses.effectorp3_fungal: False,
    Analyses.apoplastp: False,
    Analyses.localizer: False,
    Analyses.deeploc: False,
    Analyses.dbcan: True,
    Analyses.pfamscan: True,
    Analyses.pepstats: False,
    Analyses.phibase: True,
    Analyses.effectordb: True,
    Analyses.regexanalysis: True,
    Analyses.deepredeff_fungi: False,
    Analyses.deepredeff_oomycete: False,
    Analyses.deepredeff_bacteria: False,
    Analyses.kex2_cutsite: True,
    Analyses.rxlr_like_motif: True,
}
