#!/usr/bin/env python3

from typing import Type
import enum

from predector.analyses.base import Analysis, GFFAble
from predector.analyses.apoplastp import ApoplastP
from predector.analyses.deepsig import DeepSig
from predector.analyses.effectorp import EffectorP1, EffectorP2
from predector.analyses.phobius import Phobius
from predector.analyses.signalp import (
    SignalP3NN, SignalP3HMM, SignalP4, SignalP5)
from predector.analyses.targetp import TargetPNonPlant, TargetPPlant
from predector.analyses.tmhmm import TMHMM
from predector.analyses.localizer import LOCALIZER
from predector.analyses.deeploc import DeepLoc
from predector.analyses.hmmer import DBCAN

__all__ = ["Analysis", "ApoplastP", "DeepSig", "EffectorP1", "EffectorP2",
           "Phobius", "SignalP3NN", "SignalP3HMM", "SignalP4", "SignalP5",
           "TargetPNonPlant", "TargetPPlant", "TMHMM", "LOCALIZER", "DeepLoc",
           "DBCAN", "GFFAble"]


class Analyses(enum.Enum):

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
    mmseqs = 17

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


NAME_TO_ANALYSIS = {
    Analyses.signalp3_nn: SignalP3NN,
    Analyses.signalp3_hmm: SignalP3HMM,
    Analyses.signalp4: SignalP4,
    Analyses.signalp5: SignalP5,
    Analyses.deepsig: DeepSig,
    Analyses.phobius: Phobius,
    Analyses.tmhmm: TMHMM,
    Analyses.targetp_plant: TargetPPlant,
    Analyses.targetp_non_plant: TargetPNonPlant,
    Analyses.effectorp1: EffectorP1,
    Analyses.effectorp2: EffectorP2,
    Analyses.apoplastp: ApoplastP,
    Analyses.localizer: LOCALIZER,
    Analyses.deeploc: DeepLoc,
    Analyses.dbcan: DBCAN,
}

#   Analyses.pfamscan: PfamScan,
#   Analyses.hmmer_domtbl: HMMERDomTab,
#   Analyses.mmseqs: MMSeqs,
