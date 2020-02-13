#!/usr/bin/env python3

from typing import Callable
from typing import List
from typing import ClassVar
from typing import Dict
from typing import Union, Optional, Any


from predector.analyses.apoplastp import ApoplastP
from predector.analyses.deepsig import DeepSig
from predector.analyses.effectorp import EffectorP1, EffectorP2
from predector.analyses.phobius import Phobius
from predector.analyses.signalp import (
    SignalP3NN, SignalP3HMM, SignalP4, SignalP5)
from predector.analyses.targetp import TargetPNonPlant, TargetPPlant
from predector.analyses.tmhmm import TMHMM

__all__ = ["Analysis", "ApoplastP", "DeepSig", "EffectorP1", "EffectorP2",
           "Phobius", "SignalP3NN", "SignalP3HMM", "SignalP4", "SignalP5",
           "TargetPNonPlant", "TargetPPlant", "TMHMM", "GFFAble"]


def int_or_none(i: Any) -> Optional[int]:
    if i is None:
        return None
    else:
        return int(i)


def str_or_none(i: Any) -> Optional[str]:
    if i is None:
        return None
    else:
        return str(i)


def float_or_none(i: Any) -> Optional[float]:
    if i is None:
        return None
    else:
        return float(i)


class Analysis(object):

    columns: ClassVar[List[str]] = []
    coltypes: ClassVar[List[Union[
        Callable[[Any], int],
        Callable[[Any], str],
        Callable[[Any], float],
        Callable[[Any], Optional[int]],
        Callable[[Any], Optional[str]],
        Callable[[Any], Optional[float]],
    ]]] = []

    def as_dict(self) -> Dict[str, Union[str, int, float, bool, None]]:
        return {k: getattr(self, k) for k in self.columns}

    @classmethod
    def from_dict(
        cls,
        d: Dict[str, Union[str, int, float, bool, None]]
    ) -> "Analysis":
        fields = (
            type_(d.get(cname))
            for cname, type_
            in zip(cls.columns, cls.coltypes)
        )
        return cls(*fields)

    def __repr__(self) -> str:
        inner = ", ".join([repr(getattr(self, k)) for k in self.columns])
        return f"{self.__class__.__name__}({inner})"


class GFFAble(Analysis):

    def to_gff(self) -> str:
        return "to do"