#!/usr/bin/env python3

from typing import Callable
from typing import List
from typing import ClassVar
from typing import Dict
from typing import Union, Optional, Any
from typing import TextIO
from typing import Iterator


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
    types: ClassVar[List[Union[
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
        fields = tuple(
            type_(d.get(cname))
            for cname, type_
            in zip(cls.columns, cls.types)
        )
        return cls(*fields)

    def __repr__(self) -> str:
        inner = ", ".join([repr(getattr(self, k)) for k in self.columns])
        return f"{self.__class__.__name__}({inner})"

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator['Analysis']:
        raise NotImplementedError()


class GFFAble(Analysis):

    def to_gff(self) -> str:
        raise NotImplementedError()
