#!/usr/bin/env python3

from typing import Dict, List
from typing import NamedTuple, Tuple
from typing import TextIO
from typing import Any, Optional
from typing import Iterator

import json

from predectorutils.analyses import Analyses


class AnalysisTuple(NamedTuple):

    analysis: Analyses
    software_version: Optional[str]
    database_version: Optional[str]


#               (name, version) checksum: (start, end)
IndexType = Dict[AnalysisTuple, Dict[str, Tuple[int, int]]]


class IndexedResults(object):
    """ There's potential to improve disk read efficiency here. """

    def __init__(
        self,
        handle: TextIO,
        index: IndexType
    ):

        self.index = dict(index)
        self.handle = handle
        return

    def analyses(self) -> List[AnalysisTuple]:
        return list(self.index.keys())

    @classmethod
    def parse(cls, handle: TextIO) -> "IndexedResults":
        handle.seek(0, 2)  # Seeks to end of file
        EOF = handle.tell()

        handle.seek(0)

        idx: IndexType = dict()

        while handle.tell() <= EOF:
            linestart = handle.tell()
            line = handle.readline().strip()
            lineend = handle.tell()

            if line == "":
                continue

            jline = json.loads(line)
            analysis = jline["analysis"]
            software_version = jline.get("software_version", None)
            database_version = jline.get("database_version", None)

            at = AnalysisTuple(
                Analyses.from_string(analysis),
                software_version,
                database_version
            )

            checksum = jline["checksum"]

            if at in idx:
                idx[at][checksum] = (linestart, lineend)
            else:
                idx[at] = {checksum: (linestart, lineend)}

        return cls(handle, idx)

    def __getitem__(self, key: Any) -> Iterator[str]:  # noqa
        from collections.abc import Iterable

        if isinstance(key, AnalysisTuple):
            idx1 = self.index.get(key, None)
            if idx1 is None:
                return

            for chk, (start, end) in idx1.items():
                self.handle.seek(start)
                yield self.handle.read(end - start)

        elif isinstance(key, Iterable):
            for an in key:
                assert isinstance(an, AnalysisTuple), an
                idx1 = self.index.get(an, None)
                if idx1 is None:
                    return

                for chk, (start, end) in idx1.items():
                    self.handle.seek(start)
                    yield self.handle.read(end - start)

        elif isinstance(key, tuple):
            assert len(key) == 2, key
            if isinstance(key[0], AnalysisTuple):
                ans = [key[0]]
            elif isinstance(key[0], Iterable):
                ans = list(key[0])
                for k in ans:
                    assert isinstance(k, AnalysisTuple), k
            elif key[0] is None:
                ans = list(self.index.keys())
            else:
                raise KeyError("The first element of the tuple must be "
                               "an AnalysisTuple, list of AnalysisTuples or "
                               "None.")

            if isinstance(key[1], str):
                chks = [key[1]]
            elif isinstance(key[1], Iterable):
                chks = list(key[1])
                for c in chks:
                    assert isinstance(c, str), c
            else:
                raise KeyError("The second element of the tuple must be "
                               "a str or iterable of str")

            for an in ans:
                idx1 = self.index.get(an, None)
                if idx1 is None:
                    continue

                for chk in chks:
                    idx2 = idx1.get(chk, None)
                    del idx2
                    if idx2 is None:
                        continue

                    start, end = idx2

                    self.handle.seek(start)
                    yield self.handle.read(end - start)

        return
