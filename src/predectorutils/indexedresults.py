#!/usr/bin/env python3

from typing import Tuple
from typing import TextIO
from typing import Any
from typing import Iterator

import json

import pandas as pd

from predectorutils.analyses import Analyses


class IndexedResults(object):
    """ There's potential to improve disk read efficiency here. """

    def __init__(
        self,
        handle: TextIO,
        index: pd.DataFrame
    ):

        self.index = index
        self.handle = handle
        return

    @classmethod
    def parse(cls, handle: TextIO) -> "IndexedResults":
        handle.seek(0, 2)  # Seeks to end of file
        EOF = handle.tell()

        handle.seek(0)

        rows: pd.Series = []

        while handle.tell() <= EOF:
            linestart = handle.tell()
            line = handle.readline().strip()
            lineend = handle.tell()

            if line == "":
                continue

            jline = json.loads(line)
            analysis = Analyses.from_string(jline["analysis"])
            an_object = analysis.get_analysis()
            software_version = jline.get("software_version", None)
            database_version = jline.get("database_version", None)
            checksum = jline["checksum"]

            if software_version is None:
                continue

            if (
                hasattr(an_object, "database")
                and (an_object.database is not None)
                and (database_version is None)
            ):
                # Cant precompute because don't know version
                continue

            rows.append(pd.Series({
                "analysis": jline["analysis"],
                "software_version": software_version,
                "database_version": database_version,
                "checksum": checksum,
                "start": linestart,
                "end": lineend,
            }))

        df = pd.DataFrame(rows)
        df.sort_values("start", ascending=True, inplace=True)
        return cls(handle, df)

    def __getitem__(self, key: Any) -> Iterator[str]:
        return self.index[key]

    def fetch_df(self, df: pd.DataFrame) -> Iterator[Tuple[pd.Series, str]]:
        for i, row in df.iterrows():
            self.handle.seek(row.start)
            yield row, self.handle.read(row.end - row.start)

    def fetch(self, key: Any) -> Iterator[Tuple[pd.Series, str]]:
        subdf = self[key]
        return self.fetch_df(subdf)
