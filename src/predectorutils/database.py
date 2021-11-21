#!/usr/bin/env python3

from typing import NamedTuple, Tuple
from typing import Dict, Set
from typing import TextIO
from typing import Any, Optional
from typing import Iterator
from math import floor

import json
import sqlite3

from .analyses import Analysis, Analyses


def text_split(text: str, sep: str) -> str:
    return json.dumps(
        text.split(sep),
        separators=(',', ':')
    )


def load_db(
    path: str,
    mem: int = 1
) -> Tuple[sqlite3.Connection, sqlite3.Cursor]:
    sqlite3.register_converter("analyses", Analyses.from_bytes_)
    con = sqlite3.connect(path, detect_types=sqlite3.PARSE_DECLTYPES)
    con.row_factory = sqlite3.Row

    cur = con.cursor()
    # Allow it to use 1GB RAM for cache
    mem_gb: int = floor(1000000 * mem)
    cur.execute(f"PRAGMA cache_size = -{mem_gb}")
    cur.execute("PRAGMA journal_mode = WAL")
    cur.execute("PRAGMA locking_mode = EXCLUSIVE")
    cur.execute("PRAGMA synchronous = NORMAL")

    con.commit()
    return con, cur


class ResultRow(NamedTuple):

    analysis: Analyses
    software: str
    software_version: str
    database: Optional[str]
    database_version: Optional[str]
    pipeline_version: Optional[str]
    checksum: str
    data: str

    @classmethod
    def from_string(cls, s: str, replace_name: bool = False) -> "ResultRow":
        d = json.loads(s.strip())
        assert isinstance(d["analysis"], str), d
        assert isinstance(d["software"], str), d
        assert isinstance(d["software_version"], str), d

        database = d.get("database", None)
        database_version = d.get("database_version", None)
        pipeline_version = d.get("pipeline_version", None)
        assert isinstance(database, str) or database is None, d
        assert isinstance(database_version, str) or database_version is None, d
        assert isinstance(pipeline_version, str) or pipeline_version is None, d
        assert isinstance(d["checksum"], str), d

        an_enum = Analyses.from_string(d["analysis"])
        # This ensures the types are all correcy
        an_obj = (
            an_enum
            .get_analysis()
            .from_dict(d["data"])
        )

        if replace_name:
            an_obj.replace_name()

        data = an_obj.as_json_str()

        return cls(
            an_enum,
            d["software"],
            d["software_version"],
            database,
            database_version,
            pipeline_version,
            d["checksum"],
            data,
        )

    def replace_name(self):
        an = (
            self.analysis
            .get_analysis()
            .from_json_str(self.data)
            .replace_name()
        )

        return self.__class__(
            self.analysis,
            self.software,
            self.software_version,
            self.database,
            self.database_version,
            self.pipeline_version,
            self.checksum,
            an.as_json_str()
        )

    @classmethod
    def from_file(
        cls,
        handle: TextIO,
        replace_name: bool = False,
        drop_null_dbversion: bool = False,
        target_analyses: Optional[Set[Analyses]] = None
    ) -> Iterator["ResultRow"]:

        if drop_null_dbversion:
            requires_database = {
                a
                for a
                in Analyses
                if (a.get_analysis().database is not None)
            }

        for line in handle:
            sline = line.strip()
            if sline == "":
                continue
            record = cls.from_string(sline, replace_name=replace_name)

            if target_analyses is not None:
                if record.analysis not in target_analyses:
                    continue

            if drop_null_dbversion:
                if (
                    (record.database_version is None)
                    and (record.analysis in requires_database)
                ):
                    continue

            yield record
        return

    def as_dict(self) -> Dict[str, Any]:
        d = {
            "analysis": str(self.analysis),
            "software": self.software,
            "software_version": self.software_version,
            "checksum": self.checksum,
            "data": json.loads(self.data)
        }

        if self.database is not None:
            d["database"] = self.database

        if self.database_version is not None:
            d["database_version"] = self.database_version

        if self.pipeline_version is not None:
            d["pipeline_version"] = self.pipeline_version
        return d

    def as_str(self) -> str:
        return json.dumps(self.as_dict(), separators=(',', ':'))

    def as_analysis(self) -> "Analysis":
        an = (
            self.analysis
            .get_analysis()
            .from_json_str(self.data)
        )
        return an

    @classmethod
    def from_rowfactory(cls, row: sqlite3.Row) -> "ResultRow":
        return cls(
            analysis=row["analysis"],
            software=row["software"],
            software_version=row["software_version"],
            database=row["database"],
            database_version=(
                None
                if (row["database_version"] == '.')
                else row["database_version"]
            ),
            pipeline_version=row["pipeline_version"],
            checksum=row["checksum"],
            data=row["data"]
        )


class TargetRow(NamedTuple):

    analysis: Analyses
    software_version: str
    database_version: Optional[str]

    @classmethod
    def from_string(cls, s: str) -> "TargetRow":
        ss = s.strip().split("\t")

        if len(ss) == 2:
            return cls(Analyses.from_string(ss[0]), ss[1], None)
        elif len(ss) == 3:
            return cls(Analyses.from_string(ss[0]), ss[1], ss[2])
        else:
            raise ValueError("Target table is in improper format")

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["TargetRow"]:
        header = ["analysis", "software_version", "database_version"]
        for line in handle:
            sline = line.strip()
            if sline in ("\t".join(header), "\t".join(header[:2])):
                continue
            elif sline == "":
                continue

            yield cls.from_string(sline)
        return

    def as_dict(self) -> Dict[str, Optional[str]]:
        return {
            "analysis": str(self.analysis),
            "software_version": self.software_version,
            "database_version": self.database_version
        }

    @classmethod
    def from_rowfactory(cls, row: sqlite3.Row) -> "TargetRow":
        return cls(
            row["analysis"],
            row["software_version"],
            (
                None
                if (row["database_version"] == '.')
                else row["database_version"]
            )
        )


class DecoderRow(NamedTuple):

    encoded: str
    filename: str
    id: str
    checksum: str

    @classmethod
    def from_string(cls, s: str) -> "DecoderRow":
        e = s.strip().split("\t")
        return DecoderRow(e[0], e[1], e[2], e[3])

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["DecoderRow"]:
        header = ["encoded", "filename", "id", "checksum"]
        for line in handle:
            sline = line.strip()
            if sline == "\t".join(header):
                continue
            elif sline == "":
                continue

            yield cls.from_string(sline)
        return


class DecodedRow(NamedTuple):

    encoded: str
    filename: str
    id: str
    analysis: Analyses
    software: str
    software_version: str
    database: Optional[str]
    database_version: Optional[str]
    pipeline_version: Optional[str]
    checksum: str
    data: str

    def as_analysis(self) -> Analysis:
        return (
            self.analysis
            .get_analysis()
            .from_json_str(self.data)
        )

    def as_result_row(self) -> ResultRow:
        an = self.as_analysis()
        an.replace_name(self.id)
        return ResultRow(
            self.analysis,
            self.software,
            self.software_version,
            self.database,
            self.database_version,
            self.pipeline_version,
            self.checksum,
            an.as_json_str()
        )

    def as_result_string(self) -> str:
        """ This is basically a copy to avoid having to serialise data again"""
        an = self.as_analysis()
        an.replace_name(self.id)

        d = {
            "analysis": str(self.analysis),
            "software": self.software,
            "software_version": self.software_version,
            "checksum": self.checksum,
            "data": an.as_dict()
        }

        if self.database is not None:
            d["database"] = self.database

        if self.database_version is not None:
            d["database_version"] = self.database_version

        if self.pipeline_version is not None:
            d["pipeline_version"] = self.pipeline_version

        return json.dumps(d, separators=(',', ':'))

    @classmethod
    def from_rowfactory(cls, row: sqlite3.Row) -> "DecodedRow":
        return cls(
            encoded=row["encoded"],
            filename=row["filename"],
            id=row["id"],
            analysis=row["analysis"],
            software=row["software"],
            software_version=row["software_version"],
            database=row["database"],
            database_version=(
                None
                if (row["database_version"] == '.')
                else row["database_version"]
            ),
            pipeline_version=row["pipeline_version"],
            checksum=row["checksum"],
            data=row["data"]
        )


class ResultsTable(object):

    def __init__(self, con: sqlite3.Connection, cur: sqlite3.Cursor) -> None:
        self.con = con
        self.cur = cur
        return

    def create_tables(self):
        self.cur.execute(
            """
            CREATE TABLE IF NOT EXISTS results (
                analysis analyses NOT NULL,
                software text NOT NULL,
                software_version text NOT NULL,
                database text,
                database_version text NOT NULL,
                pipeline_version text,
                checksum text NOT NULL,
                data json NOT NULL,
                UNIQUE (
                    analysis,
                    software_version,
                    database_version,
                    checksum,
                    data
                )
            )
            """
        )

        self.con.commit()
        return

    def index_results(self):
        self.cur.execute(
            """
            CREATE INDEX IF NOT EXISTS analysis_version
            ON results (
                analysis,
                software_version,
                database_version,
                checksum
            )
            """
        )
        self.con.commit()
        return

    def drop_index(self):
        self.cur.execute("DROP INDEX IF EXISTS analysis_version")
        self.con.commit()
        return

    def insert_results(self, rows: Iterator[ResultRow]) -> None:
        self.cur.executemany(
            """
            INSERT INTO results
            VALUES (
                :analysis,
                :software,
                :software_version,
                :database,
                IFNULL(:database_version, '.'),
                :pipeline_version,
                :checksum,
                json(:data)
            )
            ON CONFLICT DO NOTHING
            """,
            rows
        )
        self.con.commit()
        return

    def exists_table(self, table: str) -> bool:
        query = (
            "SELECT 1 FROM sqlite_master "
            "WHERE type IN ('table', 'view') and name = ?"
        )
        main = self.cur.execute(query, (table,)).fetchone() is not None

        query = (
            "SELECT 1 FROM sqlite_temp_master "
            "WHERE type IN ('table', 'view') and name = ?"
        )
        temp = self.cur.execute(query, (table,)).fetchone() is not None
        return main or temp

    def insert_checksums(self, checksums: Set[str]) -> None:
        self.cur.execute("DROP TABLE IF EXISTS checksums")
        self.cur.execute(
            """
            CREATE TEMP TABLE checksums (checksum text NOT NULL UNIQUE)
            """
        )
        self.cur.executemany(
            "INSERT INTO checksums VALUES (?)",
            ((c,) for c in checksums)
        )

    def select_checksums(self) -> Iterator[ResultRow]:
        assert self.exists_table("checksums"), "no checksums table"

        result = self.cur.execute(
            """
            SELECT DISTINCT r.*
            FROM results r
            INNER JOIN checksums c
                ON r.checksum = c.checksum
            """
        )

        for r in result:
            yield ResultRow.from_rowfactory(r)

        self.con.commit()
        return

    def select_target(
        self,
        target: TargetRow,
        checksums: bool = False,
    ) -> Iterator[ResultRow]:

        if checksums:
            assert self.exists_table("checksums"), "no checksums table"

            result = self.cur.execute(
                """
                SELECT
                    analysis,
                    software,
                    software_version,
                    database,
                    (
                        CASE
                            WHEN database_version == '.'
                            THEN NULL
                            ELSE database_version
                        END
                    ) as database_version,
                    pipeline_version,
                    checksum,
                    data
                FROM results
                WHERE analysis = CAST(:analysis AS analyses)
                AND software_version = :software_version
                AND database_version = IFNULL(:database_version, '.')
                AND checksum IN checksums
                """,
                target
            )

        else:
            result = self.cur.execute(
                """
                SELECT
                    analysis,
                    software,
                    software_version,
                    database,
                    (
                        CASE
                            WHEN database_version == '.'
                            THEN NULL
                            ELSE database_version
                        END
                    ) as database_version,
                    pipeline_version,
                    checksum,
                    data
                FROM results
                WHERE analysis = CAST(:analysis AS analyses)
                AND software_version = :software_version
                AND database_version = IFNULL(:database_version, '.')
                """,
                target
            )

        for r in result:
            yield ResultRow.from_rowfactory(r)

        return

    def find_remaining(self, target: TargetRow) -> Iterator[str]:
        assert self.exists_table("checksums"), "no checksums table"

        result = self.cur.execute(
            """
            SELECT DISTINCT
                c.checksum AS checksum
            FROM checksums AS c
            WHERE c.checksum NOT IN (
                SELECT r.checksum
                FROM results AS r
                WHERE analysis = CAST(:analysis AS analyses)
                AND software_version = :software_version
                AND database_version = IFNULL(:database_version, '.')
            )
            """,
            target
        )

        for r in result:
            yield r["checksum"]
        return

    def select_all(self) -> Iterator[ResultRow]:
        result = self.cur.execute(
            """
            SELECT
                analysis,
                software,
                software_version,
                database,
                (
                    CASE
                        WHEN database_version == '.'
                        THEN NULL
                        ELSE database_version
                    END
                ) as database_version,
                pipeline_version,
                checksum,
                data
            FROM results
            """
        )
        for r in result:
            yield ResultRow.from_rowfactory(r)
        return

    def insert_decoder(self, rows: Iterator[DecoderRow]) -> None:
        self.cur.execute("DROP TABLE IF EXISTS decoder")
        self.cur.execute(
            """
            CREATE TEMP TABLE IF NOT EXISTS decoder (
                encoded text NOT NULL,
                filename text NOT NULL,
                id text NOT NULL,
                checksum text NOT NULL
            )
            """
        )

        self.cur.executemany(
            "INSERT INTO decoder VALUES (:encoded, :filename, :id, :checksum)",
            rows
        )

        self.cur.execute(
            """
            CREATE TEMP VIEW IF NOT EXISTS decoded
            AS
            SELECT
                d.encoded,
                d.filename,
                d.id,
                r.analysis,
                r.software,
                r.software_version,
                r.database,
                (
                    CASE
                        WHEN r.database_version == '.'
                        THEN NULL
                        ELSE r.database_version
                    END
                ) as database_version,
                r.pipeline_version,
                r.checksum,
                r.data
            FROM results r
            INNER JOIN decoder d
                ON r.checksum = d.checksum
            """
        )
        self.con.commit()
        return

    def decode(self) -> Iterator[Tuple[str, Iterator[ResultRow]]]:
        assert self.exists_table("decoded"), "table decoder doesn't exist"

        fnames = (
            self.cur.execute("SELECT DISTINCT filename FROM decoded")
            .fetchall()
        )

        for fname, in fnames:
            rows = self.cur.execute(
                "SELECT DISTINCT * FROM decoded WHERE filename = :filename",
                {"filename": fname}
            )

            gen = (
                (DecodedRow
                    .from_rowfactory(r)
                    .as_result_row())
                for r
                in rows
            )
            yield fname, gen
        return
