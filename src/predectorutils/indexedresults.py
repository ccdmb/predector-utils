#!/usr/bin/env python3

from typing import NamedTuple, Tuple
from typing import Dict, Set
from typing import TextIO
from typing import Optional
from typing import Iterator

import json
import sqlite3

from predectorutils.analyses import Analysis, Analyses


class ResultRow(NamedTuple):

    analysis: str
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

        if replace_name:
            an_cls = Analyses.from_string(d["analysis"]).get_analysis()
            an = an_cls.from_dict(d["data"])
            setattr(an, an.name_column, "dummy")
            d["data"] = an.as_dict()

        return cls(
            d["analysis"],
            d["software"],
            d["software_version"],
            database,
            database_version,
            pipeline_version,
            d["checksum"],
            json.dumps(d["data"]),
        )

    def replace_name(self):
        d = json.loads(self.data)
        an_cls = Analyses.from_string(self.analysis).get_analysis()
        an = an_cls.from_dict(d)
        setattr(an, an.name_column, "dummy")
        return self.__class__(
            self.analysis,
            self.software,
            self.software_version,
            self.database,
            self.database_version,
            self.pipeline_version,
            self.checksum,
            json.dumps(an.as_dict())
        )

    @classmethod
    def from_file(
        cls,
        handle: TextIO,
        replace_name: bool = False
    ) -> Iterator["ResultRow"]:
        for line in handle:
            sline = line.strip()
            if sline == "":
                continue
            yield cls.from_string(sline, replace_name=replace_name)
        return

    def as_str(self) -> str:
        d = {
            "analysis": self.analysis,
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

        return json.dumps(d)

    def as_analysis(self) -> "Analysis":
        cls = Analyses.from_string(self.analysis).get_analysis()
        d = json.loads(self.data)
        return cls.from_dict(d)


class TargetRow(NamedTuple):

    analysis: str
    software_version: str
    database_version: Optional[str]

    @classmethod
    def from_string(cls, s: str) -> "TargetRow":
        ss = s.strip().split("\t")

        if len(ss) == 2:
            return cls(ss[0], ss[1], None)
        elif len(ss) == 3:
            return cls(ss[0], ss[1], ss[2])
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
            "analysis": self.analysis,
            "software_version": self.software_version,
            "database_version": self.database_version
        }


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
    analysis: str
    software: str
    software_version: str
    database: Optional[str]
    database_version: Optional[str]
    pipeline_version: Optional[str]
    checksum: str
    data: str

    def as_result_row(self) -> ResultRow:
        cls = Analyses.from_string(self.analysis).get_analysis()
        d = json.loads(self.data)
        an = cls.from_dict(d)
        setattr(an, an.name_column, self.id)
        return ResultRow(
            self.analysis,
            self.software,
            self.software_version,
            self.database,
            self.database_version,
            self.pipeline_version,
            self.checksum,
            json.dumps(an.as_dict())
        )

    def as_result_string(self) -> str:
        """ This is basically a copy to avoid having to serialise data again"""
        cls = Analyses.from_string(self.analysis).get_analysis()
        d = json.loads(self.data)
        an = cls.from_dict(d)
        setattr(an, an.name_column, self.id)

        d = {
            "analysis": self.analysis,
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

        return json.dumps(d)


class ResultsTable(object):

    def __init__(self, con: sqlite3.Connection, cur: sqlite3.Cursor) -> None:
        self.con = con
        self.cur = cur
        return

    def create_tables(self):
        self.cur.execute("""
        CREATE TABLE IF NOT EXISTS results (
            analysis text NOT NULL,
            software text NOT NULL,
            software_version text NOT NULL,
            database text,
            database_version text,
            pipeline_version text,
            checksum text NOT NULL,
            data json NOT NULL
        )
        """)

        self.cur.execute("""
        CREATE INDEX IF NOT EXISTS analysis_version
        ON results (analysis, software_version, database_version, checksum);
        """)
        self.con.commit()
        return

    def deduplicate_table(self):

        multiples_ok = ", ".join([
            "'" + str(a) + "'"
            for a
            in Analyses
            if a.multiple_ok()
        ])

        requires_database = ", ".join({
            "'" + str(a) + "'"
            for a
            in Analyses
            if (a.get_analysis().database is not None)
        })

        self.cur.execute(f"""
        CREATE TEMP TABLE IF NOT EXISTS results_deduplicated
        AS
        SELECT DISTINCT *
        FROM results
        WHERE (analysis IN ({multiples_ok}))
        UNION
        SELECT *
        FROM results
        WHERE (analysis NOT IN ({multiples_ok}))
        GROUP BY analysis, software_version, database_version, checksum
        HAVING ROWID=MIN(ROWID)
        """)

        self.cur.execute("DROP TABLE results")

        self.cur.execute(f"""
        CREATE TABLE results
        AS
        SELECT *
        FROM results_deduplicated
        WHERE NOT (
            (analysis IN ({requires_database}))
            AND (database_version IS NULL)
        )
        """)

        self.cur.execute("DROP TABLE results_deduplicated")
        self.con.commit()
        return

    def insert_results(self, rows: Iterator[ResultRow]) -> None:
        buf = []

        for row in rows:
            buf.append(row)

            if len(buf) > 50000:
                self.cur.executemany(
                    (
                        "INSERT INTO results "
                        "VALUES (?, ?, ?, ?, ?, ?, ?, json(?))"
                    ),
                    buf
                )
                self.con.commit()
                buf = []

        if len(buf) > 0:
            self.cur.executemany(
                (
                    "INSERT INTO results "
                    "VALUES (?, ?, ?, ?, ?, ?, ?, json(?))"
                ),
                buf
            )

        self.con.commit()
        return

    def insert_targets(self, rows: Iterator[TargetRow]) -> None:
        self.cur.execute("DROP TABLE IF EXISTS targets")
        self.cur.execute("""
        CREATE TEMP TABLE IF NOT EXISTS targets (
            analysis text NOT NULL,
            software_version text NOT NULL,
            database_version text
        )
        """)

        self.cur.executemany("INSERT INTO targets VALUES (?, ?, ?)", rows)

        self.cur.execute("DROP TABLE IF EXISTS subset")
        self.cur.execute("""
        CREATE TEMP TABLE IF NOT EXISTS subset
        AS
        SELECT DISTINCT r.*
        FROM results r
        INNER JOIN targets t
            ON r.analysis = t.analysis
            AND r.software_version = t.software_version
            AND (r.database_version = t.database_version
                OR (r.database_version IS NULL AND t.database_version IS NULL))
        """)
        self.con.commit()
        return

    def insert_decoder(self, rows: Iterator[DecoderRow]) -> None:
        self.cur.execute("DROP TABLE IF EXISTS decoder")
        self.cur.execute("""
        CREATE TEMP TABLE IF NOT EXISTS decoder (
            encoded text NOT NULL,
            filename text NOT NULL,
            id text NOT NULL,
            checksum text NOT NULL
        )
        """)

        self.cur.executemany("INSERT INTO decoder VALUES (?, ?, ?, ?)", rows)

        self.cur.execute("""
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
            r.database_version,
            r.pipeline_version,
            r.checksum,
            r.data
        FROM results r
        INNER JOIN decoder d
            ON r.checksum = d.checksum
        """)
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

    def select_checksums(
        self,
        chk: Set[str],
        table: str
    ) -> Iterator[ResultRow]:

        assert self.exists_table(table), f"table {table} does not exist"

        self.cur.execute("DROP TABLE IF EXISTS checksums")

        self.cur.execute("""
        CREATE TEMP TABLE checksums (checksum text NOT NULL)
        """)
        self.cur.executemany(
            "INSERT INTO checksums VALUES (?)",
            [(c,) for c in chk]
        )

        result = self.cur.execute(f"""
        SELECT DISTINCT r.*
        FROM {table} r
        INNER JOIN checksums c
            ON r.checksum = c.checksum
        """)

        for r in result:
            yield ResultRow(*r)

        self.con.commit()
        return

    def select_target(
        self,
        target: TargetRow,
        table: str = "results"
    ) -> Iterator[ResultRow]:
        assert self.exists_table(table), f"table {table} does not exist"

        result = self.cur.execute(
            "SELECT * "
            f"FROM {table} "
            "WHERE analysis = :analysis "
            "AND software_version = :software_version "
            "AND ("
            "  (database_version = :database_version) "
            "  OR (database_version IS NULL AND :database_version IS NULL)"
            ")",
            target.as_dict()
        )
        for r in result:
            yield ResultRow(*r)
        return

    def select_all(self, table: str) -> Iterator[ResultRow]:
        assert self.exists_table(table), f"table {table} does not exist"

        result = self.cur.execute(f"SELECT * FROM {table}")
        for r in result:
            yield ResultRow(*r)
        return

    def find_remaining(
        self,
        checksums: Set[str],
        table: str = "subset"
    ) -> Iterator[Tuple[TargetRow, Set[str]]]:
        assert self.exists_table(table), f"table {table} does not exist"

        self.cur.execute("DROP TABLE IF EXISTS checksums")

        self.cur.execute("""
        CREATE TEMP TABLE checksums (checksum text NOT NULL)
        """)
        self.cur.executemany(
            "INSERT INTO checksums VALUES (?)",
            [(c,) for c in checksums]
        )

        self.cur.execute("DROP TABLE IF EXISTS all_tasks")
        self.cur.execute("""
        CREATE TEMP TABLE IF NOT EXISTS all_tasks
        AS
        SELECT t.analysis, t.software_version, t.database_version, c.checksum
        FROM targets t
        CROSS JOIN checksums c
        """)

        self.cur.execute("DROP TABLE IF EXISTS remaining")
        self.cur.execute("""
        CREATE TEMP TABLE IF NOT EXISTS remaining
        AS
        SELECT a.analysis, a.software_version, a.database_version, a.checksum
        FROM all_tasks a
        EXCEPT
        SELECT t.analysis, t.software_version, t.database_version, t.checksum
        from subset t
        """)

        self.con.commit()

        targets = self.con.execute("SELECT * FROM targets").fetchall()

        for target in targets:
            chks = self.cur.execute(
                """
                SELECT DISTINCT checksum FROM remaining
                WHERE analysis = :analysis
                AND software_version = :software_version
                AND (
                    (database_version = :database_version)
                    OR (database_version IS NULL and :database_version IS NULL)
                )
                """,
                dict(zip(
                    ["analysis", "software_version", "database_version"],
                    target
                ))
            ).fetchall()
            yield TargetRow(*target), set(map(lambda t: t[0], chks))
        return

    def decode(self) -> Iterator[Tuple[str, Iterator[str]]]:
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

            gen = (DecodedRow(*r).as_result_string() for r in rows)
            yield fname, gen
        return
