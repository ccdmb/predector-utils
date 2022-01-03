#!/usr/bin/env python3

from typing import NamedTuple, Tuple
from typing import Dict, Set
from typing import TextIO
from typing import Any, Optional
from typing import Iterator, Iterable
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

    checksum: str
    name: Optional[str]
    filename: Optional[str]
    analysis: Analyses
    software: str
    software_version: str
    database: Optional[str]
    database_version: Optional[str]
    pipeline_version: Optional[str]
    data: str

    @classmethod
    def from_string(cls, s: str, drop_name: bool = False) -> "ResultRow":
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
        # This ensures the types are all correct
        an_obj = (
            an_enum
            .get_analysis()
            .from_dict(d["data"])
        )
        data = an_obj.as_dict()

        if drop_name:
            if an_obj.name_column in data:
                del data[an_obj.name_column]
            name = None
        else:
            name = data.get(an_obj.name_column, None)

        if name is not None:
            assert isinstance(name, str)

        if "filename" in data:
            assert isinstance(data["filename"], str), data
            filename: Optional[str] = data["filename"]
        else:
            filename = None

        return cls(
            d["checksum"],
            name,
            filename,
            an_enum,
            d["software"],
            d["software_version"],
            database,
            database_version,
            pipeline_version,
            json.dumps(data, separators=(',', ':'))
        )

    def drop_name(self):
        return self.__class__(
            self.checksum,
            None,
            None,
            self.analysis,
            self.software,
            self.software_version,
            self.database,
            self.database_version,
            self.pipeline_version,
            self.data
        )

    @classmethod
    def from_file(
        cls,
        handle: TextIO,
        drop_name: bool = False,
        drop_null_dbversion: bool = False,
        target_analyses: Optional[Set[Analyses]] = None
    ) -> Iterator["ResultRow"]:

        for line in handle:
            sline = line.strip()
            if sline == "":
                continue
            record = cls.from_string(sline, drop_name=drop_name)

            if target_analyses is not None:
                if record.analysis not in target_analyses:
                    continue

            if drop_null_dbversion:
                if (
                    (record.database_version is None)
                    and record.analysis.needs_database()
                ):
                    continue

            yield record
        return

    def as_dict(self) -> Dict[str, Any]:
        from .higher import or_else
        d = {
            "analysis": str(self.analysis),
            "software": self.software,
            "software_version": self.software_version,
            "checksum": self.checksum,
            "data": json.loads(self.data),
        }

        d["data"][self.analysis.name_column()] = or_else(".", self.name)

        if self.database is not None:
            d["database"] = self.database

        if self.filename is not None:
            d["filename"] = self.filename,

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
            .from_dict(self.as_dict()["data"])
        )
        return an

    @classmethod
    def from_rowfactory(cls, row: sqlite3.Row) -> "ResultRow":
        drow = dict(row)
        return cls(
            checksum=drow["checksum"],
            name=drow.get("name", None),
            filename=(
                None
                if (drow.get("filename", "") == '')
                else drow.get("filename", None)
            ),
            analysis=Analyses(int(drow["analysis"])),
            software=drow["software"],
            software_version=drow["software_version"],
            database=(
                None
                if (drow.get("database", "") == '')
                else drow.get("database", None)
            ),
            database_version=(
                None
                if (drow.get("database_version", "") == '')
                else drow.get("database_version", None)
            ),
            pipeline_version=drow["pipeline_version"],
            data=drow["data"]
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
            Analyses(int(row["analysis"])),
            row["software_version"],
            (
                None
                if (row["database_version"] == '')
                else row["database_version"]
            )
        )


class DecoderRow(NamedTuple):

    checksum: str
    filename: Optional[str]
    name: str

    @classmethod
    def from_string(cls, s: str) -> "DecoderRow":
        e = s.strip().split("\t")
        if e[1] in (".", "", None):
            fname = None
        else:
            fname = e[1]
        return DecoderRow(e[3], fname, e[2])

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["DecoderRow"]:
        header = "\t".join(["encoded", "filename", "id", "checksum"])
        for line in handle:
            sline = line.strip()
            if sline == header:
                continue
            elif sline == "":
                continue

            yield cls.from_string(sline)
        return


class ResultsTable(object):

    def __init__(self, con: sqlite3.Connection, cur: sqlite3.Cursor) -> None:
        self.con = con
        self.cur = cur
        return

    def create_result_tables(self):
        from predectorutils import analyses
        for an in analyses.Analyses:
            ans = str(an)

            if an.needs_database():
                db_cols = (
                    """
                    database text NOT NULL,
                    database_version text NOT NULL,
                    """
                )
                db_index = "database_version,"
            else:
                db_cols = ""
                db_index = ""

            if an.multiple_ok():
                unique_constraint = "checksum,data"
            else:
                unique_constraint = "checksum"

            self.cur.execute(
                f"""
                CREATE TABLE IF NOT EXISTS results_{ans} (
                    software text NOT NULL,
                    software_version text NOT NULL,
                    {db_cols}
                    pipeline_version text,
                    checksum text NOT NULL,
                    data json NOT NULL,
                    UNIQUE (
                        software_version,
                        {db_index}
                        {unique_constraint}
                    )
                )
                """
            )

        self.con.commit()
        return

    def create_result_index(self):
        from predectorutils import analyses
        for an in analyses.Analyses:
            ans = str(an)

            if an.needs_database():
                db_index = "database_version,"
            else:
                db_index = ""

            self.cur.execute(
                f"""
                CREATE INDEX IF NOT EXISTS results_index_{ans}
                ON results_{ans} (
                    software_version,
                    {db_index}
                    checksum
                )
                """
            )

        self.con.commit()
        return

    def drop_result_index(self):
        from predectorutils import analyses
        for an in analyses.Analyses:
            ans = str(an)
            self.cur.execute(f"DROP INDEX IF EXISTS results_index_{ans}")

        self.con.commit()
        return

    def create_decoder_table(self) -> None:
        self.cur.execute("DROP TABLE IF EXISTS decoder")
        self.cur.execute(
            """
            CREATE TABLE IF NOT EXISTS decoder (
                checksum text NOT NULL,
                filename text NOT NULL,
                name text NOT NULL,
                UNIQUE (
                    checksum,
                    filename,
                    name
                )
            )
            """
        )
        self.con.commit()
        return

    def drop_decoder_table(self):
        self.cur.execute("DROP TABLE IF EXISTS decoder")
        self.con.commit()
        return

    def create_decoder_index(self):
        self.cur.execute(
            """
            CREATE INDEX IF NOT EXISTS decoder_index
            ON decoder (checksum)
            """
        )
        self.con.commit()
        return

    def drop_decoder_index(self):
        self.cur.execute("DROP INDEX IF EXISTS decoder_index")
        self.con.commit()
        return

    def _insert_results(self, an: Analyses, rows: Iterable[ResultRow]) -> None:
        ans = str(an)
        keys = ["software", "software_version",
                "pipeline_version", "checksum", "data"]

        if an.needs_database():
            dcols = (
                """
                :database,
                :database_version,
                """
            )
            keys.extend(["database", "database_version"])
        else:
            dcols = ""

        self.cur.executemany(
            f"""
            INSERT INTO results_{ans}
            VALUES (
                :software,
                :software_version,
                {dcols}
                :pipeline_version,
                :checksum,
                json(:data)
            )
            ON CONFLICT DO NOTHING
            """,
            map(lambda r: {k: getattr(r, k) for k in keys}, rows)
        )
        self.con.commit()
        return

    def insert_results(
        self,
        rows: Iterable[ResultRow],
        insert_names: bool = False
    ) -> None:
        from collections import defaultdict
        cache = defaultdict(list)

        decoder_rows = []

        if insert_names:
            assert self.exists_table("decoder"), "no decoder table"

        i = 0
        for row in rows:
            an = row.analysis
            cache[an].append(row)

            if insert_names:

                if row.name is not None:
                    decoder_rows.append(DecoderRow(
                        row.checksum,
                        row.filename,
                        row.name
                    ))

            if i > 100000:
                for an, an_rows in cache.items():
                    if len(an_rows) == 0:
                        continue
                    self._insert_results(an, an_rows)
                    cache[an] = []

                if len(decoder_rows) > 0:
                    self.insert_decoder(decoder_rows)
                    decoder_rows = []
                i = 0

        for an, an_rows in cache.items():
            if len(an_rows) == 0:
                continue

            self._insert_results(an, an_rows)

        if len(decoder_rows) > 0:
            self.insert_decoder(decoder_rows)
        return

    def insert_decoder(self, rows: Iterable[DecoderRow]) -> None:
        self.cur.executemany(
            "INSERT INTO decoder VALUES "
            "(:checksum, IFNULL(:filename, ''), :name) "
            "ON CONFLICT DO NOTHING",
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
        from . import analyses
        assert self.exists_table("checksums"), "no checksums table"

        query = []
        for an in analyses.Analyses:
            ans = str(an)
            anv = an.value

            if not self.exists_table(f"results_{ans}"):
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
            SELECT DISTINCT
                CAST({anv} AS analyses) as analysis,
                software,
                software_version,
                {db_cols}
                pipeline_version,
                r.checksum as checksum,
                data
            FROM result_{ans} r
            INNER JOIN checksums c
                ON r.checksum = c.checksum
            """
            query.append(q)

        qstring = "\nUNION\n".join(query)
        result = self.cur.execute(qstring)

        for r in result:
            yield ResultRow.from_rowfactory(r)

        return

    def select_target(
        self,
        target: TargetRow,
        checksums: bool = False,
    ) -> Iterator[ResultRow]:
        an = target.analysis
        anstr = str(an)
        table_name = f"results_{anstr}"

        if not self.exists_table(f"results_{anstr}"):
            return

        if checksums:
            assert self.exists_table("checksums"), "no checksums table"
            chk_where = "AND checksum IN checksums"
        else:
            chk_where = ""

        target_dict = {
            "analysis": target.analysis,
            "software_version": target.software_version
        }
        if an.needs_database():
            db_cols = (
                """
                database,
                database_version,
                """
            )
            db_where = "AND database_version = IFNULL(:database_version, '')"
            target_dict["database_version"] = target.database_version
        else:
            db_cols = ""
            db_where = ""

        result = self.cur.execute(
            f"""
            SELECT
                CAST(:analysis AS analyses) as analysis,
                software,
                software_version,
                {db_cols}
                pipeline_version,
                checksum,
                data
            FROM {table_name}
            WHERE software_version = :software_version
            {db_where}
            {chk_where}
            """,
            target_dict
        )

        for r in result:
            yield ResultRow.from_rowfactory(r)

        return

    def find_remaining(self, target: TargetRow) -> Iterator[str]:
        assert self.exists_table("checksums"), "no checksums table"
        an = target.analysis
        anstr = str(an)
        table_name = f"results_{anstr}"

        if not self.exists_table(table_name):
            return

        if an.needs_database():
            db_where = "AND database_version = IFNULL(:database_version, '')"
        else:
            db_where = ""

        result = self.cur.execute(
            f"""
            SELECT DISTINCT
                c.checksum AS checksum
            FROM checksums AS c
            WHERE c.checksum NOT IN (
                SELECT r.checksum
                FROM {table_name} AS r
                WHERE software_version = :software_version
                {db_where}
            )
            """,
            target
        )

        for r in result:
            yield r["checksum"]
        return

    def select_all(self, checksums: bool = False) -> Iterator[ResultRow]:
        from . import analyses

        if checksums:
            assert self.exists_table("checksums"), "no checksums table"
            chk_join = "INNER JOIN checksums c ON r.checksum = c.checksum"
        else:
            chk_join = ""

        query = []
        for an in analyses.Analyses:
            ans = str(an)
            anv = an.value

            if not self.exists_table(f"results_{ans}"):
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
                data
            FROM results_{ans} r
            {chk_join}
            """
            query.append(q)

        qstring = "\nUNION\n".join(query)
        result = self.cur.execute(qstring)

        for r in result:
            yield ResultRow.from_rowfactory(r)

        return

    def decode(self) -> Iterator[Tuple[str, Iterator[ResultRow]]]:
        from . import analyses
        assert self.exists_table("decoder"), "table decoder doesn't exist"

        fnames = (
            self.cur.execute("SELECT DISTINCT filename FROM decoder")
            .fetchall()
        )

        query_template = """
        SELECT DISTINCT
            CAST({anv} AS analyses) as analysis,
            r.checksum as checksum,
            d.name as name,
            d.filename as filename,
            software,
            software_version,
            {db_cols}
            pipeline_version,
            data
        FROM results_{ans} r
        INNER JOIN fname_decoder d
            ON r.checksum = d.checksum
        """

        for fname, in fnames:
            self.cur.execute(
                """
                CREATE TEMP TABLE fname_decoder
                AS
                SELECT DISTINCT * FROM decoder WHERE filename = :filename
                """,
                {"filename": fname}
            )

            query = []
            for an in analyses.Analyses:
                ans = str(an)
                anv = an.value

                if not self.exists_table(f"results_{ans}"):
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

                q = query_template.format(
                    ans=ans,
                    anv=anv,
                    db_cols=db_cols
                )

                query.append(q)

            qstring = "\nUNION\n".join(query)
            result = self.cur.execute(qstring)

            gen = (
                ResultRow.from_rowfactory(r)
                for r
                in result
            )
            yield fname, gen
            self.cur.execute("DROP TABLE IF EXISTS fname_decoder")
        return

    def fetch_targets(
        self,
    ) -> Iterator[TargetRow]:
        from . import analyses

        query = []
        for an in analyses.Analyses:
            ans = str(an)
            anv = an.value

            if not self.exists_table(f"results_{ans}"):
                continue

            if an.needs_database():
                db_cols = "database_version"
            else:
                db_cols = "'' as database_version"

            q = (
                f"""
                SELECT DISTINCT
                    CAST({anv} AS analyses) as analysis,
                    software_version,
                    {db_cols}
                FROM results_{ans}
                """
            )
            query.append(q)

        result = self.cur.execute("\nUNION\n".join(query))
        for r in result:
            yield TargetRow.from_rowfactory(r)
        return
