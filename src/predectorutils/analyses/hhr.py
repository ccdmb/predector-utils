#!/usr/bin/env python3

import re

from typing import TextIO
from typing import Iterator
from typing import Optional
from typing import Sequence, List
from typing import Mapping
from typing import Tuple
from typing import TypeVar, Callable

from predectorutils.higher import fmap
from predectorutils.analyses.base import Analysis
from predectorutils.analyses.base import float_or_none
from predectorutils.analyses.parsers import (
    ParseError,
    LineParseError,
    BlockParseError
)

from predectorutils.analyses.parsers import (
    parse_int,
    parse_float,
    split_at_eq,
    split_at_multispace,
    parse_string_not_empty,
    get_from_dict_or_err,
    is_one_of,
    MULTISPACE_REGEX,
)

T = TypeVar("T")
ALI_REGEX = re.compile(r"[QT]\s+[^\s]+\s+\d+\s+")


def get_and_parse(
    key: str,
    field_name: str,
    d: Mapping[str, str],
    func: Callable[[str, str], T]
) -> T:
    return func(
        get_from_dict_or_err(key, d, field_name),
        field_name
    )


class HHRAlignment(Analysis):

    """ """

    columns = [
        'query_id',
        'query_length',
        'query_neff',
        'template_id',
        'template_length',
        'template_info',
        'template_neff',
        'query_ali',
        'template_ali',
        'confidence',
        'query_start',
        'template_start',
        'query_end',
        'template_end',
        'probability',
        'evalue',
        'score',
        'aligned_cols',
        'identity',
        'similarity',
        'sum_probs'
    ]

    types = [
        str,
        int,
        float,
        str,
        int,
        str,
        float_or_none,
        str,
        str,
        str,
        int,
        int,
        int,
        int,
        float,
        float,
        float,
        int,
        float,
        float,
        float,
    ]
    analysis = "hhblits"
    name_column = "query_id"
    software = "HHSuite"

    def __init__(
        self,
        query_id: str,
        query_length: int,
        query_neff: float,
        template_id: str,
        template_length: int,
        template_info: str,
        template_neff: Optional[float],
        query_ali: str,
        template_ali: str,
        confidence: str,
        query_start: int,
        template_start: int,
        query_end: int,
        template_end: int,
        probability: float,
        evalue: float,
        score: float,
        aligned_cols: int,
        identity: float,
        similarity: float,
        sum_probs: float,
    ):
        self.query_id = query_id
        self.query_length = query_length
        self.query_neff = query_neff
        self.template_id = template_id
        self.template_length = template_length
        self.template_info = template_info
        self.template_neff = template_neff
        self.query_ali = query_ali
        self.template_ali = template_ali
        self.confidence = confidence
        self.query_start = query_start
        self.template_start = template_start
        self.query_end = query_end
        self.template_end = template_end
        self.probability = probability
        self.evalue = evalue
        self.score = score
        self.aligned_cols = aligned_cols
        self.identity = identity
        self.similarity = similarity
        self.sum_probs = sum_probs
        return

    @classmethod
    def from_block(cls, lines: Sequence[str]) -> Iterator["HHRAlignment"]:
        if len(lines) == 0:
            raise BlockParseError(0, "The block was empty.")

        query: Optional[str] = None
        query_length: Optional[int] = None
        query_neff: Optional[float] = None

        is_alignment = False
        alignment_block: List[str] = []

        for i, line in enumerate(lines):
            sline = line.rstrip('\n')

            if sline.startswith(">") and is_alignment:
                try:
                    yield cls._parse_alignment(
                        alignment_block,
                        query,
                        query_length,
                        query_neff
                    )
                except BlockParseError as e:
                    raise BlockParseError.from_block_error(e, i)

                alignment_block = []

            elif sline.startswith(">"):
                is_alignment = True

            if is_alignment:
                alignment_block.append(sline)
                continue

            try:
                if sline.startswith("Query"):
                    query = cls._parse_query_line(sline)
                elif sline.startswith("Match_columns"):
                    query_length = cls._parse_query_length_line(sline)
                elif sline.startswith("Neff"):
                    query_neff = cls._parse_query_neff_line(sline)
            except LineParseError as e:
                raise BlockParseError.from_line_error(e, i)

        if len(alignment_block) > 0:
            try:
                yield cls._parse_alignment(
                    alignment_block,
                    query,
                    query_length,
                    query_neff
                )
            except BlockParseError as e:
                raise BlockParseError.from_block_error(e, i)
        return

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["HHRAlignment"]:
        block: List[str] = []

        for i, line in enumerate(handle):
            sline = line.rstrip('\n')

            if sline == "":
                continue

            try:
                if sline.startswith("Query") and len(block) > 0:
                    for b in cls.from_block(block):
                        yield b
                    block = []

                block.append(sline)

            except BlockParseError as e:
                raise ParseError.from_block_error(e, i, handle)

        try:
            if len(block) > 0:
                for b in cls.from_block(block):
                    yield b

        except BlockParseError as e:
            raise ParseError.from_block_error(e, i, handle)
        return

    @staticmethod
    def _is_not_none(val: Optional[T], field_name: str) -> T:
        if val is None:
            raise LineParseError(
                f"Did not encounter {field_name} in alignment."
            )
        return val

    @staticmethod
    def _is_not_empty(val: List[T], field_name: str) -> List[T]:
        if len(val) == 0:
            raise LineParseError(
                f"Did not encounter {field_name} in alignment."
            )
        return val

    @classmethod
    def _parse_alignment(
        cls,
        block: Sequence[str],
        query_id: Optional[str],
        query_length: Optional[int],
        query_neff: Optional[float],
    ) -> "HHRAlignment":
        if ((query_id is None) or
                (query_length is None) or
                (query_neff is None)):
            raise BlockParseError(
                0,
                "Reached an alignment block before the query information."
            )

        skip_ali_tags = ("ss_dssp", "ss_pred", "Consensus")

        template_id: Optional[str] = None
        template_info: Optional[str] = None
        query_starts: List[int] = []
        query_ends: List[int] = []
        query_sequence: List[str] = []
        template_starts: List[int] = []
        template_ends: List[int] = []
        template_sequence: List[str] = []
        template_length: Optional[int] = None
        confidence_sequence: List[str] = []
        seq_begin_col: Optional[int] = None

        probability: Optional[float] = None
        evalue: Optional[float] = None
        score: Optional[float] = None
        identity: Optional[float] = None
        similarity: Optional[float] = None
        template_neff: Optional[float] = None
        sum_probs: Optional[float] = None
        aligned_cols: Optional[int] = None

        for i, line in enumerate(block):
            if line.startswith(">"):
                template_id = cls._parse_sequence_name(line)
                template_info = line

            elif line.startswith("Q"):
                # type id ali_start sequence ali_end length
                try:
                    query_line = cls._parse_alignment_line(line)
                except LineParseError as e:
                    raise BlockParseError.from_line_error(e, i)

                if query_line[1] in skip_ali_tags:
                    continue

                query_starts.append(query_line[2])
                query_ends.append(query_line[4])
                query_sequence.append(query_line[3])

                if query_line[6] is not None:
                    seq_begin_col = query_line[6]

            elif line.startswith("T"):
                try:
                    template_line = cls._parse_alignment_line(line)
                except LineParseError as e:
                    raise BlockParseError.from_line_error(e, i)

                if template_line[1] in skip_ali_tags:
                    continue

                template_starts.append(template_line[2])
                template_ends.append(template_line[4])
                template_sequence.append(template_line[3])
                template_length = template_line[5]

                if template_line[6] is not None:
                    seq_begin_col = template_line[6]

            elif line.startswith("Probab"):
                try:
                    (probability,
                     evalue,
                     score,
                     aligned_cols,
                     identity,
                     similarity,
                     sum_probs,
                     template_neff) = cls._parse_probab_line(line)
                except LineParseError as e:
                    raise BlockParseError.from_line_error(e, i)

            elif line.startswith("Confidence"):
                if seq_begin_col is None:
                    raise BlockParseError(i, (
                        "Encountered a 'Confidence' line "
                        "before the alignments themselves."
                    ))

                actual_lhs = line[:seq_begin_col].rstrip()
                rhs = line[seq_begin_col:]
                if not "Confidence" == actual_lhs:
                    print(actual_lhs)
                    raise BlockParseError(
                        i,
                        f"Could not parse the confidence line at "
                    )

                confidence_sequence.append(rhs)

        query_ali = "".join(query_sequence)
        template_ali = "".join(template_sequence)
        confidence = "".join(confidence_sequence).replace(" ", "-")
        assert len(query_ali) == len(template_ali) == len(confidence), \
            (template_id, repr(confidence))

        query_start = min(cls._is_not_empty(query_starts, "query_start"))
        query_end = max(cls._is_not_empty(query_ends, "query_ends"))

        template_start = min(cls._is_not_empty(
            template_starts,
            "template_start"
        ))
        template_end = max(cls._is_not_empty(template_ends, "template_ends"))

        return cls(
            query_id,
            query_length,
            query_neff,
            cls._is_not_none(template_id, "template_id"),
            cls._is_not_none(template_length, "template_length"),
            cls._is_not_none(template_info, "template_info"),
            template_neff,
            query_ali,
            template_ali,
            confidence,
            query_start,
            template_start,
            query_end,
            template_end,
            cls._is_not_none(probability, "probability"),
            cls._is_not_none(evalue, "evalue"),
            cls._is_not_none(score, "score"),
            cls._is_not_none(aligned_cols, "aligned_cols"),
            cls._is_not_none(identity, "identity"),
            cls._is_not_none(similarity, "similarity"),
            cls._is_not_none(sum_probs, "sum_probs")
        )

    @staticmethod
    def _parse_query_line(field: str) -> str:
        return split_at_multispace(field, "query", "Query")

    @staticmethod
    def _parse_query_length_line(field: str) -> int:
        return parse_int(
            split_at_multispace(field, "query_length", "Match_columns"),
            "query_length",
        )

    @staticmethod
    def _parse_query_neff_line(field: str) -> float:
        return parse_float(
            split_at_multispace(field, "query_neff", "Neff"),
            "query_neff",
        )

    @staticmethod
    def _parse_probab_line(
        field: str
    ) -> Tuple[float, float, float, int, float, float, float, Optional[float]]:
        sline = (s for s in MULTISPACE_REGEX.split(field.strip()))
        columns = [
            "Probab",
            "E-value",
            "Score",
            "Aligned_cols",
            "Identities",
            "Similarity",
            "Sum_probs",
            "Template_Neff",
        ]

        dline = {
            col: split_at_eq(f, col, col)
            for f, col
            in zip(sline, columns)
        }

        if "Template_Neff" in dline:
            template_neff: Optional[float] = parse_float(
                dline["Template_Neff"],
                "template_neff"
               )
        else:
            template_neff = None

        return (
            get_and_parse("Probab", "probability", dline, parse_float),
            get_and_parse("E-value", "evalue", dline, parse_float),
            get_and_parse("Score", "score", dline, parse_float),
            get_and_parse("Aligned_cols", "aligned_cols", dline, parse_int),
            get_and_parse(
                "Identities",
                "identity",
                dline,
                lambda x, y: parse_float(x.rstrip("%"), y) / 100.0
            ),
            get_and_parse("Similarity", "similarity", dline, parse_float),
            get_and_parse("Sum_probs", "sum_probs", dline, parse_float),
            template_neff,
        )

    @staticmethod
    def _parse_sequence_name(header: str) -> str:
        name = header.replace(">", "").split()[0]
        return name

    @staticmethod
    def _parse_alignment_line(
        line: str
    ) -> Tuple[str, str, int, str, int, int, Optional[int]]:
        sline = MULTISPACE_REGEX.split(line.strip(), maxsplit=5)

        columns = ["type", "id", "ali_start", "sequence", "ali_end", "length"]
        dline = dict(zip(columns, sline))

        length = fmap(
            lambda x: x.lstrip("(").rstrip(")"),
            dline.get("length", None)
        )

        if length is None:
            raise LineParseError(
                f"Missing 'length' from alignment line: '{line}'.")

        seq_begin_match = ALI_REGEX.match(line)
        if seq_begin_match is None:
            seq_begin: Optional[int] = None
        else:
            seq_begin = seq_begin_match.end()

        return (
            get_and_parse(
                "type",
                "type",
                dline,
                lambda f, fn: is_one_of(f, ["T", "Q"], fn)
            ),
            get_and_parse("id", "id", dline, parse_string_not_empty),
            get_and_parse("ali_start", "ali_start", dline, parse_int),
            get_and_parse("sequence", "sequence",
                          dline, parse_string_not_empty),
            get_and_parse("ali_end", "ali_end", dline, parse_int),
            parse_int(length, "length"),
            seq_begin
        )
