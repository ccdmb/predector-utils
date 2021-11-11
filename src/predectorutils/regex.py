#!/usr/bin/env python3

import re
from typing import Pattern, Iterator
from .analyses import RegexAnalysis

from Bio.SeqRecord import SeqRecord

FLAGS = re.ASCII | re.IGNORECASE
KEX2_OUTRAM_2021 = r"[LIJVAP][A-Z][KRTPEI]R"
KEX2_OUTRAM_2020 = r"L[A-Z][A-Z]R"
KEX2_CLASSIC = r"[KR]R"

RXLRLIKE = r"[RKH][A-Z][LMIFYW][A-Z]"


class RegexMatcher(object):

    def __init__(self, regex: Pattern, kind: str):
        self.regex = regex
        self.kind = kind
        return

    def find(self, seq: SeqRecord) -> Iterator[RegexAnalysis]:
        matches = self.regex.finditer(str(seq.seq))

        for match in matches:
            string = match.group()
            start, stop = match.span()
            yield RegexAnalysis(
                seq.id,
                self.kind,
                self.regex.pattern,
                string,
                start,
                stop
            )

        return

    def findmany(self, seqs) -> Iterator[RegexAnalysis]:
        for seq in seqs:
            for match in self.find(seq):
                yield match
        return
