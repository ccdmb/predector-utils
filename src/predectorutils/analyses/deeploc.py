#!/usr/bin/env python3

from typing import Optional
from typing import TextIO
from typing import Iterator

from predectorutils.analyses.base import Analysis
from predectorutils.analyses.parsers import ParseError, LineParseError
from predectorutils.analyses.parsers import (
    parse_string_not_empty,
    parse_float,
    is_one_of
)


class DeepLoc(Analysis):

    """ Doesn't have output format documentation yet
    """

    columns = ["name", "prediction", "membrane", "nucleus", "cytoplasm",
               "extracellular", "mitochondrion", "cell_membrane",
               "endoplasmic_reticulum", "plastid", "golgi_apparatus",
               "lysosome_vacuole", "peroxisome"]
    types = [str, str, float, float, float, float,
             float, float, float, float, float, float, float]
    analysis = "deeploc"
    software = "DeepLoc"

    def __init__(
        self,
        name: str,
        prediction: str,
        membrane: float,
        nucleus: float,
        cytoplasm: float,
        extracellular: float,
        mitochondrion: float,
        cell_membrane: float,
        endoplasmic_reticulum: float,
        plastid: float,
        golgi_apparatus: float,
        lysosome_vacuole: float,
        peroxisome: float,
    ) -> None:
        self.name = name
        self.prediction = prediction
        self.membrane = membrane
        self.nucleus = nucleus
        self.cytoplasm = cytoplasm
        self.extracellular = extracellular
        self.mitochondrion = mitochondrion
        self.cell_membrane = cell_membrane
        self.endoplasmic_reticulum = endoplasmic_reticulum
        self.plastid = plastid
        self.golgi_apparatus = golgi_apparatus
        self.lysosome_vacuole = lysosome_vacuole
        self.peroxisome = peroxisome
        return

    @classmethod
    def from_line(cls, line: str) -> "DeepLoc":
        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) != 13:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 13 but got {len(sline)}"
            )

        prediction = is_one_of(
            sline[1],
            [
                "Membrane", "Nucleus", "Cytoplasm", "Extracellular",
                "Mitochondrion", "Cell_membrane", "Endoplasmic_reticulum",
                "Plastid", "Golgi_apparatus", "Lysosome/Vacuole",
                "Peroxisome"
            ],
            "prediction"
        )

        if prediction == "noTP":
            prediction = "OTHER"

        return cls(
            parse_string_not_empty(sline[0], "name"),
            prediction,
            parse_float(sline[2], "membrane"),
            parse_float(sline[3], "nucleus"),
            parse_float(sline[4], "cytoplasm"),
            parse_float(sline[5], "extracellular"),
            parse_float(sline[6], "mitochondrion"),
            parse_float(sline[7], "cell_membrane"),
            parse_float(sline[8], "endoplasmic_reticulum"),
            parse_float(sline[9], "plastid"),
            parse_float(sline[10], "golgi_apparatus"),
            parse_float(sline[11], "lysosome_vacuole"),
            parse_float(sline[12], "peroxisome"),
        )

    @classmethod
    def from_file(
        cls,
        handle: TextIO,
    ) -> Iterator["DeepLoc"]:
        for i, line in enumerate(handle):
            sline = line.strip()
            if sline.startswith("#"):
                continue
            elif sline == "":
                continue
            elif sline.startswith("ID	Location	Membrane"):
                continue

            try:
                yield cls.from_line(sline)

            except LineParseError as e:
                if hasattr(handle, "name"):
                    filename: Optional[str] = handle.name
                else:
                    filename = None

                raise ParseError(
                    filename,
                    i,
                    e.message
                )
        return
