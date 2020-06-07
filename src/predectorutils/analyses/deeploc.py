#!/usr/bin/env python3

from typing import TextIO
from typing import Iterator

from predectorutils.analyses.base import Analysis
from predectorutils.parsers import (
    FieldParseError,
    LineParseError,
    parse_field,
    raise_it,
    parse_str,
    parse_float,
    is_one_of
)

dl_name = raise_it(parse_field(parse_str, "name"))
dl_prediction = raise_it(parse_field(
    is_one_of([
        "Membrane", "Nucleus", "Cytoplasm", "Extracellular",
        "Mitochondrion", "Cell_membrane", "Endoplasmic_reticulum",
        "Plastid", "Golgi_apparatus", "Lysosome/Vacuole",
        "Peroxisome"
    ]),
    "prediction"
))
dl_membrane = raise_it(parse_field(parse_float, "membrane"))
dl_nucleus = raise_it(parse_field(parse_float, "nucleus"))
dl_cytoplasm = raise_it(parse_field(parse_float, "cytoplasm"))
dl_extracellular = raise_it(parse_field(parse_float, "extracellular"))
dl_mitochondrion = raise_it(parse_field(parse_float, "mitochondrion"))
dl_cell_membrane = raise_it(parse_field(parse_float, "cell_membrane"))
dl_endoplasmic_reticulum = raise_it(parse_field(
    parse_float,
    "endoplasmic_reticulum"
))
dl_plastid = raise_it(parse_field(parse_float, "plastid"))
dl_golgi_apparatus = raise_it(parse_field(parse_float, "golgi_apparatus"))
dl_lysosome = raise_it(parse_field(parse_float, "lysosome_vacuole"))
dl_peroxisome = raise_it(parse_field(parse_float, "peroxisome"))


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

        return cls(
            dl_name(sline[0]),
            dl_prediction(sline[1]),
            dl_membrane(sline[2]),
            dl_nucleus(sline[3]),
            dl_cytoplasm(sline[4]),
            dl_extracellular(sline[5]),
            dl_mitochondrion(sline[6]),
            dl_cell_membrane(sline[7]),
            dl_endoplasmic_reticulum(sline[8]),
            dl_plastid(sline[9]),
            dl_golgi_apparatus(sline[10]),
            dl_lysosome(sline[11]),
            dl_peroxisome(sline[12]),
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

            except (LineParseError, FieldParseError) as e:
                raise e.as_parse_error(line=i).add_filename_from_handle(handle)
        return
