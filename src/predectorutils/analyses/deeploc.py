#!/usr/bin/env python3

from typing import TextIO, Callable, Union
from collections.abc import Iterator

from .base import Analysis, list_of_str
from ..parsers import (
    FieldParseError,
    LineParseError,
    parse_field,
    raise_it,
    parse_str,
    parse_float,
    is_one_of,
    parse_delim
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



dl2_name = raise_it(parse_field(parse_str, "Protein_ID"))
dl2_localisations = raise_it(
    parse_field(
        parse_delim(r"\|", is_one_of([
            "Cytoplasm", "Nucleus", "Extracellular",
            "Cell membrane", "Mitochondrion", "Plastid",
            "Endoplasmic reticulum", "Lysosome/Vacuole", "Golgi apparatus",
            "Peroxisome"
        ])),
        "Localizations"
    )
)

dl2_signals = raise_it(
    parse_field(
        parse_delim(r"\|", is_one_of([
            "Signal peptide", "Transmembrane domain", "Mitochondrial transit peptide",
            "Chloroplast transit peptide", "Thylakoid luminal transit peptide",
            "Nuclear localization signal", "Nuclear export signal", "Peroxisomal targeting signal"
        ])),
        "Signals"
    )
)

dl2_cytoplasm = raise_it(parse_field(parse_float, "Cytoplasm"))
dl2_nucleus = raise_it(parse_field(parse_float, "Nucleus"))
dl2_extracellular = raise_it(parse_field(parse_float, "Extracellular"))
dl2_cell_membrane = raise_it(parse_field(parse_float, "Cell membrane"))
dl2_mitochondrion = raise_it(parse_field(parse_float, "Mitochondrion"))
dl2_plastid = raise_it(parse_field(parse_float, "Plastid"))
dl2_endoplasmic_reticulum = raise_it(parse_field(parse_float, "Endoplasmic reticulum"))
dl2_lysosome = raise_it(parse_field(parse_float, "Lysosome/Vacuole"))
dl2_golgi_apparatus = raise_it(parse_field(parse_float, "Golgi apparatus"))
dl2_peroxisome = raise_it(parse_field(parse_float, "Peroxisome"))


class DeepLoc2(Analysis):

    """ Doesn't have output format documentation yet
    """

    columns = ["name", "localisations", "signals", "cytoplasm", "nucleus",
               "extracellular", "cell_membrane", "mitochondrion",
               "plastid", "endoplasmic_reticulum",
               "lysosome_vacuole", "golgi_apparatus", "peroxisome"]
    types = [str, list_of_str, list_of_str, float, float, float,
             float, float, float, float, float, float, float]
    analysis = "deeploc2"
    software = "DeepLoc2"

    def __init__(
        self,
        name: str,
        localisations: list[str],
        signals: list[str],
        cytoplasm: float,
        nucleus: float,
        extracellular: float,
        cell_membrane: float,
        mitochondrion: float,
        plastid: float,
        endoplasmic_reticulum: float,
        lysosome_vacuole: float,
        golgi_apparatus: float,
        peroxisome: float,
    ) -> None:
        self.name = name
        self.localisations = localisations
        self.signals = signals
        self.cytoplasm = cytoplasm
        self.nucleus = nucleus
        self.extracellular = extracellular
        self.cell_membrane = cell_membrane
        self.mitochondrion = mitochondrion
        self.plastid = plastid
        self.endoplasmic_reticulum = endoplasmic_reticulum
        self.lysosome_vacuole = lysosome_vacuole
        self.golgi_apparatus = golgi_apparatus
        self.peroxisome = peroxisome
        return

    @classmethod
    def from_line(cls, line: str) -> "DeepLoc2":
        if line == "":
            raise LineParseError("The line was empty.")

        sline = list(map(str.strip, line.strip().split(",")))

        if len(sline) != 13:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 13 but got {len(sline)}"
            )

        return cls(
            dl2_name(sline[0]),
            dl2_localisations(sline[1]),
            dl2_signals(sline[2]),
            dl2_cytoplasm(sline[3]),
            dl2_nucleus(sline[4]),
            dl2_extracellular(sline[5]),
            dl2_cell_membrane(sline[6]),
            dl2_mitochondrion(sline[7]),
            dl2_plastid(sline[8]),
            dl2_endoplasmic_reticulum(sline[9]),
            dl2_lysosome(sline[10]),
            dl2_golgi_apparatus(sline[11]),
            dl2_peroxisome(sline[12]),
        )

    @classmethod
    def from_file(
        cls,
        handle: TextIO,
    ) -> Iterator["DeepLoc2"]:
        for i, line in enumerate(handle):
            sline = line.strip()
            if sline.startswith("#"):
                continue
            elif sline == "":
                continue
            elif sline.startswith("Protein_ID,Localizations,Signals"):
                continue

            try:
                yield cls.from_line(sline)

            except (LineParseError, FieldParseError) as e:
                raise e.as_parse_error(line=i).add_filename_from_handle(handle)
        return
