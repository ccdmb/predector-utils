#!/usr/bin/env python3

import re
from typing import TextIO
from typing import Iterator
from typing import Dict
from typing import TypeVar
from typing import Callable
from typing import Tuple
from typing import List, Sequence
from typing import Pattern

from typing import Iterable

from predectorutils.analyses.base import Analysis
from predectorutils.analyses.parsers import (
    ParseError,
    BlockParseError,
    LineParseError
)

from predectorutils.analyses.parsers import (
    parse_regex,
)


T = TypeVar("T")


def convert_line_err(
    lineno: int,
    field: str,
    parser: Callable[[str], T]
) -> T:
    try:
        return parser(field)
    except LineParseError as e:
        raise BlockParseError.from_line_error(e, lineno)


def get_line(lines: Iterator[Tuple[int, str]]) -> Tuple[int, str]:
    i, line = next(lines)

    while line.strip() == "":
        i, line = next(lines)

    return i, line.strip()


def parse_regex_line(
    lines: Iterable[Tuple[int, str]],
    regex: Pattern,
    record: Dict[str, str]
) -> None:
    i, line = get_line(iter(lines))
    record.update(
        convert_line_err(
            i,
            line,
            lambda l: parse_regex(l, regex)
        )
    )
    return


def parse_residue_table(lines: Iterable[Tuple[int, str]]) -> Dict[str, str]:
    table: Dict[str, str] = dict()

    for i, line in lines:
        if line == "":
            break

        dline = convert_line_err(
            i,
            line,
            lambda l: parse_regex(l, RESIDUE_REGEX)
        )

        aa = dline['aa'].lower()
        table[f"residue_{aa}_number"] = dline["number"]
        table[f"residue_{aa}_mole"] = dline["mole"]
        table[f"residue_{aa}_dayhoff"] = dline["dayhoff"]

    return table


def parse_property_table(lines: Iterable[Tuple[int, str]]) -> Dict[str, str]:
    table: Dict[str, str] = dict()

    for i, line in lines:
        if line == "":
            break

        dline = convert_line_err(
            i,
            line,
            lambda l: parse_regex(l, PROP_REGEX)
        )

        kind = dline['kind'].lower().replace("-", "")
        table[f"property_{kind}_number"] = dline["number"]
        table[f"property_{kind}_mole"] = dline["mole"]

    return table


NAME_REGEX = re.compile(r"^PEPSTATS of\s+(?P<name>[^\s]+)")
MOLWGT_REGEX = re.compile(
    r"^Molecular weight\s+=\s+(?P<molecular_weight>[-+]?\d*\.?\d+)\s+"
    r"Residues\s+=\s+(?P<residues>\d+)"
)
AVWGT_REGEX = re.compile(
    r"Average\s+Residue\s+Weight\s+=\s+"
    r"(?P<average_residue_weight>[-+]?\d*\.?\d+)\s+"
    r"Charge\s+=\s+(?P<charge>[-+]?\d*\.?\d+)"
)
ISOELEC_REGEX = re.compile(
    r"Isoelectric\s+Point\s+=\s+(?P<isoelectric_point>[-+]?\d*\.?\d+)"
)
A280_MOLAR_REGEX = re.compile(
    r"A280\s+Molar\s+Extinction\s+Coefficients\s+=\s+"
    r"(?P<a280_molar_extinction_coefficients_reduced>\d+)\s+\(reduced\)\s+"
    r"(?P<a280_molar_extinction_coefficients_cysteine_bridges>\d+)\s+"
    r"\(cystine bridges\)"
)
A280_1MGML_REGEX = re.compile(
    r"A280\s+Extinction\s+Coefficients\s+1mg/ml\s+=\s+"
    r"(?P<a280_1mgml_extinction_coefficients_reduced>[-+]?\d*\.?\d+)\s+"
    r"\(reduced\)\s+"
    r"(?P<a280_1mgml_extinction_coefficients_cysteine_bridges>"
    r"[-+]?\d*\.?\d+)\s+"
    r"\(cystine bridges\)"
)
IMPROB_REGEX = re.compile(
    r"(Imp)?P?robability\s+of\s+expression\s+in\s+inclusion\s+bodies\s+=\s+"
    r"(?P<improbability_expression_inclusion_bodies>[-+]?\d*\.?\d+)"
)
RESIDUE_REGEX = re.compile(
    r"(?P<aa>[a-zA-Z])\s+=\s+.+\s+"
    r"(?P<number>\d+)\s+"
    r"(?P<mole>[-+]?\d*\.?\d+)\s+"
    r"(?P<dayhoff>[-+]?\d*\.?\d+)"
)
PROP_REGEX = re.compile(
    r"(?P<kind>[^\s]+)\s+.+\s+"
    r"(?P<number>\d+)\s+"
    r"(?P<mole>[-+]?\d*\.?\d+)"
)


class PepStats(Analysis):

    """ """

    columns = [
        "name",
        "molecular_weight",
        "residues",
        "average_residue_weight",
        "charge",
        "isoelectric_point",
        "a280_molar_extinction_coefficients_reduced",
        "a280_molar_extinction_coefficients_cysteine_bridges",
        "a280_1mgml_extinction_coefficients_reduced",
        "a280_1mgml_extinction_coefficients_cysteine_bridges",
        "improbability_expression_inclusion_bodies",
        "residue_a_number",
        "residue_b_number",
        "residue_c_number",
        "residue_d_number",
        "residue_e_number",
        "residue_f_number",
        "residue_g_number",
        "residue_h_number",
        "residue_i_number",
        "residue_j_number",
        "residue_k_number",
        "residue_l_number",
        "residue_m_number",
        "residue_n_number",
        "residue_o_number",
        "residue_p_number",
        "residue_q_number",
        "residue_r_number",
        "residue_s_number",
        "residue_t_number",
        "residue_u_number",
        "residue_v_number",
        "residue_w_number",
        "residue_x_number",
        "residue_y_number",
        "residue_z_number",
        "residue_a_mole",
        "residue_b_mole",
        "residue_c_mole",
        "residue_d_mole",
        "residue_e_mole",
        "residue_f_mole",
        "residue_g_mole",
        "residue_h_mole",
        "residue_i_mole",
        "residue_j_mole",
        "residue_k_mole",
        "residue_l_mole",
        "residue_m_mole",
        "residue_n_mole",
        "residue_o_mole",
        "residue_p_mole",
        "residue_q_mole",
        "residue_r_mole",
        "residue_s_mole",
        "residue_t_mole",
        "residue_u_mole",
        "residue_v_mole",
        "residue_w_mole",
        "residue_x_mole",
        "residue_y_mole",
        "residue_z_mole",
        "residue_a_dayhoff",
        "residue_b_dayhoff",
        "residue_c_dayhoff",
        "residue_d_dayhoff",
        "residue_e_dayhoff",
        "residue_f_dayhoff",
        "residue_g_dayhoff",
        "residue_h_dayhoff",
        "residue_i_dayhoff",
        "residue_j_dayhoff",
        "residue_k_dayhoff",
        "residue_l_dayhoff",
        "residue_m_dayhoff",
        "residue_n_dayhoff",
        "residue_o_dayhoff",
        "residue_p_dayhoff",
        "residue_q_dayhoff",
        "residue_r_dayhoff",
        "residue_s_dayhoff",
        "residue_t_dayhoff",
        "residue_u_dayhoff",
        "residue_v_dayhoff",
        "residue_w_dayhoff",
        "residue_x_dayhoff",
        "residue_y_dayhoff",
        "residue_z_dayhoff",
        "property_tiny_number",
        "property_small_number",
        "property_aliphatic_number",
        "property_aromatic_number",
        "property_nonpolar_number",
        "property_polar_number",
        "property_charged_number",
        "property_basic_number",
        "property_acidic_number",
        "property_tiny_mole",
        "property_small_mole",
        "property_aliphatic_mole",
        "property_aromatic_mole",
        "property_nonpolar_mole",
        "property_polar_mole",
        "property_charged_mole",
        "property_basic_mole",
        "property_acidic_mole",
    ]
    types = [
        str, float, int, float, float, float, int, float, float, float, float,
        int, int, int, int, int, int, int, int, int, int, int, int, int, int,
        int, int, int, int, int, int, int, int, int, int, int, int,
        float, float, float, float, float, float, float, float, float, float,
        float, float, float, float, float, float, float, float, float, float,
        float, float, float, float, float, float,
        float, float, float, float, float, float, float, float, float, float,
        float, float, float, float, float, float, float, float, float, float,
        float, float, float, float, float, float,
        int, int, int, int, int, int, int, int, int,
        float, float, float, float, float, float, float, float, float,
    ]
    analysis = "pepstats"
    software = "EMBOSS"

    def __init__(
        self,
        name: str,
        molecular_weight: float,
        residues: int,
        average_residue_weight: float,
        charge: float,
        isoelectric_point: float,
        a280_molar_extinction_coefficients_reduced: int,
        a280_molar_extinction_coefficients_cysteine_bridges: float,
        a280_1mgml_extinction_coefficients_reduced: float,
        a280_1mgml_extinction_coefficients_cysteine_bridges: float,
        improbability_expression_inclusion_bodies: float,
        residue_a_number: int,
        residue_b_number: int,
        residue_c_number: int,
        residue_d_number: int,
        residue_e_number: int,
        residue_f_number: int,
        residue_g_number: int,
        residue_h_number: int,
        residue_i_number: int,
        residue_j_number: int,
        residue_k_number: int,
        residue_l_number: int,
        residue_m_number: int,
        residue_n_number: int,
        residue_o_number: int,
        residue_p_number: int,
        residue_q_number: int,
        residue_r_number: int,
        residue_s_number: int,
        residue_t_number: int,
        residue_u_number: int,
        residue_v_number: int,
        residue_w_number: int,
        residue_x_number: int,
        residue_y_number: int,
        residue_z_number: int,
        residue_a_mole: float,
        residue_b_mole: float,
        residue_c_mole: float,
        residue_d_mole: float,
        residue_e_mole: float,
        residue_f_mole: float,
        residue_g_mole: float,
        residue_h_mole: float,
        residue_i_mole: float,
        residue_j_mole: float,
        residue_k_mole: float,
        residue_l_mole: float,
        residue_m_mole: float,
        residue_n_mole: float,
        residue_o_mole: float,
        residue_p_mole: float,
        residue_q_mole: float,
        residue_r_mole: float,
        residue_s_mole: float,
        residue_t_mole: float,
        residue_u_mole: float,
        residue_v_mole: float,
        residue_w_mole: float,
        residue_x_mole: float,
        residue_y_mole: float,
        residue_z_mole: float,
        residue_a_dayhoff: float,
        residue_b_dayhoff: float,
        residue_c_dayhoff: float,
        residue_d_dayhoff: float,
        residue_e_dayhoff: float,
        residue_f_dayhoff: float,
        residue_g_dayhoff: float,
        residue_h_dayhoff: float,
        residue_i_dayhoff: float,
        residue_j_dayhoff: float,
        residue_k_dayhoff: float,
        residue_l_dayhoff: float,
        residue_m_dayhoff: float,
        residue_n_dayhoff: float,
        residue_o_dayhoff: float,
        residue_p_dayhoff: float,
        residue_q_dayhoff: float,
        residue_r_dayhoff: float,
        residue_s_dayhoff: float,
        residue_t_dayhoff: float,
        residue_u_dayhoff: float,
        residue_v_dayhoff: float,
        residue_w_dayhoff: float,
        residue_x_dayhoff: float,
        residue_y_dayhoff: float,
        residue_z_dayhoff: float,
        property_tiny_number: int,
        property_small_number: int,
        property_aliphatic_number: int,
        property_aromatic_number: int,
        property_nonpolar_number: int,
        property_polar_number: int,
        property_charged_number: int,
        property_basic_number: int,
        property_acidic_number: int,
        property_tiny_mole: float,
        property_small_mole: float,
        property_aliphatic_mole: float,
        property_aromatic_mole: float,
        property_nonpolar_mole: float,
        property_polar_mole: float,
        property_charged_mole: float,
        property_basic_mole: float,
        property_acidic_mole: float,
    ) -> None:
        self.name = name
        self.molecular_weight = molecular_weight
        self.residues = residues
        self.average_residue_weight = average_residue_weight
        self.charge = charge
        self.isoelectric_point = isoelectric_point
        self.a280_molar_extinction_coefficients_reduced = a280_molar_extinction_coefficients_reduced  # noqa
        self.a280_molar_extinction_coefficients_cysteine_bridges = a280_1mgml_extinction_coefficients_cysteine_bridges  # noqa
        self.a280_1mgml_extinction_coefficients_reduced = a280_1mgml_extinction_coefficients_reduced  # noqa
        self.a280_1mgml_extinction_coefficients_cysteine_bridges = a280_1mgml_extinction_coefficients_cysteine_bridges  # noqa
        self.improbability_expression_inclusion_bodies = improbability_expression_inclusion_bodies  # noqa
        self.residue_a_number = residue_a_number
        self.residue_b_number = residue_b_number
        self.residue_c_number = residue_c_number
        self.residue_d_number = residue_d_number
        self.residue_e_number = residue_e_number
        self.residue_f_number = residue_f_number
        self.residue_g_number = residue_g_number
        self.residue_h_number = residue_h_number
        self.residue_i_number = residue_i_number
        self.residue_j_number = residue_j_number
        self.residue_k_number = residue_k_number
        self.residue_l_number = residue_l_number
        self.residue_m_number = residue_m_number
        self.residue_n_number = residue_n_number
        self.residue_o_number = residue_o_number
        self.residue_p_number = residue_p_number
        self.residue_q_number = residue_q_number
        self.residue_r_number = residue_r_number
        self.residue_s_number = residue_s_number
        self.residue_t_number = residue_t_number
        self.residue_u_number = residue_u_number
        self.residue_v_number = residue_v_number
        self.residue_w_number = residue_w_number
        self.residue_x_number = residue_x_number
        self.residue_y_number = residue_y_number
        self.residue_z_number = residue_z_number
        self.residue_a_mole = residue_a_mole
        self.residue_b_mole = residue_b_mole
        self.residue_c_mole = residue_c_mole
        self.residue_d_mole = residue_d_mole
        self.residue_e_mole = residue_e_mole
        self.residue_f_mole = residue_f_mole
        self.residue_g_mole = residue_g_mole
        self.residue_h_mole = residue_h_mole
        self.residue_i_mole = residue_i_mole
        self.residue_j_mole = residue_j_mole
        self.residue_k_mole = residue_k_mole
        self.residue_l_mole = residue_l_mole
        self.residue_m_mole = residue_m_mole
        self.residue_n_mole = residue_n_mole
        self.residue_o_mole = residue_o_mole
        self.residue_p_mole = residue_p_mole
        self.residue_q_mole = residue_q_mole
        self.residue_r_mole = residue_r_mole
        self.residue_s_mole = residue_s_mole
        self.residue_t_mole = residue_t_mole
        self.residue_u_mole = residue_u_mole
        self.residue_v_mole = residue_v_mole
        self.residue_w_mole = residue_w_mole
        self.residue_x_mole = residue_x_mole
        self.residue_y_mole = residue_y_mole
        self.residue_z_mole = residue_z_mole
        self.residue_a_dayhoff = residue_a_dayhoff
        self.residue_b_dayhoff = residue_b_dayhoff
        self.residue_c_dayhoff = residue_c_dayhoff
        self.residue_d_dayhoff = residue_d_dayhoff
        self.residue_e_dayhoff = residue_e_dayhoff
        self.residue_f_dayhoff = residue_f_dayhoff
        self.residue_g_dayhoff = residue_g_dayhoff
        self.residue_h_dayhoff = residue_h_dayhoff
        self.residue_i_dayhoff = residue_i_dayhoff
        self.residue_j_dayhoff = residue_j_dayhoff
        self.residue_k_dayhoff = residue_k_dayhoff
        self.residue_l_dayhoff = residue_l_dayhoff
        self.residue_m_dayhoff = residue_m_dayhoff
        self.residue_n_dayhoff = residue_n_dayhoff
        self.residue_o_dayhoff = residue_o_dayhoff
        self.residue_p_dayhoff = residue_p_dayhoff
        self.residue_q_dayhoff = residue_q_dayhoff
        self.residue_r_dayhoff = residue_r_dayhoff
        self.residue_s_dayhoff = residue_s_dayhoff
        self.residue_t_dayhoff = residue_t_dayhoff
        self.residue_u_dayhoff = residue_u_dayhoff
        self.residue_v_dayhoff = residue_v_dayhoff
        self.residue_w_dayhoff = residue_w_dayhoff
        self.residue_x_dayhoff = residue_x_dayhoff
        self.residue_y_dayhoff = residue_y_dayhoff
        self.residue_z_dayhoff = residue_z_dayhoff
        self.property_tiny_number = property_tiny_number
        self.property_small_number = property_small_number
        self.property_aliphatic_number = property_aliphatic_number
        self.property_aromatic_number = property_aromatic_number
        self.property_nonpolar_number = property_nonpolar_number
        self.property_polar_number = property_polar_number
        self.property_charged_number = property_charged_number
        self.property_basic_number = property_basic_number
        self.property_acidic_number = property_acidic_number
        self.property_tiny_mole = property_tiny_mole
        self.property_small_mole = property_small_mole
        self.property_aliphatic_mole = property_aliphatic_mole
        self.property_aromatic_mole = property_aromatic_mole
        self.property_nonpolar_mole = property_nonpolar_mole
        self.property_polar_mole = property_polar_mole
        self.property_charged_mole = property_charged_mole
        self.property_basic_mole = property_basic_mole
        self.property_acidic_mole = property_acidic_mole
        return

    @classmethod
    def from_block(cls, lines: Sequence[str]) -> "PepStats":
        record: Dict[str, str] = dict()

        if not isinstance(lines, Iterable):
            ilines = enumerate(iter(lines))
        else:
            ilines = enumerate(lines)

        parse_regex_line(ilines, NAME_REGEX, record)
        parse_regex_line(ilines, MOLWGT_REGEX, record)
        parse_regex_line(ilines, AVWGT_REGEX, record)
        parse_regex_line(ilines, ISOELEC_REGEX, record)
        parse_regex_line(ilines, A280_MOLAR_REGEX, record)
        parse_regex_line(ilines, A280_1MGML_REGEX, record)
        parse_regex_line(ilines, IMPROB_REGEX, record)

        i, line = get_line(ilines)
        if not line.startswith("Residue"):
            raise BlockParseError(
                i,
                "Expected to get the header line for the Residues table."
            )

        record.update(parse_residue_table(ilines))

        i, line = get_line(ilines)
        if not line.startswith("Property"):
            raise ParseError(
                None,
                i,
                "Expected to get the header line for the property table."
            )

        record.update(parse_property_table(ilines))

        fields = tuple(
            type_(record.get(cname))
            for cname, type_
            in zip(cls.columns, cls.types)
        )

        return cls(*fields)

    @classmethod  # noqa
    def from_file(cls, handle: TextIO) -> Iterator["PepStats"]:
        block: List[str] = []

        for i, line in enumerate(handle):
            sline = line.strip()

            if sline.startswith("PEPSTATS") and len(block) > 0:
                try:
                    yield cls.from_block(block)

                except BlockParseError as e:
                    raise ParseError.from_block_error(
                        e,
                        i - len(block),
                        handle
                    )

                block = [sline]

            elif (len(block) == 0) and (sline == ""):
                continue

            else:
                block.append(sline)

        if len(block) > 0:
            try:
                yield cls.from_block(block)

            except BlockParseError as e:
                raise ParseError.from_block_error(
                    e,
                    i - len(block),
                    handle
                )

        return
