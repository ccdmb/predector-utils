#!/usr/bin/env python3

import re
from typing import Optional
from typing import Sequence
from typing import Dict

MULTISPACE_REGEX = re.compile(r"\s+")


class ParseError(Exception):
    """ Some aspect of parsing failed. """

    def __init__(
        self,
        filename: Optional[str],
        line: Optional[int],
        message: str
    ):
        self.filename = filename
        self.line = line
        self.message = message
        return


class LineParseError(Exception):

    def __init__(self, message: str):
        self.message = message
        return


def parse_regex(field: str, regex: re.Pattern) -> Dict[str, str]:
    matches = regex.match(field)
    if matches is None:
        raise LineParseError(
            f"Expected a line matching the regular expression: {regex}, "
            f"but got '{field}'.")
    else:
        return matches.groupdict()


def split_at_eq(field: str, field_name: str, expected_lhs: str) -> str:
    """ Parse a field of a key=value type, returning the value. """
    sfield = field.split("=", maxsplit=1)
    if len(sfield) != 2:
        raise LineParseError(
            f"Expected a '{expected_lhs}=value' field value in column "
            f"{field_name} but got {field}."
        )
    if sfield[0] != expected_lhs:
        raise LineParseError(
            f"Expected a '{expected_lhs}=value' field value in column "
            f"{field_name} but got {field}."
        )

    return sfield[1]


def parse_int(field: str, field_name: str) -> int:
    """ Parse a string as an integer, raising a custom error if fails. """

    try:
        return int(field)
    except ValueError:
        raise LineParseError(
            f"Could not parse value in {field_name} column as an integer. "
            f"The offending value was: '{field}'."
        )


def parse_float(field: str, field_name: str) -> float:
    """ Parse a string as a float, raising a custom error if fails. """

    try:
        return float(field)
    except ValueError:
        raise LineParseError(
            f"Could not parse value in {field_name} column as a float. "
            f"The offending value was: '{field}'."
        )


def parse_bool(
    field: str,
    field_name: str,
    true_value: str,
    false_value: str
) -> bool:
    """ """

    if field == true_value:
        return True
    elif field == false_value:
        return False
    else:
        raise LineParseError(
            f"Invalid value: '{field}' in the column: '{field_name}'. "
            f"Must be either {true_value} or {false_value}."
        )


def is_value(field: str, field_name: str, value: str) -> bool:
    """ """
    return field == value


def parse_string_not_empty(field: str, field_name: str) -> str:
    """ """

    if field.strip() == "":
        raise LineParseError(f"The value in column: '{field_name}' was empty.")
    else:
        return field


def is_one_of(field: str, options: Sequence[str], field_name: str) -> str:
    """ """

    if field in options:
        return field
    else:
        raise LineParseError(
            f"Invalid value: '{field}' in the column: '{field_name}'. "
            f"Must be one of {options}."
        )
