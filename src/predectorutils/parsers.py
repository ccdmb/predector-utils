#!/usr/bin/env python3

import re

from collections.abc import Sequence, Mapping
from typing import TypeVar
from typing import TextIO
from typing import Pattern
from typing import Callable

from .higher import or_else

MULTISPACE_REGEX = re.compile(r"\s+")
T = TypeVar("T")


class ValueParseError(Exception):

    def __init__(
        self,
        got: str,
        expected: str,
        template: str = "Could not parse value {got} as {expected}."
    ):
        self.got = got
        self.expected = expected
        self.template = template

    def __str__(self) -> str:
        return self.template.format(got=self.got, expected=self.expected)

    def as_field_error(
        self,
        field: str,
        kind: str = "field"
    ) -> "FieldParseError":
        return FieldParseError(field, str(self), kind=kind)

    def as_line_error(
        self,
    ) -> "LineParseError":
        return LineParseError(str(self))

    def as_block_error(self, line: int | None = None) -> "BlockParseError":
        return BlockParseError(line, str(self))

    def as_parse_error(
        self,
        filename: str | None = None,
        line: int | None = None
    ) -> "ParseError":
        return ParseError(filename, line, str(self))


class FieldParseError(Exception):

    def __init__(self, field: str, message: str, kind: str = "column"):
        self.field = field
        self.message = message
        self.kind = kind
        return

    def __str__(self) -> str:
        return f"In {self.kind} '{self.field}': {self.message}"

    def as_line_error(self) -> "LineParseError":
        return LineParseError(str(self))

    def as_block_error(self, line: int | None = None) -> "BlockParseError":
        return BlockParseError(line, str(self))

    def as_parse_error(
        self,
        filename: str | None = None,
        line: int | None = None,
    ) -> "ParseError":
        return ParseError(filename, line, str(self))


class LineParseError(Exception):
    def __init__(self, message: str):
        self.message = message
        return

    def __str__(self) -> str:
        return self.message

    def as_block_error(self, line: int | None = None) -> "BlockParseError":
        return BlockParseError(line, str(self))

    def as_parse_error(
        self,
        filename: str | None = None,
        line: int | None = None,
    ) -> "ParseError":
        return ParseError(filename, line, str(self))


class BlockParseError(Exception):

    def __init__(self, line: int | None, message: str):
        self.line = line
        self.message = message
        return

    def as_block_error(self, line: int | None = None) -> "BlockParseError":
        if (line is None) and (self.line is None):
            nline = None
        else:
            nline = or_else(0, line) + or_else(0, self.line)

        return BlockParseError(
            nline,
            self.message
        )

    def as_parse_error(
        self,
        filename: str | None = None,
        line: int | None = None,
    ) -> "ParseError":
        if (line is None) and (self.line is None):
            line = None
        else:
            line = or_else(0, line) + or_else(0, self.line)

        return ParseError(
            filename,
            line,
            self.message
        )


class ParseError(Exception):
    """ Some aspect of parsing failed. """

    def __init__(
        self,
        filename: str | None,
        line: int | None,
        message: str
    ):
        self.filename = filename
        self.line = line
        self.message = message
        return

    def add_filename_from_handle(self, handle: TextIO) -> "ParseError":
        if hasattr(handle, "name"):
            self.filename = handle.name
        else:
            self.filename = None

        return self


def convert_line_err(
    lineno: int,
    field: str,
    parser: Callable[[str], T]
) -> T:
    try:
        return parser(field)
    except LineParseError as e:
        raise e.as_block_error(lineno)


def split_at_eq(
    fn: Callable[[str], T],
    expected_lhs: str
) -> Callable[[str], T | ValueParseError]:

    def inner(value: str) -> T | ValueParseError:
        """ Parse a field of a "key value" type, returning the value. """
        sfield = value.split("=", maxsplit=1)
        if len(sfield) != 2:
            return ValueParseError(
                value,
                f"{expected_lhs}=<value>",
                "Expected a string of the form '{expected}', but got '{got}'",
            )
        elif sfield[0] != expected_lhs:
            return ValueParseError(
                value,
                f"{expected_lhs} <value>",
                "Expected a string of the form '{expected}', but got '{got}'",
            )
        else:
            return fn(sfield[1])

    return inner


def split_at_multispace(
    fn: Callable[[str], T],
    expected_lhs: str
) -> Callable[[str], T | ValueParseError]:

    def inner(value: str) -> T | ValueParseError:
        """ Parse a field of a "key value" type, returning the value. """
        sfield = MULTISPACE_REGEX.split(value, maxsplit=1)
        if len(sfield) != 2:
            return ValueParseError(
                value,
                f"{expected_lhs} <value>",
                "Expected a string of the form '{expected}', but got '{got}'",
            )
        elif sfield[0] != expected_lhs:
            return ValueParseError(
                value,
                f"{expected_lhs} <value>",
                "Expected a string of the form '{expected}', but got '{got}'",
            )
        else:
            return fn(sfield[1])

    return inner


def is_value(expected: str) -> Callable[[str], bool]:
    """ """

    def inner(value: str) -> bool:
        return value == expected

    return inner


def is_not_empty_string(value: str) -> str:
    if value.strip() == "":
        raise ValueError("String cannot be empty")
    else:
        return value


def is_one_of(
    options: Sequence[str]
) -> Callable[[str], ValueParseError | str]:

    def inner(value: str) -> ValueParseError | str:
        if value in options:
            return value
        else:
            return ValueParseError(
                value,
                "{" + ", ".join(options) + "}",
                "Invalid value '{got}', must be one of {expected}."
            )

    return inner


def parse_bool(
    true_value: str,
    false_value: str
) -> Callable[[str], ValueParseError | bool]:
    """ """

    def inner(value: str) -> ValueParseError | bool:
        if value == true_value:
            return True
        elif value == false_value:
            return False
        else:
            return ValueParseError(
                value,
                f"{true_value} or {false_value}",
                "Invalid value '{got}', must be either {expected}."
            )

    return inner


def parse_bool_options(
    true_values: Sequence[str],
    false_values: Sequence[str]
) -> Callable[[str], ValueParseError | bool]:
    """ """
    trues = set(true_values)
    falses = set(false_values)

    def inner(value: str) -> ValueParseError | bool:
        if value in trues:
            return True
        elif value in falses:
            return False
        else:
            return ValueParseError(
                value,
                f"{repr(true_values)} or {repr(false_values)}",
                "Invalid value '{got}', must be one of {expected}."
            )

    return inner


def parse_regex(
    regex: Pattern
) -> Callable[[str], ValueParseError | dict[str, str]]:

    def inner(value: str) -> ValueParseError | dict[str, str]:
        matches = regex.match(value)
        if matches is None:
            return ValueParseError(
                value,
                str(regex),
                ("Value '{got}' does not match the "
                 "regular expression {expected}.")
            )
        else:
            return matches.groupdict()

    return inner


def parse_or_none(
    fn: Callable[[str], ValueParseError | T],
    none_value: str,
) -> Callable[[str], ValueParseError | T | None]:
    """ If the value is the same as the none value, will return None.
    Otherwise will attempt to run the fn with field and field name as the
    first and 2nd arguments.
    """

    def inner(value: str) -> ValueParseError | T | None:
        if value == none_value:
            return None

        val = fn(value)
        if isinstance(val, ValueParseError):
            val.template = (
                f"{val.template}. The value may also be "
                f"'{none_value}', which will be interpreted as None."
            )

        return val

    return inner


def parse_delim(
    delim: str,
    fn: Callable[[str], ValueParseError | T],
    ignore_empty: bool = True,
) -> Callable[[str], ValueParseError | list[T]]:
    """ Split the string and apply a function to each element.
    """
    import re
    splitter = re.compile(delim)

    def inner(value: str) -> ValueParseError | list[T]:
        svalue = splitter.split(value)
        out: list[T] = []
        errors: list[ValueParseError] = []
        for si in svalue:
            if ignore_empty and (si == ''):
                continue

            sj = fn(si)

            if isinstance(sj, ValueParseError):
                errors.append(sj)
            else:
                out.append(sj)

        if len(errors) == 1:
            return errors[0]
        elif len(errors) > 1:
            return ValueParseError(
                got="\n".join(map(str, errors)),
                expected="",
                template="\n{got}"
            )
        return out

    return inner


def parse_sequence(
    options: Sequence[str],
) -> Callable[[str], ValueParseError | str]:
    """ Convenience method to check that all characters in a string are valid.
    """

    checker = is_one_of(options)

    def inner(value: str) -> ValueParseError | str:
        errors: list[ValueParseError] = []
        for char in value:

            check = checker(char)

            if isinstance(check, ValueParseError):
                errors.append(check)

        if len(errors) == 1:
            return errors[0]
        elif len(errors) > 1:
            return ValueParseError(
                got="\n".join(map(str, errors)),
                expected="",
                template="\n{got}"
            )

        return value

    return inner



def get_from_dict_or_err(
    key: str
) -> Callable[[Mapping[str, T]], ValueParseError | T]:

    def inner(d: Mapping[str, T]) -> ValueParseError | T:
        if key in d:
            return d[key]
        else:
            return ValueParseError(key, "", "{got} was missing.")

    return inner


def parse_value(
    fn: Callable[[str], T],
    expected: str,
    etemplate: str | None = None,
) -> Callable[[str], ValueParseError | T]:

    def inner(value: str) -> ValueParseError | T:
        try:
            return fn(value)
        except ValueParseError as e:
            return e
        except ValueError:
            if etemplate is not None:
                return ValueParseError(value, expected, etemplate)
            else:
                return ValueParseError(value, expected)

    return inner




parse_int = parse_value(int, "an integer")
parse_float = parse_value(float, "a float")
parse_str = parse_value(is_not_empty_string, "a non-empty string")


def parse_field(
    fn: Callable[[str], ValueParseError | T],
    field: str,
    kind: str = "field",
) -> Callable[[str], FieldParseError | T]:

    def inner(value: str) -> FieldParseError | T:
        result = fn(value)
        if isinstance(result, ValueParseError):
            return result.as_field_error(field, kind=kind)
        else:
            return result

    return inner


def raise_it(fn: Callable[[str], Exception | T]) -> Callable[[str], T]:

    def inner(value: str) -> T:
        fval = fn(value)
        if isinstance(fval, Exception):
            raise fval
        else:
            return fval

    return inner
