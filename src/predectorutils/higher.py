""" Some higher order functions for dealing with static analysis. """

from typing import TypeVar
from typing import Callable


T = TypeVar("T")
U = TypeVar("U")


def fmap(function: Callable[[T], U], option: T | None) -> U | None:
    if option is None:
        return None
    else:
        return function(option)


def applicative(
    function: Callable[[T], U | None],
    option: T | None,
) -> U | None:
    """ Same as fmap except for the type signature. """
    if option is None:
        return None
    else:
        return function(option)
    return


def or_else(default: U, option: T | None) -> T | U:
    """ Replaces None with some default value. """
    if option is None:
        return default
    else:
        return option
