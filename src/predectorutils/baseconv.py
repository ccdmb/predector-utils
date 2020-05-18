# Copyright (c) 2010-2017 Guilherme Gondim. All rights reserved.
# Copyright (c) 2009 Simon Willison. All rights reserved.
# Copyright (c) 2002 Drew Perttula. All rights reserved.
#
# License:
#   Python Software Foundation License version 2
#
# See the file "LICENSE" for terms & conditions for usage, and a DISCLAIMER OF
# ALL WARRANTIES.
#
# This Baseconv distribution contains no GNU General Public Licensed (GPLed)
# code so it may be used in proprietary projects just like prior ``baseconv``
# distributions.
#
# All trademarks referenced herein are property of their respective holders.
#
# This code except IdConverter is from
# <https://github.com/semente/python-baseconv>
# __version__ = '1.2.2'


__all__ = ["IdConverter"]


class BaseConverter(object):
    decimal_digits = '0123456789'

    def __init__(self, digits, sign='-'):
        self.sign = sign
        self.digits = digits
        if sign in self.digits:
            raise ValueError('sign character found in converter base digits')
        if len(self.digits) <= 1:
            raise ValueError('converter base digits length too short')

    def __repr__(self):
        data = (self.__class__.__name__, self.digits, self.sign)
        return "%s(%r, sign=%r)" % data

    def _convert(self, number, from_digits, to_digits):
        # make an integer out of the number
        x = 0
        for digit in str(number):
            try:
                x = x * len(from_digits) + from_digits.index(digit)
            except ValueError:
                raise ValueError('invalid digit "%s"' % digit)

        # create the result in base 'len(to_digits)'
        if x == 0:
            res = to_digits[0]
        else:
            res = ''
            while x > 0:
                digit = x % len(to_digits)
                res = to_digits[digit] + res
                x = int(x // len(to_digits))
        return res

    def encode(self, number):
        if str(number)[0] == '-':
            neg = True
            number = str(number)[1:]
        else:
            neg = False

        value = self._convert(number, self.decimal_digits, self.digits)
        if neg:
            return self.sign + value
        return value

    def decode(self, number):
        if str(number)[0] == self.sign:
            neg = True
            number = str(number)[1:]
        else:
            neg = False

        value = self._convert(number, self.digits, self.decimal_digits)
        if neg:
            return '-' + value
        return value


class IdConverter(object):

    """
    Generates an id based on an integer.
    Essentially it's a small hashing algorithm that allows us to
    give unique ids that aren't too long.
    """

    def __init__(
            self,
            state: int = 0,
            prefix: str = "",
            length: int = 4,
            alphabet: str = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ',
    ) -> None:
        """ Create a new id converter with a given state and alphabet. """

        self.state = state
        self.prefix = prefix
        self.length = length
        self.alphabet = alphabet

        self.converter = BaseConverter(alphabet)
        return

    def encode(self, number: int) -> str:
        """ Given an integer get the string representation. """

        template = "{pre}{{p:{first}>{length}}}".format(
            pre=self.prefix,
            first=self.converter.digits[0],
            length=self.length,
            )
        return template.format(p=self.converter.encode(number))

    def decode(self, pattern: str) -> int:
        """ Given a string representation get the integer that produced it. """

        pattern = pattern[len(self.prefix):]
        return int(self.converter.decode(pattern))

    def __next__(self):
        string = self.encode(self.state)
        self.state += 1
        return string

    def __iter__(self):
        return self
