#!/usr/bin/env python3

import sys
import argparse
import enum
import datetime


class ValidResults(enum.Enum):

    signalp3_nn = 1
    signalp3_hmm = 2
    signalp4 = 3
    signalp5 = 4
    deepsig = 5
    phobius = 6
    tmhmm = 7
    deeploc = 8
    targetp = 9
    effectorp1 = 10
    effectorp2 = 11
    apoplastp = 12
    localizer = 13
    pfamscan = 14
    hmmer_domtbl = 15
    mmseqs = 16

    def __str__(self) -> str:
        return self.name

    @classmethod
    def from_string(cls, s: str) -> "ValidResults":
        try:
            return cls[s]
        except KeyError:
            raise ValueError(f"{s} is not a valid result type to parse.")


def cli(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "format",
        dest="file_format",
        type=ValidResults.from_string,
        choices=list(ValidResults),
        help="The file results to parse into a line delimited JSON format."
    )

    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help="The text file to parse as input. Use '-' for stdin."
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Where to write the output to. Default: stdout"
    )

    parser.add_argument(
        "-r", "--run-name",
        dest="run_name",
        type=str,
        default=None,
        help="Add a run name to the output."
    )

    parser.add_argument(
        "-s", "--session-id",
        dest="session_id",
        type=str,
        default=None,
        help="Add a session id to the output."
    )

    # Should I convert this to a proper time object?
    parser.add_argument(
        "-t", "--start",
        dest="start",
        type=datetime.fromisoformat,
        default=None,
        help=("Add the run start-time to the output. Should be in the "
              "iso format i.e. '2011-11-04 00:05:23.283'")
    )

    return


def runner(args: argparse.Namespace) -> None:
    return
