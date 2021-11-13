#!/usr/bin/env python3

import sys
import traceback
import argparse

from typing import List

from predectorutils import __version__, __email__
from predectorutils.subcommands.r2js import cli as r2js_cli
from predectorutils.subcommands.r2js import runner as r2js_runner

from predectorutils.subcommands.encode import cli as encode_cli
from predectorutils.subcommands.encode import runner as encode_runner

from predectorutils.subcommands.split_fasta import cli as split_fasta_cli
from predectorutils.subcommands.split_fasta import runner as split_fasta_runner

from predectorutils.subcommands.analysis_tables import cli as tables_cli
from predectorutils.subcommands.analysis_tables import runner as tables_runner

from predectorutils.subcommands.gff import cli as gff_cli
from predectorutils.subcommands.gff import runner as gff_runner

from predectorutils.subcommands.rank import cli as rank_cli
from predectorutils.subcommands.rank import runner as rank_runner

from predectorutils.subcommands.decode import cli as decode_cli
from predectorutils.subcommands.decode import runner as decode_runner

from predectorutils.subcommands.regex import cli as regex_cli
from predectorutils.subcommands.regex import runner as regex_runner

from predectorutils.subcommands.precomputed import cli as precomp_cli
from predectorutils.subcommands.precomputed import runner as precomp_runner

from predectorutils.subcommands.load_db import cli as load_db_cli
from predectorutils.subcommands.load_db import runner as load_db_runner

from predectorutils.subcommands.dump_db import cli as dump_db_cli
from predectorutils.subcommands.dump_db import runner as dump_db_runner

from predectorutils.parsers import ParseError
from predectorutils.exceptions import (
    EXIT_VALID, EXIT_KEYBOARD, EXIT_UNKNOWN, EXIT_CLI, EXIT_INPUT_FORMAT,
    EXIT_INPUT_NOT_FOUND, EXIT_SYSERR, EXIT_CANT_OUTPUT
)


class MyArgumentParser(argparse.ArgumentParser):

    def error(self, message: str):
        """ Override default to have more informative exit codes. """
        self.print_usage(sys.stderr)
        raise MyArgumentError("{}: error: {}".format(self.prog, message))


class MyArgumentError(Exception):

    def __init__(self, message: str):
        self.message = message
        self.errno = EXIT_CLI

        # This is a bit hacky, but I can't figure out another way to do it.
        if "No such file or directory" in message:
            if "infile" in message:
                self.errno = EXIT_INPUT_NOT_FOUND
            elif "outfile" in message:
                self.errno = EXIT_CANT_OUTPUT
        return


def cli(prog: str, args: List[str]) -> argparse.Namespace:
    parser = MyArgumentParser(
        prog=prog,
        description=(
            "Examples:"
        ),
        epilog=(
            "Exit codes:\n\n"
            f"{EXIT_VALID} - Everything's fine\n"
            f"{EXIT_KEYBOARD} - Keyboard interrupt\n"
            f"{EXIT_CLI} - Invalid command line usage\n"
            f"{EXIT_INPUT_FORMAT} - Input format error\n"
            f"{EXIT_INPUT_NOT_FOUND} - Cannot open the input\n"
            f"{EXIT_SYSERR} - System error\n"
            f"{EXIT_CANT_OUTPUT} - Can't create output file\n"
            f"{EXIT_UNKNOWN} - Unhandled exception, please file a bug!\n"
        )
    )

    parser.add_argument(
        '--version',
        action='version',
        version='{}'.format(__version__),
        help="Print the version of %(prog)s and exit"
    )

    subparsers = parser.add_subparsers(dest="subparser_name")
    r2js_subparser = subparsers.add_parser(
        "r2js",
        help=(
            "Parse the result of some analysis as a common line-delimited "
            "json format."
        )
    )

    r2js_cli(r2js_subparser)

    encode_subparser = subparsers.add_parser(
        "encode",
        help=(
            "Remove duplicate sequences and give them new names."
        )
    )

    encode_cli(encode_subparser)

    split_fasta_subparser = subparsers.add_parser(
        "split_fasta",
        help=(
            "Split a fasta file into chunks."
        )
    )

    split_fasta_cli(split_fasta_subparser)

    tables_subparser = subparsers.add_parser(
        "tables",
        help=(
            "Split line-delimited json into tsvs for each analysis."
        )
    )

    tables_cli(tables_subparser)

    gff_subparser = subparsers.add_parser(
        "gff",
        help=(
            "Write analyses with location information as a GFF3 file."
        )
    )

    gff_cli(gff_subparser)

    rank_subparser = subparsers.add_parser(
        "rank",
        help=(
            "Rank effector candidates."
        )
    )

    rank_cli(rank_subparser)

    decode_subparser = subparsers.add_parser(
        "decode",
        help=(
            "Decode and reduplicate encoded names."
        )
    )

    decode_cli(decode_subparser)

    regex_subparser = subparsers.add_parser(
        "regex",
        help=(
            "Search for a regular expression in sequences."
        )
    )

    regex_cli(regex_subparser)

    precomp_subparser = subparsers.add_parser(
        "precomputed",
        help=(
            "Fetch precomputed results and output remaining sequences."
        )
    )

    precomp_cli(precomp_subparser)

    load_db_subparser = subparsers.add_parser(
        "load_db",
        help=(
            "Load ldjson results into an SQlite database."
        )
    )

    load_db_cli(load_db_subparser)

    dump_db_subparser = subparsers.add_parser(
        "dump_db",
        help=(
            "Dump an SQLite database into ldjson results."
        )
    )

    dump_db_cli(dump_db_subparser)

    parsed = parser.parse_args(args)

    if parsed.subparser_name is None:
        parser.print_help()
        sys.exit(0)

    return parsed


def main():  # noqa
    try:
        args = cli(prog=sys.argv[0], args=sys.argv[1:])
    except MyArgumentError as e:
        print(e.message, file=sys.stderr)
        sys.exit(e.errno)

    try:
        if args.subparser_name == "r2js":
            r2js_runner(args)
        elif args.subparser_name == "encode":
            encode_runner(args)
        elif args.subparser_name == "split_fasta":
            split_fasta_runner(args)
        elif args.subparser_name == "tables":
            tables_runner(args)
        elif args.subparser_name == "gff":
            gff_runner(args)
        elif args.subparser_name == "rank":
            rank_runner(args)
        elif args.subparser_name == "decode":
            decode_runner(args)
        elif args.subparser_name == "regex":
            regex_runner(args)
        elif args.subparser_name == "precomputed":
            precomp_runner(args)
        elif args.subparser_name == "load_db":
            load_db_runner(args)
        elif args.subparser_name == "dump_db":
            dump_db_runner(args)

    except ParseError as e:
        if e.line is not None:
            header = "Failed to parse file <{}> at line {}.\n".format(
                e.filename, e.line)
        else:
            header = "Failed to parse file <{}>.\n".format(e.filename)

        print("{}\n{}".format(header, e.message), file=sys.stderr)
        sys.exit(EXIT_INPUT_FORMAT)

    except OSError as e:
        msg = (
            "Encountered a system error.\n"
            "We can't control these, and they're usually related to your OS.\n"
            "Try running again.\n"
        )
        raise e
        print(msg, file=sys.stderr)
        print(e.strerror, file=sys.stderr)
        sys.exit(EXIT_SYSERR)

    except MemoryError:
        msg = (
            "Ran out of memory!\n"
            "Catastrophy shouldn't use much RAM, so check other "
            "processes and try running again."
        )
        print(msg, file=sys.stderr)
        sys.exit(EXIT_SYSERR)

    except KeyboardInterrupt:
        print("Received keyboard interrupt. Exiting.", file=sys.stderr)
        sys.exit(EXIT_KEYBOARD)

    except Exception as e:
        msg = (
            "I'm so sorry, but we've encountered an unexpected error.\n"
            "This shouldn't happen, so please file a bug report with the "
            "authors.\nWe will be extremely grateful!\n\n"
            "You can email us at {}.\n"
            "Alternatively, you can file the issue directly on the repo "
            "<https://github.com/ccdmb/predector-utils>\n\n"
            "Please attach a copy of the following message:"
        ).format(__email__)
        print(e, file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(EXIT_UNKNOWN)

    return


if __name__ == '__main__':
    main()
