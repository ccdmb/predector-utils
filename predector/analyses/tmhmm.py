#!/usr/bin/env python3

class TMHMM(NamedTuple):

    """ .
    """

    name: str
    length: int
    exp_aa: float
    first_60: float
    pred_hel: int
    topology: str

    def as_dict(self) -> Dict[str, Union[str, int, float, bool]]:
        return {k: getattr(self, k) for k in self._fields}

    @classmethod
    def from_dict(cls, d: Dict[str, Union[str, int, float, bool]]) -> "TMHMM":
        # assertions appease the typechecker
        name = d["name"]
        assert isinstance(name, str)

        length = d["length"]
        assert isinstance(length, int)

        exp_aa = d["exp_aa"]
        assert isinstance(exp_aa, float)

        first_60 = d["first_60"]
        assert isinstance(first_60, float)

        pred_hel = d["pred_hel"]
        assert isinstance(pred_hel, int)

        topology = d["topology"]
        assert isinstance(topology, str)
        return cls(name, length, exp_aa, first_60, pred_hel, topology)

    @classmethod
    def from_short_line(cls, line: str) -> "TMHMM":
        """ Parse a tmhmm line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) != 6:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 6 but got {len(sline)}"
            )

        return cls(
            parse_string_not_empty(sline[0], "name"),
            parse_int(split_at_eq(sline[1], "length", "len"), "length"),
            parse_float(split_at_eq(sline[2], "exp_aa", "ExpAA"), "exp_aa"),
            parse_float(
                split_at_eq(sline[3], "first_60", "First60"),
                "first_60"
            ),
            parse_int(
                split_at_eq(sline[4], "pred_hel", "PredHel"),
                "pred_hel"
            ),
            parse_string_not_empty(
                split_at_eq(sline[5], "topology", "Topology"),
                "topology"
            )
        )

    @classmethod
    def from_short_file(cls, handle: TextIO) -> Iterator["TMHMM"]:
        for i, line in enumerate(handle):
            sline = line.strip()
            if sline.startswith("#"):
                continue
            elif sline == "":
                continue

            try:
                yield cls.from_short_line(sline)

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

