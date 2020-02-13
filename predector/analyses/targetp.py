#!/usr/bin/env python3

class TargetP(NamedTuple):

    """ Doesn't have output format documentation yet
    """

    name: str
    prediction: str
    other: float
    sp: float
    mtp: float
    ctp: Optional[float]
    lutp: Optional[float]
    cs_pos: Optional[str]

    def as_dict(self) -> Dict[str, Union[str, int, float, bool]]:
        return {k: getattr(self, k) for k in self._fields}

    @classmethod
    def from_dict(
        cls,
        d: Dict[str, Union[str, int, float, bool]]
    ) -> "TargetP":
        # assertions appease the typechecker
        name = d["name"]
        assert isinstance(name, str)

        prediction = d["prediction"]
        assert isinstance(prediction, str)

        other = d["other"]
        assert isinstance(other, float)

        sp = d["sp"]
        assert isinstance(sp, float)

        mtp = d["mtp"]
        assert isinstance(mtp, float)

        ctp = d.get("ctp", None)
        assert isinstance(ctp, float) or ctp is None

        lutp = d.get("lutp", None)
        assert isinstance(lutp, float) or lutp is None

        cs_pos = d.get("cs_pos", None)
        assert isinstance(cs_pos, str) or cs_pos is None
        return cls(name, prediction, other, sp, mtp, ctp, lutp, cs_pos)

    @classmethod
    def from_short_line_plant(cls, line: str) -> "TargetP":
        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) == 8:
            cs_pos: Optional[str] = str(sline[7])
        elif len(sline) == 7:
            cs_pos = None
        else:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 7 or 8 but got {len(sline)}"
            )

        return cls(
            parse_string_not_empty(sline[0], "name"),
            is_one_of(
                sline[1],
                ["OTHER", "SP", "mTP", "cTP", "luTP"],
                "prediction"
            ),
            parse_float(sline[2], "OTHER"),
            parse_float(sline[3], "SP"),
            parse_float(sline[4], "mTP"),
            parse_float(sline[5], "cTP"),
            parse_float(sline[6], "luTP"),
            cs_pos,
        )

    @classmethod
    def from_short_line_non_plant(cls, line: str) -> "TargetP":
        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) == 6:
            cs_pos: Optional[str] = str(sline[5])
        elif len(sline) == 5:
            cs_pos = None
        else:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 5 or 6 but got {len(sline)}"
            )

        prediction = is_one_of(
            sline[1],
            ["noTP", "SP", "mTP"],
            "prediction"
        )

        if prediction == "noTP":
            prediction = "OTHER"

        return cls(
            parse_string_not_empty(sline[0], "name"),
            prediction,
            parse_float(sline[2], "OTHER"),
            parse_float(sline[3], "SP"),
            parse_float(sline[4], "mTP"),
            ctp=None,
            lutp=None,
            cs_pos=cs_pos,
        )

    @classmethod
    def from_short_file(
        cls,
        handle: TextIO,
        plant: bool = False
    ) -> Iterator["TargetP"]:
        for i, line in enumerate(handle):
            sline = line.strip()
            if sline.startswith("#"):
                continue
            elif sline == "":
                continue

            try:
                if plant:
                    yield cls.from_short_line_plant(sline)
                else:
                    yield cls.from_short_line_non_plant(sline)

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

