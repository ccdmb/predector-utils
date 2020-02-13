#!/usr/bin/env python3

class Phobius(NamedTuple):

    """ .
    """

    name: str
    tm: int
    sp: bool
    topology: str

    def as_dict(self) -> Dict[str, Union[str, int, float, bool]]:
        return {k: getattr(self, k) for k in self._fields}

    @classmethod
    def from_dict(
        cls,
        d: Dict[str, Union[str, int, float, bool]]
    ) -> "Phobius":
        # assertions appease the typechecker
        name = d["name"]
        assert isinstance(name, str)

        tm = d["tm"]
        assert isinstance(tm, int)

        sp = d["sp"]
        assert isinstance(sp, bool)

        topology = d["topology"]
        assert isinstance(topology, str)
        return cls(name, tm, sp, topology)

    @classmethod
    def from_short_line(cls, line: str) -> "Phobius":
        """ Parse a phobius line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = MULTISPACE_REGEX.split(line.strip())

        if len(sline) != 4:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 4 but got {len(sline)}"
            )

        if sline == ["SEQENCE", "ID", "TM", "SP", "PREDICTION"]:
            raise LineParseError("The line appears to be the header line")

        return cls(
            parse_string_not_empty(sline[0], "name"),
            parse_int(sline[1], "tm"),
            parse_bool(sline[2], "sp", "Y", "0"),
            parse_string_not_empty(sline[3], "topology")
        )

    @classmethod
    def from_short_file(cls, handle: TextIO) -> Iterator["Phobius"]:
        for i, line in enumerate(handle):
            sline = line.strip()
            if sline.startswith("#"):
                continue
            elif i == 0 and sline.startswith("SEQENCE"):
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

