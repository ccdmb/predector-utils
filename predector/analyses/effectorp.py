#!/usr/bin/env python3


class EffectorP1(NamedTuple):

    """ """

    name: str
    prediction: str
    prob: float

    def as_dict(self) -> Dict[str, Union[str, int, float, bool]]:
        return {k: getattr(self, k) for k in self._fields}

    @classmethod
    def from_dict(
        cls,
        d: Dict[str, Union[str, int, float, bool]]
    ) -> "EffectorP1":
        # assertions appease the typechecker
        name = d["name"]
        assert isinstance(name, str)

        prediction = d["prediction"]
        assert isinstance(prediction, str)

        prob = d["prob"]
        assert isinstance(prob, float)
        return cls(name, prediction, prob)

    @classmethod
    def from_line(cls, line: str) -> "EffectorP1":
        """ Parse an EffectorP1 line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) != 3:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 3 but got {len(sline)}."
            )

        return cls(
            parse_string_not_empty(sline[0], "name"),
            is_one_of(
                sline[1],
                ["Effector", "Non-effector"],
                "prediction"
            ),
            parse_float(sline[2], "prob"),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["EffectorP1"]:
        comment = False
        for i, line in enumerate(handle):
            sline = line.strip()
            if comment and line.startswith("---------"):
                comment = False
                continue
            elif comment and line.startswith("# Identifier"):
                comment = False
                continue
            elif comment:
                continue
            elif (i == 0) and line.startswith("---------"):
                comment = True
                continue

            if sline.startswith("#"):
                continue
            elif sline == "":
                continue

            try:
                yield cls.from_line(sline)

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


class EffectorP2(NamedTuple):

    """ """

    name: str
    prediction: str
    prob: float

    def as_dict(self) -> Dict[str, Union[str, int, float, bool]]:
        return {k: getattr(self, k) for k in self._fields}

    @classmethod
    def from_dict(
        cls,
        d: Dict[str, Union[str, int, float, bool]]
    ) -> "EffectorP2":
        # assertions appease the typechecker
        name = d["name"]
        assert isinstance(name, str)

        prediction = d["prediction"]
        assert isinstance(prediction, str)

        prob = d["prob"]
        assert isinstance(prob, float)
        return cls(name, prediction, prob)

    @classmethod
    def from_line(cls, line: str) -> "EffectorP2":
        """ Parse an EffectorP2 line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) != 3:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 3 but got {len(sline)}."
            )

        return cls(
            parse_string_not_empty(sline[0], "name"),
            is_one_of(
                sline[1],
                ["Effector", "Unlikely effector", "Non-effector"],
                "prediction"
            ),
            parse_float(sline[2], "prob"),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["EffectorP2"]:
        comment = False
        for i, line in enumerate(handle):
            sline = line.strip()
            if comment and line.startswith("---------"):
                comment = False
                continue
            elif comment and line.startswith("# Identifier"):
                comment = False
                continue
            elif comment:
                continue
            elif (i == 0) and line.startswith("---------"):
                comment = True
                continue

            if sline.startswith("#"):
                continue
            elif sline == "":
                continue

            try:
                yield cls.from_line(sline)

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
class EffectorP1(NamedTuple):

    """ """

    name: str
    prediction: str
    prob: float

    def as_dict(self) -> Dict[str, Union[str, int, float, bool]]:
        return {k: getattr(self, k) for k in self._fields}

    @classmethod
    def from_dict(
        cls,
        d: Dict[str, Union[str, int, float, bool]]
    ) -> "EffectorP1":
        # assertions appease the typechecker
        name = d["name"]
        assert isinstance(name, str)

        prediction = d["prediction"]
        assert isinstance(prediction, str)

        prob = d["prob"]
        assert isinstance(prob, float)
        return cls(name, prediction, prob)

    @classmethod
    def from_line(cls, line: str) -> "EffectorP1":
        """ Parse an EffectorP1 line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) != 3:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 3 but got {len(sline)}."
            )

        return cls(
            parse_string_not_empty(sline[0], "name"),
            is_one_of(
                sline[1],
                ["Effector", "Non-effector"],
                "prediction"
            ),
            parse_float(sline[2], "prob"),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["EffectorP1"]:
        comment = False
        for i, line in enumerate(handle):
            sline = line.strip()
            if comment and line.startswith("---------"):
                comment = False
                continue
            elif comment and line.startswith("# Identifier"):
                comment = False
                continue
            elif comment:
                continue
            elif (i == 0) and line.startswith("---------"):
                comment = True
                continue

            if sline.startswith("#"):
                continue
            elif sline == "":
                continue

            try:
                yield cls.from_line(sline)

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


class EffectorP2(NamedTuple):

    """ """

    name: str
    prediction: str
    prob: float

    def as_dict(self) -> Dict[str, Union[str, int, float, bool]]:
        return {k: getattr(self, k) for k in self._fields}

    @classmethod
    def from_dict(
        cls,
        d: Dict[str, Union[str, int, float, bool]]
    ) -> "EffectorP2":
        # assertions appease the typechecker
        name = d["name"]
        assert isinstance(name, str)

        prediction = d["prediction"]
        assert isinstance(prediction, str)

        prob = d["prob"]
        assert isinstance(prob, float)
        return cls(name, prediction, prob)

    @classmethod
    def from_line(cls, line: str) -> "EffectorP2":
        """ Parse an EffectorP2 line as an object. """

        if line == "":
            raise LineParseError("The line was empty.")

        sline = line.strip().split("\t")

        if len(sline) != 3:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 3 but got {len(sline)}."
            )

        return cls(
            parse_string_not_empty(sline[0], "name"),
            is_one_of(
                sline[1],
                ["Effector", "Unlikely effector", "Non-effector"],
                "prediction"
            ),
            parse_float(sline[2], "prob"),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["EffectorP2"]:
        comment = False
        for i, line in enumerate(handle):
            sline = line.strip()
            if comment and line.startswith("---------"):
                comment = False
                continue
            elif comment and line.startswith("# Identifier"):
                comment = False
                continue
            elif comment:
                continue
            elif (i == 0) and line.startswith("---------"):
                comment = True
                continue

            if sline.startswith("#"):
                continue
            elif sline == "":
                continue

            try:
                yield cls.from_line(sline)

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

