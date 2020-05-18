from predectorutils.baseconv import IdConverter


def test_encodes():
    id_conv = IdConverter(prefix="TS", length=4)
    assert id_conv.encode(0) == "TS0000"
    assert id_conv.encode(10) == "TS000A"
    return


def test_decodes():
    id_conv = IdConverter(prefix="TS", length=4)
    assert id_conv.decode("TS0000") == 0
    assert id_conv.decode("TS000A") == 10

    id_conv = IdConverter(prefix="S", length=6)
    assert id_conv.decode("S000000") == 0
    assert id_conv.decode("S00000A") == 10
    return


def test_next():
    id_conv = IdConverter(prefix="TS", length=4)

    assert next(id_conv) == "TS0000"
    assert next(id_conv) == "TS0001"
    assert next(id_conv) == "TS0002"
    assert next(id_conv) == "TS0003"
    assert next(id_conv) == "TS0004"
    assert next(id_conv) == "TS0005"
    assert next(id_conv) == "TS0006"
    assert next(id_conv) == "TS0007"
    assert next(id_conv) == "TS0008"
    return
