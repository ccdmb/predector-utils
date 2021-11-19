#!/usr/bin/env bash

import hashlib
import base64


def checksum(s: bytes) -> bytes:
    sha256 = hashlib.sha256(s).digest()
    return base64.b64encode(sha256).rstrip(b"=")
