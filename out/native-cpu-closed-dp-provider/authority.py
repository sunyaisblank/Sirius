"""Read an exact predecessor script; require every declared substitution."""
from pathlib import Path
import hashlib

R = Path.cwd()


def source(name, digest):
    data = (R / name).read_bytes()
    assert hashlib.sha256(data).hexdigest() == digest, name
    return data.decode()


def replace(text, before, after):
    assert text.count(before) == 1, before
    return text.replace(before, after, 1)
