"""Read an exact predecessor script; require every declared substitution."""
from pathlib import Path
import hashlib

R = Path.cwd()


def source(name, digest):
    path = R / name
    data = path.read_bytes()
    assert hashlib.sha256(data).hexdigest() == digest, path
    return data.decode()


def replace(text, before, after):
    assert text.count(before) == 1, before
    return text.replace(before, after, 1)
