"""Audit the exact RGBA32F scanline ZIP files produced by the live disk tests."""

from hashlib import sha256
import json
from math import fsum, isfinite
from pathlib import Path
import struct
import sys
import zlib


def field_name(data, start):
    finish = data.index(0, start)
    return data[start:finish].decode('ascii'), finish + 1


def uncompress(data, length):
    if len(data) == length:
        return data
    if len(data) > length:
        raise ValueError('Incorrect packed byte count')
    shuffled = bytearray(zlib.decompress(data))
    if len(shuffled) != length:
        raise ValueError('Incorrect inflated byte count')
    for index in range(1, length):
        shuffled[index] = (shuffled[index - 1] + shuffled[index] - 128) & 255
    even = (length + 1) // 2
    interleaved = bytearray(length)
    interleaved[::2] = shuffled[:even]
    interleaved[1::2] = shuffled[even:]
    return interleaved


def read_frame(path):
    payload = path.read_bytes()
    if struct.unpack_from('<II', payload) != (20000630, 2):
        raise ValueError(f'Expected single-part OpenEXR scanline version 2 in {path}')
    attributes = {}
    start = 8
    while payload[start]:
        name, start = field_name(payload, start)
        kind, start = field_name(payload, start)
        size, = struct.unpack_from('<I', payload, start)
        start += 4
        raw = payload[start:start + size]
        if name in attributes or len(raw) != size:
            raise ValueError(f'Invalid header field {name} in {path}')
        attributes[name] = (kind, raw)
        start += size
    start += 1

    if attributes['compression'] != ('compression', b'\x03'):
        raise ValueError(f'Expected lossless ZIP compression in {path}')
    kind, window = attributes['dataWindow']
    if kind != 'box2i' or struct.unpack('<iiii', window) != (0, 0, 39, 23):
        raise ValueError(f'Expected the original 40 by 24 pixel region in {path}')
    kind, raw_channels = attributes['channels']
    ordered_channels = []
    offset = 0
    while raw_channels[offset]:
        name, offset = field_name(raw_channels, offset)
        channels = struct.unpack_from('<i4xii', raw_channels, offset)
        if channels != (2, 1, 1):
            raise ValueError(f'Expected unsampled FLOAT channel {name} in {path}')
        ordered_channels.append(name)
        offset += 16
    if kind != 'chlist' or sorted(ordered_channels) != list('ABGR'):
        raise ValueError(f'Expected exactly the four RGBA channels in {path}')

    by_row = {}
    for offset in struct.unpack_from('<QQ', payload, start):
        y, size = struct.unpack_from('<ii', payload, offset)
        if y not in (0, 16) or y in by_row or size <= 0:
            raise ValueError(f'Malformed ZIP chunk in {path}')
        first = offset + 8
        rows = min(16, 24 - y)
        packed = payload[first:first + size]
        if len(packed) != size:
            raise ValueError(f'Truncated ZIP chunk in {path}')
        data = uncompress(packed, rows * 4 * 40 * 4)
        floats = struct.unpack('<' + 'f' * (rows * 4 * 40), data)
        for index in range(rows):
            by_row[y + index] = floats[index * 160:(index + 1) * 160]
    values = tuple(value for y in range(24) for value in by_row[y])
    if not all(map(isfinite, values)):
        raise ValueError(f'Nonfinite output in {path}')
    return ordered_channels, values, {'path': str(path), 'bytes': len(payload),
                                       'sha256': sha256(payload).hexdigest()}


def compare(left, right):
    channels, a, left_pin = read_frame(left)
    actual_channels, b, right_pin = read_frame(right)
    if channels != actual_channels:
        raise ValueError('Channel order differs')
    deltas = [abs(x - y) for x,y in zip(a,b,strict=True)]
    return {'format': 'OpenEXR ZIP scanline FLOAT RGBA', 'width': 40, 'height': 24,
            'compared_samples': len(deltas), 'changed_samples': sum(x != 0 for x in deltas),
            'absolute_difference': fsum(deltas), 'maximum_absolute_difference': max(deltas),
            'inputs': [left_pin,right_pin]}


if __name__ == '__main__':
    left, right = (Path(arg) for arg in sys.argv[1:])
    print(json.dumps(compare(left, right), sort_keys=True))
