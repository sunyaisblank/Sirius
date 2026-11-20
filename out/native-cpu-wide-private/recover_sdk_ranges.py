"""Restore the pinned LLVM SDK members from verified release-archive ranges."""

from pathlib import Path, PurePosixPath
import hashlib
import io
import lzma
import shutil
import struct
import sys
import tarfile
import tempfile
import urllib.request
import zlib

from resume_raw import R, check_space, pin, read_json
from run_strict_scope import HISTORICAL_SDK


CATALOG = R / 'attestations/qualification/432d65a9/native-cpu-simd16-value/sdk-ranges.json'
SDK_ROOT = 'bin/linux-gcc/native-cpu-strict-ir-linking/toolchain'
URL = 'https://github.com/llvm/llvm-project/releases/download/llvmorg-23.1.2/LLVM-23.1.2-Linux-X64.tar.xz'
SUPPORT = Path('lib/libLLVMSupport.a')
HEADERS = (Path('include/llvm'), Path('include/llvm-c'), Path('lib/clang/23/include'))


def integer(value):
    encoded = bytearray()
    while value >= 128:
        encoded.append((value & 127) | 128)
        value >>= 7
    encoded.append(value)
    return encoded


def trailer(header, records):
    # XZ blocks are independent. Build an index and footer describing only
    # these original blocks. Their bytes and integrity checks stay exact.
    index = b'\0' + integer(len(records)) + b''.join(integer(packed) + integer(raw) for packed, raw in records)
    index += b'\0' * (-len(index) % 4)
    index += struct.pack('<I', zlib.crc32(index))
    footer = struct.pack('<I', len(index) // 4 - 1) + header[6:8]
    return index + struct.pack('<I', zlib.crc32(footer)) + footer + b'YZ'


class PackedRange:
    def __init__(self, file, header, span):
        self.file = file
        self.parts = [io.BytesIO(header), None, io.BytesIO(trailer(header, span['index_records']))]
        self.left = span['compressed_bytes']
        self.digest = hashlib.sha256()

    def read(self, size):
        while size and self.parts:
            part = self.parts[0]
            if part is None:
                data = self.file.read(min(size, self.left))
                if self.left and not data:
                    raise RuntimeError('Truncated SDK range')
                self.left -= len(data)
                self.digest.update(data)
            else:
                data = part.read(size)
            if data:
                return data
            self.parts.pop(0)
        return b''


class TarRegion:
    def __init__(self, payload, length):
        self.payload = payload
        self.left = length
        self.padding = 1024

    def read(self, size):
        if self.left:
            data = self.payload.read(min(size, self.left))
            if not data:
                raise RuntimeError('Truncated SDK tar region')
            self.left -= len(data)
            return data
        size = min(size, self.padding)
        self.padding -= size
        return b'\0' * size


def selected_members(file, header, span, archive_root, selected, consume):
    packed = PackedRange(file, header, span)
    with lzma.LZMAFile(packed) as payload:
        skip = span['tar_skip']
        while skip:
            data = payload.read(min(skip, 1024 * 1024))
            if not data:
                raise RuntimeError('Truncated SDK tar prefix')
            skip -= len(data)
        with tarfile.open(mode='r|', fileobj=TarRegion(payload, span['tar_bytes'])) as bundle:
            for member in bundle:
                path = PurePosixPath(member.name)
                if not path.parts or path.is_absolute() or '..' in path.parts or path.parts[0] != archive_root:
                    raise RuntimeError(f'Unexpected SDK tar path {path}')
                relative = Path(*path.parts[1:])
                if member.isfile() and (relative in selected or any(relative.is_relative_to(p) for p in HEADERS)):
                    with bundle.extractfile(member) as stream:
                        consume(relative, member, stream)
        while payload.read(1024 * 1024):
            pass
    if packed.left or packed.digest.hexdigest() != span['sha256']:
        raise RuntimeError('SDK range SHA-256 differs')


def download_span(span, archive_bytes, path):
    start = span['compressed_offset']
    size = span['compressed_bytes']
    end = start + size - 1
    request = urllib.request.Request(URL, headers={
        'Range': f'bytes={start}-{end}', 'User-Agent': 'Sirius-SDK-recovery'})
    with urllib.request.urlopen(request, timeout=90) as response, path.open('xb') as output:
        if response.status != 206 or response.headers.get('Content-Range') != f'bytes {start}-{end}/{archive_bytes}':
            raise RuntimeError('Unexpected LLVM SDK range response')
        remaining = size
        while data := response.read(min(remaining + 1, 1024 * 1024)):
            remaining -= len(data)
            if remaining < 0:
                raise RuntimeError('SDK range exceeded its pinned extent')
            output.write(data)
    if remaining or pin(path) != {'bytes': size, 'sha256': span['sha256']}:
        raise RuntimeError('SDK download differs from its pinned range')


def descriptor():
    sdk = read_json(HISTORICAL_SDK)['setup']['official_sdk']
    catalog = read_json(CATALOG)
    if catalog['archive'] != sdk['archive']:
        raise RuntimeError('SDK slice catalog has the wrong archive identity')
    entries = {Path(item['resolved']).relative_to(SDK_ROOT): item for item in sdk['tools'].values()}
    return sdk, catalog, entries


def recover(scratch):
    check_space()
    if scratch == R / 'out' or not scratch.is_relative_to(R / 'out'):
        raise RuntimeError('Use a per-task directory under out/')
    sdk, catalog, entries = descriptor()
    header = bytes.fromhex(catalog['stream_header_hex'])
    destination = scratch / 'sdk'
    destination.mkdir()
    selected = set(entries) | {SUPPORT}
    seen = set()

    def extract(relative, member, stream):
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('xb') as output:
            shutil.copyfileobj(stream, output, length=1024 * 1024)
        target.chmod(member.mode & 0o777)
        seen.add(relative)

    try:
        with tempfile.TemporaryDirectory(prefix='sdk-ranges-', dir=scratch) as directory:
            compressed = Path(directory) / 'span.xz-blocks'
            for span in catalog['ranges']:
                check_space()
                download_span(span, catalog['archive']['bytes'], compressed)
                with compressed.open('rb') as file:
                    selected_members(file, header, span, sdk['extracted_root'], selected, extract)
                compressed.unlink()
        if not selected <= seen or len(seen) != catalog['member_count']:
            raise RuntimeError('Selected SDK files are incomplete')
        for relative, item in entries.items():
            if pin(destination / relative) != {key: item[key] for key in ('bytes', 'sha256')}:
                raise RuntimeError(f'Recovered compiler fingerprint differs: {relative}')
        return destination
    except BaseException:
        shutil.rmtree(destination)
        raise


if __name__ == '__main__':
    location, = sys.argv[1:]
    print(recover(Path(location).resolve()))
