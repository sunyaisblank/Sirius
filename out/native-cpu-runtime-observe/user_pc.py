"""One USER-only Linux x86_64 sample of an already-owned live process.

Usage: user_pc.py PID START_TICKS SESSION PHASE_DIR
PHASE_DIR/owner.json must identify the live owned session leader. The fixed
user-pc-observation.json is exclusively reserved; never overwrite or retry it.
This instruments the phase: its timing is not uninstrumented benefit evidence.
"""
from pathlib import Path
import ctypes
import fcntl
import hashlib
import json
import mmap
import os
import platform
import signal
import struct
import sys
import time


def proc(pid, tid=None):
    path = Path('/proc') / str(pid)
    fields = (path / ('stat' if tid is None else f'task/{tid}/stat')).read_text().rsplit(')', 1)[1].split()
    return dict(pid=pid, tid=pid if tid is None else tid, state=fields[0], ppid=int(fields[1]),
                session=int(fields[3]), start_ticks=fields[19], cpu_ticks=int(fields[11]) + int(fields[12]))


def threads(pid):
    result = {}
    for path in (Path('/proc') / str(pid) / 'task').iterdir():
        try:
            result[int(path.name)] = proc(pid, int(path.name))
        except (FileNotFoundError, ProcessLookupError):
            pass
    return result


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def stopped(number, _frame):
    raise InterruptedError(f'observer signal {number}')


if len(sys.argv) != 5:
    raise SystemExit(__doc__)
phase = Path(sys.argv[4]).resolve(strict=True)
assert phase.is_dir() and phase.is_relative_to(Path(__file__).resolve().parents[2] / 'out')
with (phase / 'user-pc-observation.json').open('x') as journal:
    record = dict(status='started', attempted=False, period_ns=10_000_000, seconds=3.0,
                  method='perf_event_open298 attr64 software CPU-clock USER IP/TID only', samples=[])
    json.dump(record, journal); journal.flush()
    fd, buffer, enabled = -1, None, False
    try:
        signal.signal(signal.SIGTERM, stopped); signal.signal(signal.SIGINT, stopped)
        assert platform.system() == 'Linux' and platform.machine() == 'x86_64' and ctypes.sizeof(ctypes.c_void_p) == 8
        pid, start, session = int(sys.argv[1]), sys.argv[2], int(sys.argv[3])
        owner_path = phase / 'owner.json'
        owner = json.loads(owner_path.read_text())
        boot = Path('/proc/sys/kernel/random/boot_id').read_text().strip()
        assert owner['status'] != 'terminal' and owner['boot_id'] == boot
        leader, leader_start = owner['child_pid'], owner['child_start_ticks']

        def check_identity():
            current = proc(pid); leading = proc(leader)
            assert current['start_ticks'] == start and current['session'] == session
            assert leading['start_ticks'] == leader_start and leading['session'] == session
            assert os.stat(f'/proc/{pid}').st_uid == os.getuid()
            ancestor = current
            for _ in range(64):
                if ancestor['pid'] == leader:
                    return current
                ancestor = proc(ancestor['ppid'])
                assert ancestor['session'] == session
            raise ValueError('target is not an owned-session descendant')

        record['target_before'] = check_identity()
        value = owner['declared_output_environment']['VK_DRIVER_FILES']
        environment = dict(part.split(b'=', 1) for part in Path(f'/proc/{pid}/environ').read_bytes().split(b'\0') if b'=' in part)
        assert os.fsdecode(environment[b'VK_DRIVER_FILES']) == value and ':' not in value
        icd = Path(value)
        if not icd.is_absolute():
            icd = Path(os.readlink(f'/proc/{pid}/cwd')) / icd
        icd = icd.resolve(strict=True)
        library = Path(json.loads(icd.read_text())['ICD']['library_path'])
        library = (library if library.is_absolute() else icd.parent / library).resolve(strict=True)
        before_hash = pin(library); record['provider'] = dict(path=str(library), **before_hash)
        record['ICD'] = dict(path=str(icd), **pin(icd))

        def maps():
            lines = Path(f'/proc/{pid}/maps').read_text().splitlines(); metadata = library.stat()
            matches = [line.split(None, 5) for line in lines if line.split(None, 5)[-1] == str(library)]
            assert any('x' in line[1] for line in matches)
            assert all(int(line[4]) == metadata.st_ino and tuple(int(x, 16) for x in line[3].split(':')) ==
                       (os.major(metadata.st_dev), os.minor(metadata.st_dev)) for line in matches)
            return lines

        record['maps_before'] = maps()
        first = threads(pid); time.sleep(0.2); second = threads(pid)
        stable = {tid: second[tid]['cpu_ticks'] - item['cpu_ticks'] for tid, item in first.items()
                  if tid in second and second[tid]['start_ticks'] == item['start_ticks']}
        tid = max(stable, key=stable.get); assert stable[tid] > 0
        thread_start = second[tid]['start_ticks']; record['thread_cpu_deltas'] = stable
        record['selected_thread'] = second[tid]; check_identity()
        attr = ctypes.create_string_buffer(struct.pack('<IIQQQQQIIQ', 1, 64, 0, 10_000_000, 3, 0, 97, 1, 0, 0))
        libc = ctypes.CDLL(None, use_errno=True); libc.syscall.restype = ctypes.c_long
        record['attempted'] = True
        fd = libc.syscall(298, ctypes.byref(attr), tid, -1, -1, 8)
        if fd < 0:
            raise OSError(ctypes.get_errno(), os.strerror(ctypes.get_errno()))
        assert proc(pid, tid)['start_ticks'] == thread_start; check_identity()
        buffer = mmap.mmap(fd, mmap.PAGESIZE * 9, flags=mmap.MAP_SHARED, prot=mmap.PROT_READ | mmap.PROT_WRITE)
        fcntl.ioctl(fd, 0x2403, 0); fcntl.ioctl(fd, 0x2400, 0); enabled = True
        began = time.monotonic(); time.sleep(3.0)
        fcntl.ioctl(fd, 0x2401, 0); enabled = False; record['sample_wall_seconds'] = time.monotonic() - began
        head, tail, offset, size = struct.unpack_from('<QQQQ', buffer, 1024)
        assert offset >= mmap.PAGESIZE and offset + size <= len(buffer) and size and head - tail <= size
        data = bytes(buffer[offset:offset + size]); position = tail
        while position < head:
            index = position % size; header = (data + data)[index:index + 8]
            kind, misc, length = struct.unpack('<IHH', header); assert 8 <= length <= head - position
            entry = (data + data)[index:index + length]; position += length
            if kind == 9:
                ip, sample_pid, sample_tid = struct.unpack_from('<QII', entry, 8)
                assert misc & 7 == 2 and sample_pid == pid and sample_tid == tid
                record['samples'].append(dict(ip=ip, pid=sample_pid, tid=sample_tid))
            else:
                raise ValueError(f'non-sample perf record {kind}; no retry')
        record['target_after'] = check_identity(); assert proc(pid, tid)['start_ticks'] == thread_start
        record['maps_after'] = maps(); record['provider_after'] = pin(library)
        assert record['provider_after'] == before_hash
        record['status'] = 'complete'
    except BaseException as error:
        record.update(status='error', error=f'{type(error).__name__}: {error}')
    finally:
        if enabled:
            fcntl.ioctl(fd, 0x2401, 0)
        if buffer is not None:
            buffer.close()
        if fd >= 0:
            os.close(fd)
        journal.seek(0); journal.truncate(); json.dump(record, journal, indent=2); journal.write('\n')
    raise SystemExit(0 if record['status'] == 'complete' else 1)
