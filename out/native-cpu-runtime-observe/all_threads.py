"""One fixed3s USER-only Linux x86_64 capture of a live owned process's threads.

Usage: all_threads.py PID START_TICKS SESSION PHASE_DIR
Reserve PHASE_DIR/all-threads-observation.json exclusively; never retry/overwrite.
The live owner.json must identify the same clean source, session child and ICD.
At most256 stable threads, one100Hz software CPU-clock event and9-page ring each.
This instruments the existing phase; it provides no uninstrumented benefit claim.
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
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
PERIOD, SECONDS, THREAD_CAP = 10_000_000, 3.0, 256


def need(condition, reason):
    if not condition:
        raise ValueError(reason)


def pin(path):
    before = path.stat()
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    after = path.stat()
    need((before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) ==
         (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns), 'file changed while hashing')
    return dict(bytes=after.st_size, sha256=digest)


def proc(pid, tid=None):
    path = Path('/proc') / str(pid)
    fields = (path / ('stat' if tid is None else f'task/{tid}/stat')).read_text().rsplit(')', 1)[1].split()
    return dict(pid=pid, tid=pid if tid is None else tid, state=fields[0], ppid=int(fields[1]),
                session=int(fields[3]), start_ticks=fields[19], cpu_ticks=int(fields[11]) + int(fields[12]))


def threads(pid):
    paths = sorted((Path('/proc') / str(pid) / 'task').iterdir(), key=lambda p: int(p.name))
    need(0 < len(paths) <= THREAD_CAP, 'empty thread set or thread cap exceeded')
    result = {}
    for path in paths:
        tid = int(path.name)
        fields = dict(line.split(':', 1) for line in (path / 'status').read_text().splitlines() if ':' in line)
        need(int(fields['Tgid']) == pid and int(fields['Pid']) == tid, 'task/TGID mismatch')
        result[tid] = proc(pid, tid)
    return result


def snapshot(value):
    normalized = {str(tid): row for tid, row in value.items()}
    raw = json.dumps(normalized, sort_keys=True, separators=(',', ':')).encode()
    return dict(threads=normalized, sha256=hashlib.sha256(raw).hexdigest())


def births(value):
    return {tid: row['start_ticks'] for tid, row in value.items()}


def source():
    environment = dict(os.environ, GIT_OPTIONAL_LOCKS='0')
    def git(*args):
        return subprocess.check_output(['git', '-C', str(ROOT), *args], env=environment, text=True, timeout=5).strip()
    return dict(revision=git('rev-parse', '--verify', 'HEAD'), status=git('status', '--porcelain', '--untracked-files=no'))


def stopped(number, _frame):
    raise InterruptedError(f'observer signal {number}')


def decode(event, pid):
    buffer, row = event['buffer'], event['record']
    head, tail, offset, size = struct.unpack_from('<QQQQ', buffer, 1024)
    row['ring'] = dict(head=head, tail=tail, data_offset=offset, data_size=size)
    need(offset >= mmap.PAGESIZE and offset + size <= len(buffer) and size > 0 and
         0 <= tail <= head and head - tail <= size, 'invalid or overrun perf ring')
    data = bytes(buffer[offset:offset + size]); doubled = data + data
    position = tail
    while position < head:
        index = position % size
        kind, misc, length = struct.unpack('<IHH', doubled[index:index + 8])
        need(8 <= length <= min(size, head - position), 'invalid perf record size')
        entry = doubled[index:index + length]; position += length
        if kind != 9:
            row['non_sample_record'] = dict(type=kind, misc=misc, size=length, prefix_hex=entry[:64].hex())
            raise ValueError(f'non-sample/lost/throttled record {kind}; no retry')
        need(length == 24, 'unexpected IP/TID record layout')
        ip, sample_pid, sample_tid = struct.unpack_from('<QII', entry, 8)
        need((misc & 7) == 2 and sample_pid == pid and sample_tid == event['tid'], 'non-USER or wrong PID/TID sample')
        row['samples'].append(dict(ip=ip, pid=sample_pid, tid=sample_tid))
    row['decoded'] = True


if len(sys.argv) != 5:
    raise SystemExit(__doc__)
phase = Path(sys.argv[4]).resolve(strict=True)
need(phase.is_dir() and phase.is_relative_to(ROOT / 'out'), 'phase must be inside checkout out/')
with (phase / 'all-threads-observation.json').open('x') as journal:
    record = dict(status='started', attempted=False, period_ns=PERIOD, requested_seconds=SECONDS,
                  thread_cap=THREAD_CAP, attr_bytes=64, ring_pages_per_thread=9,
                  exclude_kernel=True, exclude_hypervisor=True, inherit=False,
                  method='perf_event_open298 attr64 software CPU-clock USER IP/TID',
                  threads={}, errors=[], incomplete_reasons=[], cleanup_errors=[])
    json.dump(record, journal); journal.flush()
    events, postcheck, completed_window = [], None, False
    try:
        signal.signal(signal.SIGTERM, stopped); signal.signal(signal.SIGINT, stopped)
        need(platform.system() == 'Linux' and platform.machine() == 'x86_64' and ctypes.sizeof(ctypes.c_void_p) == 8,
             'only Linux x86_64 is supported')
        pid, start, session = int(sys.argv[1]), sys.argv[2], int(sys.argv[3])
        need(pid > 0 and session > 0 and start.isdecimal(), 'invalid process identity')
        owner_path = phase / 'owner.json'
        need(not owner_path.is_symlink() and owner_path.stat().st_uid == os.getuid(), 'unowned/symlink owner file')
        owner = json.loads(owner_path.read_text()); boot = Path('/proc/sys/kernel/random/boot_id').read_text().strip()
        expected_source = source(); need(not expected_source['status'], 'tracked source tree is dirty')
        need(owner['source'] == expected_source and owner['boot_id'] == boot and owner['status'] != 'terminal',
             'owner source/boot/live identity mismatch')
        leader, leader_start = owner['child_pid'], owner['child_start_ticks']
        value = owner['declared_output_environment']['VK_DRIVER_FILES']
        record.update(source_before=expected_source, observer=pin(Path(__file__)), owner_before=pin(owner_path))

        def identity():
            current, leading = proc(pid), proc(leader)
            fresh = json.loads(owner_path.read_text())
            need(fresh['status'] != 'terminal' and fresh['boot_id'] == boot and fresh['source'] == expected_source and
                 fresh['child_pid'] == leader and fresh['child_start_ticks'] == leader_start and
                 fresh['declared_output_environment']['VK_DRIVER_FILES'] == value, 'live owner identity changed')
            need(current['start_ticks'] == start and current['session'] == session and current['state'] not in ('Z', 'X'),
                 'target birth/session/liveness changed')
            need(leading['start_ticks'] == leader_start and leading['session'] == session and
                 leading['state'] not in ('Z', 'X') and os.stat(f'/proc/{pid}').st_uid == os.getuid(), 'owner child identity changed')
            ancestor = current
            for _ in range(64):
                if ancestor['pid'] == leader:
                    return current
                ancestor = proc(ancestor['ppid'])
                need(ancestor['session'] == session, 'target left the owned session')
            raise ValueError('target is not an owned-session descendant')

        record['target_before'] = identity()
        environment = dict(part.split(b'=', 1) for part in Path(f'/proc/{pid}/environ').read_bytes().split(b'\0') if b'=' in part)
        need(os.fsdecode(environment[b'VK_DRIVER_FILES']) == value and ':' not in value, 'target/owner ICD mismatch')
        icd = Path(value); icd = (icd if icd.is_absolute() else Path(os.readlink(f'/proc/{pid}/cwd')) / icd).resolve(strict=True)
        library = Path(json.loads(icd.read_text())['ICD']['library_path'])
        library = (library if library.is_absolute() else icd.parent / library).resolve(strict=True)
        executable = Path(os.readlink(f'/proc/{pid}/exe')).resolve(strict=True)
        need(os.path.samefile(f'/proc/{pid}/exe', executable), 'mapped executable differs from path')
        expected_pins = {str(p): pin(p) for p in (icd, library, executable)}
        record['inputs'] = expected_pins

        def maps():
            lines = Path(f'/proc/{pid}/maps').read_text().splitlines(); metadata = library.stat()
            matches = [line.split(None, 5) for line in lines if line.split(None, 5)[-1] == str(library)]
            need(any('x' in row[1] for row in matches), 'provider has no executable mapping')
            need(all(int(row[4]) == metadata.st_ino and tuple(int(x, 16) for x in row[3].split(':')) ==
                     (os.major(metadata.st_dev), os.minor(metadata.st_dev)) for row in matches), 'provider map inode mismatch')
            return lines, matches

        record['maps_before'], provider_maps = maps()
        first = threads(pid); time.sleep(0.2); stable = threads(pid)
        record['threads_before_first'] = snapshot(first); record['threads_before'] = snapshot(stable)
        need(births(first) == births(stable), 'thread set changed during stability snapshots')
        record['threads'] = {str(tid): dict(identity=row, samples=[]) for tid, row in stable.items()}

        def after():
            try:
                record['target_after'] = identity()
            except Exception as error:
                record['incomplete_reasons'].append(f'target/owner after: {type(error).__name__}: {error}')
            try:
                last = threads(pid); record['threads_after'] = snapshot(last)
                if births(stable) != births(last):
                    record['incomplete_reasons'].append('thread births/exits/reuse between boundary snapshots')
            except Exception as error:
                record['incomplete_reasons'].append(f'threads after: {type(error).__name__}: {error}')
            record['source_after'] = source(); need(record['source_after'] == expected_source, 'source changed')
            record['inputs_after'] = {name: pin(Path(name)) for name in expected_pins}
            need(record['inputs_after'] == expected_pins, 'provider/ICD/executable bytes changed')
            need(os.path.samefile(f'/proc/{pid}/exe', executable), 'target executable identity changed')
            record['maps_after'], later_maps = maps()
            need(later_maps == provider_maps, 'provider mappings changed')

        postcheck = after
        attr = (ctypes.c_ubyte * 64).from_buffer_copy(struct.pack('<IIQQQQQIIQ', 1, 64, 0, PERIOD, 3, 0, 97, 1, 0, 0))
        libc = ctypes.CDLL(None, use_errno=True); libc.syscall.restype = ctypes.c_long
        record['attempted'] = True
        for tid, row in stable.items():
            fd = libc.syscall(ctypes.c_long(298), ctypes.byref(attr), ctypes.c_long(tid),
                              ctypes.c_long(-1), ctypes.c_long(-1), ctypes.c_ulong(8))
            if fd < 0:
                raise OSError(ctypes.get_errno(), f'TID{tid}: {os.strerror(ctypes.get_errno())}')
            event = dict(tid=tid, fd=fd, buffer=None, enabled=False, record=record['threads'][str(tid)])
            events.append(event)
            need(proc(pid, tid)['start_ticks'] == row['start_ticks'], 'TID birth changed after open')
            event['buffer'] = mmap.mmap(fd, mmap.PAGESIZE * 9, flags=mmap.MAP_SHARED, prot=mmap.PROT_READ | mmap.PROT_WRITE)
            fcntl.ioctl(fd, 0x2403, 0)
        identity(); need(births(threads(pid)) == births(stable), 'thread set changed while preparing events')
        began = time.monotonic(); record['window_begin_monotonic'] = began
        for event in events:
            fcntl.ioctl(event['fd'], 0x2400, 0); event['enabled'] = True
            event['record']['enable_offset_seconds'] = time.monotonic() - began
        remaining = began + SECONDS - time.monotonic()
        need(remaining > 0, 'enabling events exhausted the single3s window')
        time.sleep(remaining); completed_window = True
    except BaseException as error:
        record['errors'].append(f'{type(error).__name__}: {error}')
    finally:
        signal.signal(signal.SIGTERM, signal.SIG_IGN); signal.signal(signal.SIGINT, signal.SIG_IGN)
        for event in events:
            if event['enabled']:
                try:
                    fcntl.ioctl(event['fd'], 0x2401, 0); event['enabled'] = False
                    event['record']['disable_offset_seconds'] = time.monotonic() - record['window_begin_monotonic']
                except BaseException as error:
                    record['cleanup_errors'].append(f'disable TID{event["tid"]}: {error}')
        if 'window_begin_monotonic' in record:
            record['actual_window_wall_seconds'] = time.monotonic() - record['window_begin_monotonic']
        for event in events:
            try:
                if event['buffer'] is not None and not event['enabled']:
                    decode(event, pid)
                    count = os.read(event['fd'], 8); need(len(count) == 8, 'invalid counter read')
                    event['record']['raw_CPU_clock_counter'] = struct.unpack('<Q', count)[0]
            except BaseException as error:
                record['errors'].append(f'decode TID{event["tid"]}: {type(error).__name__}: {error}')
            finally:
                for label, close in [('ring', event['buffer'].close if event['buffer'] is not None else None),
                                     ('fd', lambda fd=event['fd']: os.close(fd))]:
                    if close is not None:
                        try:
                            close()
                        except BaseException as error:
                            record['cleanup_errors'].append(f'close {label} TID{event["tid"]}: {error}')
        if postcheck is not None:
            try:
                postcheck()
            except BaseException as error:
                record['errors'].append(f'postcheck: {type(error).__name__}: {error}')
        record['sample_count'] = sum(len(row['samples']) for row in record['threads'].values())
        record['event_count_opened'] = len(events)
        record['completed_requested_window'] = completed_window
        record['status'] = ('error' if record['errors'] or record['cleanup_errors'] or not completed_window else
                            'incomplete' if record['incomplete_reasons'] else 'complete')
        record['limits'] = ('USER IP samples only, approximate100Hz per-thread CPU-clock exposure; no callchains, '
                            'kernel attribution or whole-frame causality. Event enable/disable is sequential; per-TID '
                            'offsets preserve skew and scheduler overshoot. Boundary snapshots cannot detect threads '
                            'born and dead entirely between them; no inheritance. Empty per-TID rings may mean asleep, '
                            'off-CPU, kernel execution or insufficient sampling exposure. Partial samples survive errors '
                            'or changed thread sets; only complete denotes unchanged observed boundary coverage. No retry.')
        journal.seek(0); journal.truncate(); json.dump(record, journal, indent=2); journal.write('\n'); journal.flush()
    raise SystemExit(0 if record['status'] == 'complete' else 1)
