"""Decode existing USER samples only after terminal owner; usage: SCRIPT PHASE_DIR."""
from pathlib import Path
from collections import Counter, defaultdict
import hashlib
import json
import os
import struct
import subprocess
import sys
ROOT = Path(__file__).resolve().parents[3]
PROOF = ROOT / 'out/native-cpu-accumulate-inline/provider/linked-native-proof.json'
def need(condition, message):
    if not condition:
        raise ValueError(message)
def pin(path):
    a = path.stat()
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    b = path.stat()
    need((a.st_dev, a.st_ino, a.st_size, a.st_mtime_ns) ==
         (b.st_dev, b.st_ino, b.st_size, b.st_mtime_ns), 'file changed while hashing')
    return {'bytes': b.st_size, 'sha256': digest}

def document(path):
    return json.loads(path.read_text())

def mapping(line):
    words = line.split(None, 5)
    lo, hi = (int(v, 16) for v in words[0].split('-'))
    return lo, hi, words

def load_address(path, offset):
    with path.open('rb') as stream:
        header = stream.read(64)
        need(header[:7] == b'\x7fELF\x02\x01\x01', 'not little-endian ELF64')
        need(struct.unpack_from('<H', header, 18)[0] == 62, 'not x86_64 ELF')
        table = struct.unpack_from('<Q', header, 32)[0]
        size, count = struct.unpack_from('<HH', header, 54)
        need(size == 56 and 0 < count < 1024, 'unsupported ELF program headers')
        stream.seek(table)
        segments = [struct.unpack('<IIQQQQQQ', stream.read(size)) for _ in range(count)]
    matches = [vaddr + offset - start for kind, flags, start, vaddr, _, length, _, _ in segments
               if kind == 1 and flags & 1 and start <= offset < start + length]
    need(len(matches) == 1, 'file offset not in one executable PT_LOAD')
    return matches[0]

need(len(sys.argv) == 2, __doc__)
phase = Path(sys.argv[1]).resolve(strict=True)
need(phase.is_dir() and phase.is_relative_to(ROOT / 'out'), 'phase outside checkout out/')
output = phase / 'all-threads-symbols.json'
with output.open('x') as journal:
    result = {'status': 'started', 'errors': [], 'addresses': {}, 'threads': {}, 'aggregate': {}}
    json.dump(result, journal); journal.flush()
    try:
        raw_path, owner_path = phase / 'all-threads-observation.json', phase / 'owner.json'
        before_path, after_path = phase / 'inputs-before.json', phase / 'inputs-after.json'
        result['inputs'] = {str(p): pin(p) for p in (raw_path, owner_path, before_path, after_path, PROOF, Path(__file__))}
        raw, owner = document(raw_path), document(owner_path)
        need(owner['status'] == 'terminal', 'owner must be terminal; no live postprocessing')
        need(owner['source_unchanged'] and owner['whole_inputs_unchanged'], 'owner inputs/source changed')
        need(not owner['cleanup_errors'] and not owner['remaining_owned_processes'], 'owner cleanup incomplete')
        if raw.get('sample_count', 0):
            target = raw['target_before']
            need(target['session'] == owner['child_pid'] and any(row['pid'] == target['pid'] and
                 row['saved_start_ticks'] == target['start_ticks'] and row['absent']
                 for row in owner['observed_births']), 'captured target not a reaped owned birth')
        before = document(before_path)
        need(before == document(after_path), 'frozen owner ledgers differ')
        need(raw.get('source_before') == owner['source'], 'capture/owner source mismatch')
        need(raw.get('source_after', owner['source']) == owner['source'], 'capture source changed')
        ledger = {str(Path(v.get('resolved_path', ROOT / k)).resolve()):
                  {field: v[field] for field in ('bytes', 'sha256')} for k, v in before.items()}
        pinned = {}
        def sealed(path):
            path = path.resolve(strict=True)
            if str(path) not in pinned:
                need(str(path) in ledger, 'mapped file/proof absent from frozen owner ledger')
                pinned[str(path)] = pin(path)
                need(pinned[str(path)] == ledger[str(path)], 'current file differs from frozen owner ledger')
            return path
        for name, expected in raw.get('inputs', {}).items():
            path = sealed(Path(name))
            need(pinned[str(path)] == expected, 'capture/current file identity mismatch')
        proof_path = sealed(PROOF)
        proof = document(proof_path)
        library = sealed(ROOT / proof['library'])
        need(proof['passed'] and pinned[str(library)] == proof['library_identity'], 'provider proof identity mismatch')
        ranges = []
        for obj in proof['objects']:
            need(obj['all_retained_function_opcodes_match'], 'unverified native object code')
            ranges += [dict(object=obj['object']['path'], **{k: f[k] for k in ('linked_address', 'bytes', 'symbol')})
                       for f in obj['functions']]
        maps = [mapping(line) for line in raw.get('maps_before', [])]
        later = set(raw.get('maps_after', []))
        unique_ips = sorted({sample['ip'] for row in raw['threads'].values() for sample in row['samples']})
        batches = defaultdict(list)
        for ip in unique_ips:
            entry = {'ip': ip, 'resolved': False}
            result['addresses'][hex(ip)] = entry
            try:
                matches = [(lo, hi, w) for lo, hi, w in maps if lo <= ip < hi and 'x' in w[1]]
                need(len(matches) == 1, 'no unique executable map')
                lo, _, words = matches[0]
                need(len(words) == 6 and words[5].startswith('/') and not words[5].endswith(' (deleted)'), 'anonymous/deleted map')
                need(any(mapping(line)[2] == words for line in later), 'mapping not confirmed at both boundaries')
                path = sealed(Path(words[5])); stat = path.stat()
                need(stat.st_ino == int(words[4]) and (os.major(stat.st_dev), os.minor(stat.st_dev)) ==
                     tuple(int(v, 16) for v in words[3].split(':')), 'map/current device-inode mismatch')
                offset = int(words[2], 16) + ip - lo
                address = load_address(path, offset)
                native = [r for r in ranges if r['linked_address'] <= address < r['linked_address'] + r['bytes']] if path == library else []
                entry.update(file=str(path), file_offset=offset, elf_vaddr=address, native_ranges=native)
                batches[path].append(entry)
            except (OSError, ValueError, struct.error) as error:
                entry['unresolved_reason'] = str(error)
        tool = Path('/usr/bin/addr2line').resolve(strict=True)
        result['addr2line'] = {'path': str(tool), **pin(tool)}
        for path, entries in batches.items():
            need(len(entries) <= 100000, 'per-ELF sample-address cap exceeded')
            decoded = subprocess.run([str(tool), '-f', '-C', '-e', str(path)],
                input=''.join(hex(row['elf_vaddr']) + '\n' for row in entries), text=True,
                capture_output=True, check=True, timeout=30).stdout.splitlines()
            need(len(decoded) == 2 * len(entries), 'addr2line response count mismatch')
            for i, row in enumerate(entries):
                row.update(function=decoded[2*i], location=decoded[2*i+1], resolved=True, symbol_resolved=decoded[2*i] != '??')
        total = Counter()
        for tid, thread in raw['threads'].items():
            counts = Counter()
            for sample in thread['samples']:
                row = result['addresses'][hex(sample['ip'])]
                key = json.dumps({k: row.get(k) for k in ('file', 'function', 'location', 'unresolved_reason')} |
                    {'native_objects': sorted({r['object'] for r in row.get('native_ranges', [])}),
                     'native_symbols': sorted({r['symbol'] for r in row.get('native_ranges', [])})}, sort_keys=True)
                counts[key] += 1
            result['threads'][tid] = {'identity': thread['identity'], 'samples': sum(counts.values()),
                'functions': [dict(json.loads(k), samples=v) for k, v in counts.most_common()]}
            total.update(counts)
        need(all(pin(Path(name)) == expected for name, expected in pinned.items()), 'mapped/proof bytes changed while decoding')
        need(all(pin(Path(name)) == expected for name, expected in result['inputs'].items()), 'input record changed while decoding')
        result.update(status='decoded' if raw['status'] == 'complete' else 'decoded_partial',
            capture_status=raw['status'], capture_errors=raw['errors'], capture_cleanup_errors=raw['cleanup_errors'],
            capture_incomplete_reasons=raw['incomplete_reasons'], capture_limits=raw.get('limits'),
            aggregate={'samples': sum(total.values()), 'functions': [dict(json.loads(k), samples=v) for k, v in total.most_common()]},
            native_function_ranges=len(ranges), verified_current_files=pinned,
            limits='Recorded USER IPs only; partial capture stays partial. Boundary maps do not prove no transient remapping. Anonymous/non-ELF/unsealed maps stay unresolved; no JIT inference, callchain, whole-frame attribution or causality.')
        need(sum(total.values()) == raw['sample_count'], 'raw sample count mismatch')
    except Exception as error:
        result['status'] = 'error'; result['errors'].append(f'{type(error).__name__}: {error}')
    journal.seek(0); journal.truncate(); json.dump(result, journal, indent=2); journal.write('\n')
raise SystemExit(0 if result['status'] in ('decoded', 'decoded_partial') else 1)
