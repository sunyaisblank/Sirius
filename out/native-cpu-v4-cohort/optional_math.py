"""Read bounded optional-target evidence without changing original C8 verdicts."""
from pathlib import Path
import hashlib
import json
import subprocess


def check(root, action, observation):
    def pin(path):
        with path.open('rb') as stream:
            return {'bytes': path.stat().st_size,
                    'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}
    path = root / 'out/native-cpu-v4-provider/provider/linked-native-proof.json'
    proof = json.loads(path.read_text())
    library = Path(action['runtime_provider'])
    assert proof['passed'] and root / proof['library'] == library
    assert pin(library) == proof['library_identity']
    symbols = {}
    for line in subprocess.check_output(['/usr/bin/nm', '-S', '-a', str(library)], text=True).splitlines():
        words = line.split()
        if len(words) == 4 and (words[3].startswith('sirius_native_v4_') or words[3] == '__cpu_features2'):
            assert words[3] not in symbols
            symbols[words[3]] = int(words[0], 16)
    selected = {row['function'][len('sirius_simd_'):-len('_packet')]
                for row in observation['simd_calls']}
    intervals = observation['optional_intervals']
    assert set(intervals) == selected
    assert all(row['closed'] and row['raw_entry_sequence'] < row['raw_return_sequence']
               for row in intervals.values())
    rows = observation['optional_math_calls']
    observed = {row['stage'] for row in rows}
    assert len(rows) == len(observed) and observed <= selected
    raw = dict(zip(('transport', 'endpoint', 'dopri_phase'), proof['objects'][4:7]))
    exports = {'sirius_native_v4_' + kind + '_' + operation + '8'
               for kind in ('dense', 'tagged')
               for operation in ('add', 'multiply', 'subtract', 'divide', 'sqrt')}
    for row in rows:
        assert row['function'] in exports and row['stage'] == row['expected_stage']
        assert row['returned'] and row['completed']
        interval = intervals[row['stage']]
        assert interval['raw_entry_sequence'] < row['entry_sequence'] < row['return_sequence'] < interval['raw_return_sequence']
        assert row['thread'] == interval['raw_thread']
        assert row['thread'] == row['return_thread']
        assert row['mxcsr'] & 0xe040 == row['return_mxcsr'] & 0xe040 == 0
        assert row['library'] == str(library)
        bases = {int(line.split()[0].split('-')[0], 16) - int(line.split()[2], 16)
                 for line in row['provider_maps']}
        base = row['pc'] - symbols[row['function']]
        assert base == row['provider_base'] and base in bases
        assert row['cpu_features2_address'] == base + symbols['__cpu_features2']
        data = bytes.fromhex(row['cpu_features2_first12'])
        assert len(data) == 12 and data[8] & 4
        assert row['scratch'] and 0 < row['mask'] <= 255
        assert row['raw_caller'] in row['caller_frames'] and len(row['caller_frames']) <= 16
        pc = row['raw_caller']['pc'] - base
        assert any(item['linked_address'] <= pc < item['linked_address'] + item['bytes'] or
                   item['linked_address'] <= pc - 1 < item['linked_address'] + item['bytes']
                   for item in raw[row['stage']]['functions'])
    return {'selected_simd_stages': sorted(selected),
            'observed_optional_math_stages': sorted(observed),
            'unobserved_optional_math_stages': sorted(selected - observed),
            'all_selected_optional_math_observed': observed == selected,
            'provider_linked_proof': pin(path),
            'scope': 'First optional arithmetic entry/return within each first observed selected raw '
                     'SIMD interval only. Exact loaded provider/CPU-cache/caller-PC/thread/FP-control '
                     'joins; no inferior calls or data/feature/FP edits. An unobserved interval is '
                     'reported without resampling. Original scalar-selected Transport may have no '
                     'optional arithmetic. This does not observe older-hardware fallback, measure '
                     'benefit, complete a frame or establish scientific/release qualification.'}
