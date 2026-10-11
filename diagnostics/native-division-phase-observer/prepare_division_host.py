"""Prepare diagnostic native host sources without changing the production checkout."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, re, struct, subprocess
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
BASELINE = '2053ef79c90752d37b0075b8137cfb9b571b6820'
ANCESTOR = '24a8e145d1da5dca75a4fa26ec8e639020eecf37'
NAMES = ('src/sirius/backend/vulkan/vulkan_device.cpp', 'tests/backend/retained_compute_test.cpp',
         'tests/backend/shader_clock_observer.h')
def git(*a):
    return subprocess.check_output(['git', *a], timeout=30)
def read(p):
    return json.loads(p.read_text())
def seal(p):
    return {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
def one(s, a, b):
    assert s.count(a) == 1, (a, s.count(a))
    return s.replace(a, b)
assert ROOT.resolve() == Path('/home/astra/.project/Sirius')
assert WORK == ROOT / 'out/native-division-phase-observer'
assert git('rev-parse', 'HEAD').decode().strip() == BASELINE == git('rev-parse', '@{u}').decode().strip()
assert not git('status', '--porcelain')
placement = read(WORK / 'actual-placement-review.json')
assert placement['accepted'] and placement['baseline_revision'] == BASELINE
assert placement['original_finite_runtime_controls_admitted']
required = (WORK / 'placement/plan.json', WORK / 'placement/placement-query-bridge.json',
            WORK / 'placement/clock-placement.json', WORK / 'placement-analysis.json',
            WORK / 'placement/native-retirement.json', WORK / 'placement/retirement-run.json',
            WORK / 'compiler.json', WORK / 'compiled-review.json',
            WORK / 'phase/retained_transport.spv', WORK / 'phase/retained_transport.slang')
assert all(str(p) in placement['reviewed_seals'] for p in required)
for name, expected in placement['reviewed_seals'].items():
    assert seal(Path(name)) == expected
compiled = read(WORK / 'compiler.json')
assert compiled['accepted'] and seal(WORK / 'phase/retained_transport.spv') == compiled['product']['spv']
assert read(WORK / 'placement/plan.json')['source_revision'] == BASELINE
bridge = read(WORK / 'placement/placement-query-bridge.json')
assert bridge['accepted'] and bridge['clock_placement']['dispatches'] == 0
assert read(WORK / 'placement-analysis.json')['stages']['Transport']['module'] == compiled['product']['spv']
module = (WORK / 'phase/retained_transport.spv').read_bytes()
assert len(module) % 4 == 0
header = git('show', ANCESTOR + ':' + NAMES[2]).decode()
source = header
words = struct.unpack('<' + 'I' * (len(module) // 4), module)
array = 'inline constexpr std::uint32_t kTransportClock[] = {\n'
for i in range(0, len(words), 7):
    array += '    ' + ', '.join(f'0x{v:08x}u' for v in words[i:i + 7]) + ',\n'
array += '};'
header, n = re.subn(r'inline constexpr std::uint32_t kTransportClock\[\] = \{.*?\n\};', array, header, count=1, flags=re.S)
assert n == 1
header = one(header, 'e3a1589d852a4254312bb0e9e4b0343ab6aea32855c02c4f0f93b9fe1e18e5a3', compiled['product']['spv']['sha256'])
header = one(header, '// Temporary diagnostic: original native-default modules with lane-zero clock\n// writes. No arithmetic, output, refusal, barrier, or governor changes.',
             '// Temporary diagnostic: original native-default arithmetic with selected\n// division and phase observations. Original outputs, guards and barriers remain.')
header = one(header, '// 16 phase cuts, followed by seven stages of 117 bounded layer records.\n// Each layer retains its node range, a three-cut prefix mask and uint2 clocks.',
             '// 16 phase cuts, then seven stages of twelve 56-word division records.\n// The existing allocation stride includes cleared unused row padding.')
# Invert the only module/comment changes to bind the reusable host implementation.
recovered = re.sub(r'inline constexpr std::uint32_t kTransportClock\[\] = \{.*?\n\};',
                   re.search(r'inline constexpr std::uint32_t kTransportClock\[\] = \{.*?\n\};', source, re.S)[0], header, count=1, flags=re.S)
recovered = one(recovered, compiled['product']['spv']['sha256'], 'e3a1589d852a4254312bb0e9e4b0343ab6aea32855c02c4f0f93b9fe1e18e5a3')
recovered = one(recovered, '// Temporary diagnostic: original native-default arithmetic with selected\n// division and phase observations. Original outputs, guards and barriers remain.',
                '// Temporary diagnostic: original native-default modules with lane-zero clock\n// writes. No arithmetic, output, refusal, barrier, or governor changes.')
recovered = one(recovered, '// 16 phase cuts, then seven stages of twelve 56-word division records.\n// The existing allocation stride includes cleared unused row padding.',
                '// 16 phase cuts, followed by seven stages of 117 bounded layer records.\n// Each layer retains its node range, a three-cut prefix mask and uint2 clocks.')
assert recovered == source
host = WORK / 'native-host'
basis = WORK / 'native-host-basis'
assert not host.exists() and not basis.exists()
host.mkdir();basis.mkdir()
record = {'accepted': False, 'baseline_revision': BASELINE, 'ancestor_revision': ANCESTOR,
          'input_seals': {str(p): seal(p) for p in (Path(__file__), WORK / 'actual-placement-review.json',
              WORK / 'compiler.json', WORK / 'phase/retained_transport.spv')}, 'basis': {}, 'output': {},
          'reusable_header_inverse_byte_exact': True,
          'scope': 'Diagnostic source preparation only. Reuse reviewed device-clock feature chain and original three-test observer adapter from exact Git ancestor, replacing only embedded Transport observation module/comments. No checkout mutation/build/runtime. Existing modes, 8MiB actual sidecar charge, clears/captures and all original assertions retained.'}
record['input_seals'].update(placement['reviewed_seals'])
for name in NAMES:
    original = git('show', ANCESTOR + ':' + name)
    b = basis / name; b.parent.mkdir(parents=True, exist_ok=True); b.write_bytes(original)
    p = host / name; p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(header.encode() if name == NAMES[2] else original)
    record['basis'][name] = {'git_blob': git('rev-parse', ANCESTOR + ':' + name).decode().strip(), 'seal': seal(b)}
    record['output'][str(p)] = seal(p)
    if (ROOT / name).exists(): record['input_seals'][str(ROOT / name)] = seal(ROOT / name)
assert git('rev-parse', 'HEAD').decode().strip() == BASELINE and not git('status', '--porcelain')
for name, expected in record['input_seals'].items(): assert seal(Path(name)) == expected
record['accepted'] = True
(WORK / 'host-preparation.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({'output': record['output'], 'reusable_header_inverse_byte_exact': True}))
