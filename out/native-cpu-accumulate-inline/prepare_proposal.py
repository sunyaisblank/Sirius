"""Prepare one unapplied qualifier-only proposal from sealed L ISPC sources."""
from pathlib import Path
from collections import Counter
import difflib
import hashlib
import json
import re
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
previous = root / 'out/native-cpu-simd-endpoint-regions'
q = root / 'out/native-cpu-transport-worker-occupancy/provider'
frame = root / 'out/native-cpu-transport-worker-cohort/cold-frame/frame'


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def bytes_pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def publish(name, data):
    path = work / name
    if path.exists():
        assert path.read_bytes() == data, path
    else:
        path.write_bytes(data)


receipt = json.loads((previous / 'ispc-derivation.json').read_text())
compiled = json.loads((previous / 'ispc-uniform-compile-result.json').read_text())
linked = json.loads((q / 'linked-native-proof.json').read_text())
assert receipt['instructions'] == 9380 and compiled['all_objects_compiled'] and linked['passed']
library = root / linked['library']
assert pin(library) == linked['library_identity']
guide = root / 'bin/linux-gcc/native-cpu-simd-rows/compiler/ispc.html'
acquisition = json.loads((previous / 'compiler-acquisition.json').read_text())
assert pin(guide) == {k: acquisition['extracted_members'][str(guide.relative_to(root))][k]
                      for k in ('bytes', 'sha256')}
assert 'inline</tt> is not a hint, but forces inlining.' in guide.read_text()
patch, stages, calls = '', {}, 0
inputs = {str(p.relative_to(root)): pin(p) for p in (previous / 'ispc-derivation.json',
    previous / 'ispc-uniform-compile-result.json', previous / 'compiler-acquisition.json',
    previous / 'derive_ispc.py', guide, q / 'linked-native-proof.json',
    frame / 'all-threads-observation.json', frame / 'all-threads-symbols.json')}
inputs[str(library.relative_to(root))] = linked['library_identity']
inputs['/usr/bin/objdump'] = pin(Path('/usr/bin/objdump'))
for slug, entry in receipt['stages'].items():
    path = root / entry['derived_path']
    assert pin(path) == entry['derived']
    old = path.read_text()
    header = 'static RetainedTriple_0 RTAccumulateTerm_0('
    replacement = 'static inline RetainedTriple_0 RTAccumulateTerm_0('
    assert old.count(header) == 1 and replacement not in old
    new = old.replace(header, replacement, 1)
    assert new.replace(replacement, header, 1) == old
    instruction_calls = [line for line in old.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
    assert instruction_calls == [line for line in new.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
    assert len(instruction_calls) == entry['literal_instruction_calls_preserved']
    assert len(re.findall(r'\bRTAccumulateTerm_0\s*\(', old)) == 7
    calls += len(instruction_calls)
    name = str(path.relative_to(root))
    patch += ''.join(difflib.unified_diff(old.splitlines(True), new.splitlines(True),
        fromfile='a/' + name, tofile='b/' + name, n=3))
    stages[slug] = {'basis_path': name, 'basis': pin(path), 'candidate': bytes_pin(new.encode()),
        'one_inline_qualifier_only': True, 'exact_inverse': True,
        'all_body_and_call_bytes_unchanged': True, 'source_call_sites': 6,
        'literal_instruction_calls_preserved': len(instruction_calls)}
    inputs[name] = pin(path)
assert calls == 9380
native = linked['objects'][4:7]
assert [row['object'] for row in native] == [row['object'] for row in compiled['objects']]
machine, excerpts = {}, []
for obj in native:
    stage = Path(obj['object']['path']).stem
    relevant = [f for f in obj['functions'] if any(f['symbol'].startswith(name) for name in
                ('RTAccumulateTerm_0', 'RTAccumulateProduct_0', 'RTAccumulateRoundedProduct_0', 'RTAdd_0'))]
    machine[stage] = []
    for f in relevant:
        argv = ['/usr/bin/objdump', '-d', '--no-show-raw-insn',
            '--start-address=' + str(f['linked_address']),
            '--stop-address=' + str(f['linked_address'] + f['bytes']), str(library)]
        lines = subprocess.check_output(argv, text=True, timeout=30).splitlines()
        ins = [line for line in lines if re.match(r'^\s*[0-9a-f]+:\s+\S', line)]
        term_calls = [i for i, line in enumerate(ins) if 'call' in line and '<RTAccumulateTerm_0' in line]
        item = {**{key: f[key] for key in ('symbol', 'linked_address', 'bytes')},
                'instruction_count': len(ins), 'term_call_sites': len(term_calls),
                'outgoing_calls': sum(bool(re.search(r':\s+call\b', line)) for line in ins),
                'stack_reference_instructions': sum('%rsp' in line or '%rbp' in line for line in ins)}
        if f['symbol'].startswith('RTAccumulateTerm_0'):
            assert not item['stack_reference_instructions'] and not item['outgoing_calls']
            assert sum('vzeroupper' in line for line in ins) == 1
        machine[stage].append(item)
        if stage == 'transport':
            excerpts.append('\n' + f['symbol'] + '\n')
            selected = set(range(min(7, len(ins)))) | set(range(max(0, len(ins)-10), len(ins)))
            for index in term_calls[:2]:
                selected.update(range(max(0, index-6), min(len(ins), index+17)))
            excerpts.extend(ins[i] + '\n' for i in sorted(selected))
symbols = json.loads((frame / 'all-threads-symbols.json').read_text())
object_hits = Counter()
for row in symbols['aggregate']['functions']:
    for name in row['native_objects']:
        object_hits[name] += row['samples']
proof = {
    'scope': 'Unapplied three-line ISPC source qualifier proposal only; no shader/native compiler, workload, provider, object, host, flags or Git mutation.',
    'source_inputs': inputs, 'stages': stages, 'instructions': calls,
    'patch_applied': False, 'new_native_compilations': 0,
    'unchanged_compile_flags': compiled['objects'][0]['argv'][2:-3],
    'official_semantics': {'url': 'https://ispc.github.io/ispc.html#functions-and-function-calls',
        'installed_release': acquisition['release_tag'], 'pinned_installed_guide': pin(guide),
        'meaning': 'ISPC inline forces inlining; static visibility and all function types/masks remain unchanged. No unmasked qualifier or compiler option is proposed.'},
    'mechanism': 'Term is a leaf with no own stack references. Its private varying aggregate is five32-byte vectors: input through RSI, hidden return pointer RDI, x in YMM0 and execution mask in YMM1. Its common epilogue stores five vectors then vzeroupper/ret. Product, RoundedProduct and Add callers materialize/copy these aggregates on stack between the original2/1/3 ordered term calls. Inlining exposes those private temporaries to scalar replacement and removes that function ABI boundary; optimization magnitude is unproven.',
    'source_semantics': 'All three RetainedSum32 calls, upward error update, input/output bound checks, invalid returns, masks, return fields and caller order remain byte-identical. All9380 literal instruction-call lines and every other source byte remain identical. Scalar D, F admission, I adapter, Q scheduler/ISA policy, P host, guards, fields/RHS/trials/quartics/limits and original gates are unchanged.',
    'existing_Q_machine_code': machine,
    'existing_frame_window': {'status': symbols['capture_status'], 'samples': symbols['aggregate']['samples'],
        'native_object_hits': dict(object_hits),
        'limits': 'Existing single3s USER window only. RET/skid concentration and object samples do not isolate call cost, prove whole-frame causality or predict speedup.'},
    'future_finishline': 'Before any original gate, inspect all three actual compiled callers: no remaining RTAccumulateTerm calls/sret boundary, record stack/register/code-size changes, and retain ISA/FP/noFMA/noFPcontrol checks. Original C8/all32/51/54, one ABBA22+Decimal, all8 and original frame stay required. Inlining may grow code or increase register spills; source annotation alone is not machine-code or benefit acceptance.'}
publish('accumulate-inline.patch', patch.encode())
publish('excerpts.asm', ''.join(excerpts).encode())
publish('source-proof.json', (json.dumps(proof, indent=2) + '\n').encode())
assert all(pin(root / name) == value for name, value in inputs.items())
print(json.dumps({'unapplied_source_changes': 3, 'instruction_calls_unchanged': calls,
                  'compiled_or_executed_shader': False}))
