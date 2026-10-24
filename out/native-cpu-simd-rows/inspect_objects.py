"""Inspect sealed object code; no provider loading, timing or shader execution."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
record = json.loads((work / 'ispc-uniform-compile-result.json').read_text())
owner = json.loads((work / 'compile/ispc-uniform/owner.json').read_text())
assert record['all_objects_compiled'] and owner['passed'] and owner['whole_inputs_unchanged']
objects = {}
for item in record['objects']:
    path = root / item['object']['path']
    assert hashlib.file_digest(path.open('rb'), 'sha256').hexdigest() == item['object']['sha256']
    symbols = subprocess.check_output(['nm', '-S', '--defined-only', str(path)], text=True)
    required = 'sirius_simd_' + item['stage'] + '_packet'
    assert re.search(r' T ' + required + r'$', symbols, re.M)
    undefined = subprocess.check_output(['nm', '-u', str(path)], text=True)
    assert 'sirius_simd_fmaxf' in undefined
    asm = subprocess.check_output(['objdump', '-d', '--no-show-raw-insn', str(path)], text=True)
    opcodes = {}
    for line in asm.splitlines():
        match = re.match(r'^\s*[0-9a-f]+:\s+(\w+)', line)
        if match:
            opcodes[match[1]] = opcodes.get(match[1], 0) + 1
    forbidden = [op for op in opcodes if op.startswith(('vfm', 'vfnm')) or
        op in ('ldmxcsr', 'stmxcsr', 'fldcw', 'fnstcw', 'fadd', 'faddp', 'fmul', 'fmulp', 'fdiv', 'fdivp', 'fsub', 'fsubp')]
    assert not forbidden, forbidden
    vector = {op: n for op, n in opcodes.items() if op in ('addps', 'subps', 'mulps', 'divps', 'sqrtps')}
    assert vector.get('mulps', 0) and vector.get('addps', 0)
    # ReadValueUniform has an explicit uniform private column. Keep the
    # ordinary varying reader separately for the original generic fallback.
    readers = []
    for function in re.split(r'\n(?=[0-9a-f]+ <)', asm):
        if re.match(r'[0-9a-f]+ <[^>]*ReadValueUniform_0', function):
            counts = {}
            for op in re.findall(r'^\s*[0-9a-f]+:\s+(\w+)', function, re.M):
                counts[op] = counts.get(op, 0) + 1
            readers.append({'symbol': function.splitlines()[0], 'opcodes': counts,
                            'instruction_count': sum(counts.values())})
    objects[item['stage']] = {'object': item['object'], 'export': required,
        'undefined_symbols': undefined.splitlines(), 'packed_float_instructions': vector,
        'forbidden_fma_fp_control_x87_found': forbidden,
        'total_instructions': sum(opcodes.values()), 'uniform_private_readers': readers}
result = {'revision': record['revision'], 'objects': objects,
    'scope': 'Actual object code inspection only. No numerical equality, runtime gain, frame or scientific/release qualification claim.'}
(work / 'object-code-inspection-uniform.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
