"""Inspect actual guarded v3/AVX2 objects; reject FMA and further extensions."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-simd-rows/inspect_objects.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '2f4572906c6f266ba0b20fe553ebb0844cf08c8c8e933ea3902554206e5e212e'
source = basis.read_text()
raw_decode = '''    raw = subprocess.check_output(['objdump', '-d', '--insn-width=16', str(path)], text=True)
    normalized, evex = [], []
    for line in raw.splitlines():
        fields = line.split('\\t')
        if len(fields) >= 3 and re.match(r'^\\s*[0-9a-f]+:$', fields[0]):
            encoded = fields[1].split()
            assert encoded and all(re.fullmatch(r'[0-9a-f]{2}', byte) for byte in encoded)
            if encoded[0] == '62':
                evex.append(line.strip())
            normalized.append(fields[0] + '\\t' + '\\t'.join(fields[2:]))
        else:
            normalized.append(line)
    asm = '\\n'.join(normalized)
    assert not evex, ('AVX512/EVEX instruction encoding', evex)
'''
extensions = '''    bmi = re.compile(r'^(?:andn|bextr|blsi|blsmsk|blsr|bzhi|mulx|pdep|pext|rorx|sarx|shlx|shrx)(?:[bwlq])?$')
    v3_extensions = [op for op in opcodes if bmi.fullmatch(op) or
        op in ('lzcnt', 'lzcntl', 'lzcntq', 'tzcnt', 'tzcntl', 'tzcntq',
               'popcnt', 'popcntl', 'popcntq', 'movbe', 'movbel', 'movbeq',
               'vcvtph2ps', 'vcvtps2ph', 'crc32')]
    extra_isa = [op for op in opcodes if
        op in ('rdrand', 'rdrandl', 'rdrandq', 'rdseed', 'rdseedl', 'rdseedq',
               'adcx', 'adox', 'rdpid', 'rdtscp') or
        op.startswith(('aes', 'vaes', 'pclmul', 'vpclmul', 'sha', 'vgf2p8', 'vpdp', 'vpmadd52', 'tile'))]
    assert not extra_isa, ('ISA requirements beyond guarded v3/AVX2', extra_isa)
    assert not re.search(r'%(?:zmm\\d+|k[0-7]|mm[0-7]|st)\\b', asm)
'''
changes = [
    ("    asm = subprocess.check_output(['objdump', '-d', '--no-show-raw-insn', str(path)], text=True)\n", raw_decode),
    ('    assert not forbidden, forbidden\n', '    assert not forbidden, forbidden\n' + extensions),
    ("'ldmxcsr', 'stmxcsr', 'fldcw', 'fnstcw'", "'ldmxcsr', 'stmxcsr', 'vldmxcsr', 'vstmxcsr', 'fldcw', 'fnstcw'"),
    ("if op in ('addps', 'subps', 'mulps', 'divps', 'sqrtps')", "if op in ('vaddps', 'vsubps', 'vmulps', 'vdivps', 'vsqrtps')"),
    ("assert vector.get('mulps', 0) and vector.get('addps', 0)", "assert vector.get('vmulps', 0) and vector.get('vaddps', 0)"),
    ("'forbidden_fma_fp_control_x87_found': forbidden,", "'forbidden_fma_fp_control_x87_found': forbidden, 'guarded_v3_extensions': v3_extensions, 'forbidden_non_v3_extensions': extra_isa, 'EVEX_encodings': len(evex), 'instruction_mnemonics': sorted(opcodes),"),
]
for before, after in changes:
    assert source.count(before) == 1, before
    source = source.replace(before, after)
source = source.replace('ispc-uniform-compile-result.json', 'ispc-zero-compile-result.json').replace('compile/ispc-uniform/owner.json', 'compile/ispc-zero/owner.json').replace('object-code-inspection-uniform.json', 'object-code-inspection-zero.json')
exec(compile(source, str(basis), 'exec'))
