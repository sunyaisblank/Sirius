"""Check actual C++ helper ISA and associate every function with the test ELF."""
from pathlib import Path
import json
import re
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
sys.path.insert(0, str(R / 'out/native-cpu-field-compatibility'))
from elf_proof import pin, associate

build = json.loads((W / 'cpp-build-result.json').read_text())
owner = json.loads((W / 'cpp/build/owner.json').read_text())
assert build['passed'] and owner['passed'] and owner['status'] == 'terminal'
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
row, = [item['object'] for item in build['compilations'] if item['name'] == 'retained8']
path = R / row['path']
assert pin(path) == {k: row[k] for k in ('bytes', 'sha256')}
executable = R / build['executable']['path']
assert pin(executable) == {k: build['executable'][k] for k in ('bytes', 'sha256')}
raw = subprocess.check_output(['/usr/bin/objdump', '-d', '--insn-width=16', str(path)], text=True)
mnemonics = set()
evex = []
for line in raw.splitlines():
    fields = line.split('\t')
    if len(fields) < 3 or not re.match(r'^\s*[0-9a-f]+:$', fields[0]):
        continue
    encoded = fields[1].split()
    assert encoded and all(re.fullmatch(r'[0-9a-f]{2}', byte) for byte in encoded)
    if encoded[0] == '62':
        evex.append(line)
    mnemonic = fields[2].strip().split()[0]
    mnemonics.add(mnemonic)
assert not evex
assert not re.search(r'%(?:zmm\d+|k[0-7]|mm[0-7]|st)\b', raw)
forbidden = {m for m in mnemonics if m.startswith(('vfm', 'vfnm', 'f', 'aes', 'vaes',
    'pclmul', 'vpclmul', 'sha', 'vgf2p8', 'vpdp', 'vpmadd52', 'tile')) or
    m in ('ldmxcsr','stmxcsr','vldmxcsr','vstmxcsr','rdrand','rdseed','adcx','adox',
          'rdpid','rdtscp','prefetchw','prefetchwt1','clflushopt','clwb','clzero',
          'movdiri','movdir64b','serialize','pconfig')}
assert not forbidden, forbidden
assert {'vaddps','vsubps','vmulps'} <= mnemonics
proof = associate(executable, [{**row, 'path': str(path)}])
assert proof['passed'] and all(item['all_retained_function_opcodes_match'] for item in proof['objects'])
result = {'passed': True, 'helper_object': row, 'executable': build['executable'],
          'EVEX': len(evex), 'forbidden_FMA_FP_control_x87_extensions': sorted(forbidden),
          'mnemonics': sorted(mnemonics), 'linked_code': proof,
          'scope': 'Actual guarded-v3 helper code/ISA and every retained function opcode join to the focused test ELF. No execution, numerical agreement, runtime gain, whole-renderer or scientific/release claim.'}
with (W / 'cpp-code-inspection.json').open('x') as stream:
    json.dump(result, stream, indent=2)
    stream.write('\n')
print(json.dumps({'passed': True, 'objects_joined': len(proof['objects']),
                  'retained_helper_functions': len(proof['objects'][0]['functions'])}))
