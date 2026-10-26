"""Reject extra ISA/FMA and verify the raw-v4 link closure before execution.

The explicit ISA additions below are AVX/AVX2/base-x86 or AVX512 F/CD/DQ/BW/VL,
all covered by the sealed OS-aware x86-64-v4 guard. A new mnemonic requires a
fresh source review; the skx CPU option is not blanket permission for its extras.
"""
from pathlib import Path
from collections import Counter
import hashlib
import json
import re
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
O = R / 'bin/linux-gcc/native-cpu-v4-raw-packets/objects'
sys.dont_write_bytecode = True
sys.path.insert(0, str(R / 'out/native-cpu-field-compatibility'))
from elf_proof import Elf

ADDITIONS = set('''
cmove kandb kandnb kmovb kmovd kmovq kmovw knotb korb kortestb kshiftrb
ktestb kxnorb kxorb vblendmps vcmpgeps vcmpgtps vextracti32x4 vinsertf128
vinserti32x4 vmovdqa32 vmovdqa64 vmovdqu32 vmovdqu8 vpandd vpblendmd
vpbroadcastq vpcmpgtw vpcmpleud vpcmpltd vpcmpltud vpcmpneqd vpcmpnleud
vpcmpnltud vplzcntd vpmovd2m vpmovm2d vpmovm2w vpmovw2m vpord
vpscatterqd vpternlogd vptestmb vptestmd vptestnmb vptestnmd vpxord
'''.split())


def doc(p):
    return json.loads(p.read_text())


def pin(p):
    with p.open('rb') as stream:
        return dict(bytes=p.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def listing(p):
    result = []
    text = subprocess.check_output(['objdump', '-d', '-M', 'intel',
        '--no-show-raw-insn', str(p)], text=True)
    for line in text.splitlines():
        match = re.match(r'^\s*[0-9a-f]+:\s+(\S+)(?:\s+(.*))?$', line)
        if not match:
            continue
        op, operands = match[1], match[2] or ''
        if op in ('{evex}', '{vex}'):
            op, _, operands = operands.partition(' ')
        assert op not in ('(bad)', '.byte')
        assert not op.startswith(('vfm', 'vfnm')), op
        result.append((op, operands))
    assert result
    return result


proof = doc(W / 'derivation.json')
built = doc(W / 'build-result.json')
assert built['passed'] and built['fixture_copy_verification_completed']
owner = doc(W / 'finish-build/owner.json')
assert owner['passed'] and owner['status'] == 'terminal'
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
for key in ('source_inputs', 'derived'):
    for name, identity in proof[key].items():
        assert pin(R / name) == identity, name
for row in built['products'] + proof['reused_objects']:
    assert pin(R / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}

rows = []
for stage in ('transport', 'endpoint', 'dopri_phase'):
    original = R / 'bin/linux-gcc/native-cpu-local-identities/objects' / (stage + '.o')
    candidate = O / (stage + '.o')
    old, new = listing(original), listing(candidate)
    old_names, new_names = {x for x, _ in old}, {x for x, _ in new}
    assert not new_names - old_names - ADDITIONS, sorted(new_names - old_names - ADDITIONS)
    elf = Elf(candidate)
    exports = {s['name'] for s in elf.symbols if s['section'] and s['info'] >> 4 in (1, 2)}
    assert len(exports) == 2
    assert all(name.startswith('sirius_simd_v4_' + stage + '_packet') for name in exports)
    undefined = {s['name'] for s in elf.symbols if not s['section'] and s['name']}
    allowed = {'memset', 'sirius_simd_fmaxf'} | {
        'sirius_native_v4_' + kind + '_' + op + '8'
        for kind in ('dense', 'tagged')
        for op in ('add', 'subtract', 'multiply', 'divide', 'sqrt')}
    assert undefined <= allowed and {'memset', 'sirius_simd_fmaxf'} <= undefined
    assert any('ymm16' in operands or re.search(r'\b[xyz]mm(?:[12][0-9]|3[01])\b', operands)
               for _, operands in new)
    rows.append(dict(stage=stage, original=dict(path=str(original.relative_to(R)), **pin(original)),
        candidate=dict(path=str(candidate.relative_to(R)), **pin(candidate)),
        original_instructions=len(old), candidate_instructions=len(new),
        original_stack_address_operands=sum('rsp' in x or 'rbp' in x for _, x in old),
        candidate_stack_address_operands=sum('rsp' in x or 'rbp' in x for _, x in new),
        additions=sorted(new_names - old_names), all_mnemonics=sorted(new_names),
        exports=sorted(exports), undefined=sorted(undefined), no_FMA=True,
        direct_optional_helpers_only=True, upper_vector_registers_used=True))

adapter = O / 'rows_adapter.o'
assert not any(op.startswith(('v', 'k')) for op, _ in listing(adapter))
elf = Elf(adapter)
for stage in ('transport', 'endpoint', 'dopri_phase'):
    name = 'sirius_native_' + stage + '_rows_try_packet'
    symbol, = [s for s in elf.symbols if s['name'] == name and s['section']]
    targets = {target['name'] for section, offset, _, _, target in elf.relocations
               if section == symbol['section'] and symbol['value'] <= offset < symbol['value'] + symbol['size']}
    assert {'sirius_native_simd_supported', 'sirius_native_v4_supported',
        'sirius_native_' + stage + '_try_packet', 'sirius_native_' + stage + '_admit_packet',
        'sirius_simd_' + stage + '_packet', 'sirius_simd_v4_' + stage + '_packet'} <= targets

result = dict(passed=True, objects=rows, baseline_adapter=dict(path=str(adapter.relative_to(R)), **pin(adapter)),
    no_optional_ISA_in_adapter=True, original_scalar_v3_and_admission_edges_present=True,
    all_emitted_instruction_mnemonics_reviewed_against_original_and_explicit_v4_additions=True,
    arithmetic_helpers_unchanged=True,
    limitations='Static whole-object code/closure check, no elapsed-speed/spill/memory-traffic claim. Stack-address counts include LEA and frame accesses. No provider/workload/scientific/release acceptance.')
with (W / 'code-proof.json').open('x') as stream:
    stream.write(json.dumps(result, indent=2) + '\n')
print(json.dumps(dict(passed=True, raw_objects=3, baseline_adapter=True,
                     no_FMA=True, direct_guarded_optional_helpers=True)))
