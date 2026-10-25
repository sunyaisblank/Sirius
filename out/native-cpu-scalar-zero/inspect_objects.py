"""Inspect the actual strict-baseline scalar objects and eliminated Term calls."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

R=Path.cwd();W=Path(__file__).resolve().parent
compiled=json.loads((W/'cpp-compile-result.json').read_text())
owner=json.loads((W/'compile/cpp/owner.json').read_text())
assert compiled['all_objects_compiled'] and owner['passed'] and owner['whole_inputs_unchanged']
rows={}
for item in compiled['objects']:
    path=R/item['object']['path']
    assert hashlib.file_digest(path.open('rb'),'sha256').hexdigest()==item['object']['sha256']
    symbols=subprocess.check_output(['/usr/bin/nm','-S','--defined-only',str(path)],text=True)
    assert 'RTAccumulateTerm_0' not in symbols
    assert re.search(r' T sirius_native_'+item['stage']+'_try_packet$',symbols,re.M)
    asm=subprocess.check_output(['/usr/bin/objdump','-d','--insn-width=16',str(path)],text=True)
    opcodes={};evex=[]
    for line in asm.splitlines():
        fields=line.split('\t')
        if len(fields)>=3 and re.match(r'^\s*[0-9a-f]+:$',fields[0]):
            encoded=fields[1].split();opcode=fields[2].split()[0]
            opcodes[opcode]=opcodes.get(opcode,0)+1
            if encoded[0]=='62':evex.append(line)
    forbidden=[op for op in opcodes if op.startswith(('vfm','vfnm','v')) or
        op in ('ldmxcsr','stmxcsr','fldcw','fnstcw','fadd','faddp','fmul','fmulp','fdiv','fdivp','fsub','fsubp')]
    assert not forbidden and not evex,(item['stage'],forbidden,evex)
    assert not re.search(r'%(?:ymm\d+|zmm\d+|k[0-7]|mm[0-7]|st)\b',asm)
    rows[item['stage']]={'object':item['object'],'Term_symbol_and_ABI_calls_absent':True,
        'forbidden_FMA_AVX_FP_control_x87':forbidden,'EVEX_encodings':len(evex),
        'instruction_mnemonics':sorted(opcodes),'instruction_count':sum(opcodes.values())}
result={'passed':True,'objects':rows,'source_revision':compiled['revision'],
    'scope':'Actual scalar object/ISA inspection only. Original strict x86-64 compiler recipe retained, no AVX/FMA/FP-control/x87 instructions, Term ABI removed. This does not establish all-input math equivalence, runtime benefit, original frame, science/release or adoption.'}
(W/'object-code-inspection.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'actual_scalar_objects':3,'Term_ABI_calls':0,'forbidden_instructions':0,'runtime_executed':False}))
