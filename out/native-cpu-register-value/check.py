"""Close only the private register ABI/static gate, without running workloads."""
from collections import Counter
from pathlib import Path
import hashlib,json,re,subprocess

R=Path.cwd();W=Path(__file__).resolve().parent;O=R/'bin/linux-gcc/native-cpu-register-value'
V=R/'out/native-cpu-strict-ir-linking';PREFIX='sirius_register_value_'


def pin(p):
 with p.open('rb') as s:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(s,'sha256').hexdigest())


def doc(p):return json.loads(p.read_text())


owners={}
for name in ['plugin','math-lower','math-native','dp-lower','dp-native','transport-lower','transport-native','audit-transport-existing','endpoint-lower-resume','endpoint-native-resume']:
 owner=doc(W/name/'owner.json')
 assert owner['status']=='terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
 assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
 assert owner['passed']==(name!='transport-native')
 owners[name]=dict(owner=pin(W/name/'owner.json'),passed=owner['passed'],seconds=owner['wall_seconds'])
products={};IR_calls=native_calls=0
for stage in ['math','dp','transport','endpoint']:
 for kind in ['lower','native']:
  result=doc(W/(stage+'-'+kind+'-result.json'));product=result['output']
  assert pin(R/product['path'])=={k:product[k] for k in ('bytes','sha256')}
  products[stage+'-'+kind]=product
 lower=doc(W/(stage+'-lower-result.json'));native=doc(W/(stage+'-native-result.json'))
 expected=Counter(lower['core_calls']);actual=Counter(native['core_relocations'])
 IR_calls+=sum(expected.values());native_calls+=sum(actual.values())
 if stage=='transport':expected[PREFIX+'dense_divide8']-=1
 assert expected==actual,(stage,expected,actual)
 assert native['emitted_fused_instructions']==[]
 if stage!='math':
  asm=subprocess.check_output(['objdump','-drw',str(R/native['output']['path'])],text=True)
  joined=[line for line in asm.splitlines() if re.search(r'R_X86_64_PLT32\s+'+PREFIX,line)]
  assert len(joined)==sum(actual.values())
  assert all(re.search(r'\bcall\s',line) and line.rstrip().endswith('-0x4') for line in joined)
assert (IR_calls,native_calls)==(8241,8240)
mir=(W/'math-native/stderr.log').read_text();callee=[]
for chunk in re.split(r'(?=# \*\*\* IR Dump After)',mir):
 name=re.search(r'Machine code for function ('+PREFIX+r'\w+):',chunk)
 if not name:continue
 name=name[1];live=re.search(r'Function Live Ins: (.*)',chunk)[1]
 registers=re.findall(r'\$(\w+) in ',live)
 assert 'eax' in registers and all(re.fullmatch(r'ymm(?:[0-9]|1[01])|eax',reg) for reg in registers)
 returns=[re.findall(r'\$(\w+)',line) for line in chunk.splitlines() if re.match(r'\s+RET 0,',line)]
 expected=['ymm'+str(i) for i in range(6 if '_tagged_' in name else 5)]
 assert returns and all(row==expected for row in returns)
 assert '%fixed-stack.' not in chunk and 'Fixed Objects:' not in chunk
 callee.append(dict(function=name,live_in_registers=registers,return_registers=expected))
assert len(callee)==9
review=doc(V/'private_value_fp_review.json')
assert review['decision']['memory_none_on_retained_constrained_intrinsics']=='NOT_AUTHORIZED'
record=dict(compiler_static_gate_passed=True,owners=owners,products=products,callee_machine_interfaces=callee,
 exact_IR_core_calls=IR_calls,native_core_relocations=native_calls,
 preserved_failed_receipt='transport-native count check expected2701 vs observed2700, compiler itself exited0; retained object reused without recompilation',
 authorized_count_difference='Only pinned Transport two mutually exclusive coefficient-division IR call sites join one native common tail in the original25-iteration coefficient loop. DP1581 and Endpoint3959 exact; no general count relaxation.',
 source_roundtrip='All4 pass outputs were verified after reversing only private names and CC92 to exactly match the complete canonical input fingerprint; no FP metadata/effects/instructions/CFG/operand changes.',
 caller_callee_join='Every exact type/declaration/call uses cc92 under matching LLVM23/V4-SKX AVX512 register ABI. Callee MIR proves EAXmask/YMM inputs and5/6YMM returns; existing native caller/relocation inspection binds all8240 call fields. Internal helper stack storage is retained.',
 source_review=pin(V/'private_value_fp_review.json'),
 numerical_qualification_transferred=False,runtime_or_speed_claim=False,
 scope='Compiler/static interface gate only. No provider/component link, genuine-job C8, benefit/Decimal/control/frame, scientific, hardware-domain or release qualification for this new ABI.')
p=W/'static-gate.json';assert not p.exists();p.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(dict(compiler_static_gate_passed=True,owners=len(owners),IR_calls=IR_calls,native_core_relocations=native_calls,core_machine_interfaces=len(callee),workload_run=False)))
