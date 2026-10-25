"""Describe recorded instruction locations, without another execution/sample."""
from pathlib import Path
from collections import Counter
import hashlib,json,re,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;J=lambda p:json.loads(p.read_text())
decoded=J(W/'observe/all-threads-symbols.json');raw=J(W/'observe/all-threads-observation.json')
assert decoded['status']=='decoded' and not decoded['errors'] and raw['status']=='complete'
library=R/'bin/linux-gcc/native-cpu-local-identities/provider/libvulkan_lvp.so'
pin=lambda p:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
asm=subprocess.check_output(['objdump','-d','--no-show-raw-insn',str(library)],text=True)
instructions={int(m[1],16):m[2].strip() for line in asm.splitlines() if (m:=re.match(r'^\s*([0-9a-f]+):\s+(.+)',line))}
counts=Counter();examples=Counter();total=0
for thread in raw['threads'].values():
 for sample in thread['samples']:
  a=decoded['addresses'][hex(sample['ip'])]
  if not any(r['object'].endswith('native-cpu-value-fact-stores/objects/rows_adapter.o') for r in a.get('native_ranges',[])):continue
  total+=1;instruction=instructions.get(a['elf_vaddr']);assert instruction
  kind='stack operand' if '%rsp' in instruction or '%rbp' in instruction else 'other';counts[kind]+=1
  if kind=='stack operand':examples[instruction]+=1
flags=[line.split(':',1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('flags')]
assert flags and len(set(flags))==1
result=dict(source=J(W/'plan.json')['source'],recorded_CPP_samples=total,sampled_instruction_kinds=dict(counts),stack_operand_examples=examples.most_common(20),provider=dict(path=str(library.relative_to(R)),**pin(library)),cpu_flags=flags[0],cache_sizes={p.parent.name:p.read_text().strip() for p in Path('/sys/devices/system/cpu/cpu0/cache').glob('index*/size')},
 scope='Decode only existing USER-PC records into exact current ELF instructions. A sampled stack operand is an observed instruction address, not a measured cache miss, spill latency, causal bottleneck, whole-frame share or speedup. Source inspection distinguishes saved locals from other stack traffic. CPU flags/cache metadata describe this machine, not a portable ISA admission check. No new capture/execution or ISA/flag sweep.')
path=W/'instruction-locations.json';assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(CPP_samples=total,instruction_kinds=dict(counts))))
