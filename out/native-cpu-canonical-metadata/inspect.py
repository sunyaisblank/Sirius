"""Read the first emitted canonical DP fast path and stop on unchanged escape."""
from pathlib import Path
from collections import Counter
import hashlib,json,re

R=Path.cwd();W=Path(__file__).resolve().parent
p=R/'bin/linux-gcc/native-cpu-canonical-metadata/dopri_phase.ll'
text=p.read_text()
f=re.search(r'^define [^\n]*@EvaluateProgram_0[^\n]*\{\n(.*?)^}',text,re.M|re.S)[1]
blocks={};current=None
for line in f.splitlines():
 match=re.match(r'^([^ ;\n]+):',line)
 if match:current=match[1];blocks[current]=[]
 else:
  assert current is not None
  blocks[current].append(line)
assert 'if_then35' in blocks and 'if_exit' in blocks
assert 'label %if_then35, label %if_exit' in '\n'.join(blocks['if_then'])
visited=set();pending=['if_then35']
while pending:
 name=pending.pop()
 if name in visited:continue
 visited.add(name)
 pending+=re.findall(r'label %([^,\s]+)','\n'.join(blocks[name]))
assert 'if_exit' not in visited
fast='\n'.join(line for name,lines in blocks.items() if name in visited for line in lines)
calls=Counter(re.findall(r'\bcall\s+.*?@([^ (]+)\(',fast))
layers={name:count for name,count in calls.items() if name.startswith('EvaluateStaticLayer_0_')}
assert len(layers)==45 and set(layers.values())=={1}
assert {int(re.search(r'EvaluateStaticLayer_0_(\d+)_0',name)[1]) for name in layers}==set(range(45))
assert not any('gather' in name or name.startswith('ReadValue_0') for name in calls)
publications={name:count for name,count in calls.items() if name.startswith('WriteStored_0')}
assert sum(publications.values())==200
context_escapes=[line.strip() for line in fast.splitlines() if re.search(r'\bcall\b.*@EvaluateStaticLayer_',line) and 'ptr nonnull %kernelContext_8' in line]
assert len(context_escapes)==45
root_temporaries=re.findall(r'^\s+%SiriusRoot\d+ = alloca %v8_varying_RetainedTriple_0, align 32$',f,re.M)
assert len(root_temporaries)==200
with p.open('rb') as stream:identity=dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(stream,'sha256').hexdigest())
result=dict(status='STOP_BEFORE_COUPLED_EMISSION_OR_RUNTIME',IR=dict(path=str(p.relative_to(R)),**identity),
 compiled=True,source_inverse_and_all9380_tuple_check_preserved=True,
 static_fast_path=dict(entry='if_then35',reachable_blocks=len(visited),fallback_entry_unreachable=True,
  dynamic_root_gathers=0,dynamic_ReadValue_calls=0,ordered_layer_calls=45,whole_context_layer_calls=45,
  root_publication_calls=200,escaping_root_temporaries=200,root_temporary_representation_bytes=32000,
  root_addressing_improved=True,layer_bank_interface_narrowed=False),
 required_gate_met=False,
 reason='Literal root addressing was emitted, but every layer still receives the complete KernelContext and all200 publication temporaries escape through WriteStored calls. This does not establish the required narrower bank/SSA interface or sufficient execution benefit. Stop before emitting T/E or new acceptance/frame work; address layer interface first.',
 scope='Source/ordinary-IR structural evidence only; no strict conversion/native object/numerical/provider/speed/scientific/release claim.')
target=W/'structural-gate.json';assert not target.exists();target.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
