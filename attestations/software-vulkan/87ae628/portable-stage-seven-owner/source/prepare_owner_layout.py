"""Generate optional owner layout/maps from the accepted program authority.

The immutable device prefix copies existing program arrays at provisioning;
the C++ header emits only maps, avoiding duplicate static program payloads.
No literal table aggregates are emitted into Slang.
"""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, re

ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
SOURCE=ROOT/'src/sirius/kernels/retained_program.py'
HEADER=ROOT/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
assert hashlib.sha256(HEADER.read_bytes()).hexdigest()=='3f7aca3372bb872934e978fb0270f2375ce8d2361555b35c8b891e42e07161cc'
namespace={'__name__':'owner_layout_producer'}
exec(compile(SOURCE.read_bytes(),str(SOURCE),'exec'),namespace)
families=namespace['build_stage_seven_owner_plan']()
arrays={}
for count,name,body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',HEADER.read_text(),re.S):
    words=tuple(int(v[:-1],0) for v in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',body))
    assert len(words)==int(count) and name not in arrays;arrays[name]=words
assert len(arrays)==40
def encode(p):
    return tuple([p['instructions'],p['registers'],len(p['outputs'])]+p['outputs']+p['operations']+p['layer_offsets'])
references={};table_offsets={};maps={}
for stage in ('transport','endpoint'):
    reference=()
    for family in ('general','schwarzschild'):
        program=families[family]['programs'][stage]
        table_offsets[(family,stage)]=len(reference)
        reference+=encode(program)
        index=0 if stage=='transport' else 1
        slots=[0xffffffff]*program['instructions']
        for slot,pair in enumerate(families[family]['writes']):
            ordinal=pair[index];assert slots[ordinal]==0xffffffff;slots[ordinal]=slot
        maps[(family,stage)]=tuple(slots)
    assert reference==arrays['k'+stage.title()+'Program']
    references[stage]=reference

# Mutable control is disjoint from the once-verified immutable prefix and
# per-row storage. High-row authorization occupies its own row-header word.
constants={'ControlWords':16,'RowHeaderWords':8,'KeyCount':9,'ValueWords':5,
           'MaximumValues':max(len(f['writes']) for f in families.values())}
cursor=constants['ControlWords']
for stage in ('transport','endpoint'):
    constants[stage.title()+'ReferenceOffset']=cursor
    constants[stage.title()+'ReferenceWords']=len(references[stage]);cursor+=len(references[stage])
for family in ('general','schwarzschild'):
    for stage in ('transport','endpoint'):
        label=family.title()+stage.title()
        constants[label+'MapOffset']=cursor
        constants[label+'Instructions']=len(maps[(family,stage)])
        constants[label+'TableOffset']=table_offsets[(family,stage)]
        constants[label+'TableWords']=len(encode(families[family]['programs'][stage]))
        cursor+=len(maps[(family,stage)])
    constants[family.title()+'Values']=len(families[family]['writes'])
constants['PrefixWords']=cursor
constants['KeyOffset']=constants['RowHeaderWords']
constants['ValuesOffset']=constants['KeyOffset']+constants['KeyCount']*constants['ValueWords']
constants['RowWords']=constants['ValuesOffset']+constants['MaximumValues']*constants['ValueWords']
assert constants['MaximumValues']==71 and constants['SchwarzschildValues']==26
assert constants['RowWords']==408

cpp=['// Generated optional owner layout/maps; original program arrays remain authoritative.',
     '#pragma once','#include <array>','#include <cstddef>','#include <cstdint>',
     'namespace sirius::backend::retained_owner {']
slang=['// Generated scalar owner layout only; references/maps reside in the device buffer.']
for name,value in constants.items():
    cpp.append(f'inline constexpr std::uint32_t k{name} = {value}u;')
    slang.append(f'static const uint S7{name} = {value}u;')
for (family,stage),values in maps.items():
    name=family.title()+stage.title()+'Map'
    cpp.append(f'inline constexpr std::array<std::uint32_t, {len(values)}> k{name}{{{{')
    for start in range(0,len(values),16):cpp.append(','.join(str(v)+'u' for v in values[start:start+16])+',')
    cpp.append('}};')
cpp.append('} // namespace sirius::backend::retained_owner')
destination=WORK/'owner-layout';assert not destination.exists();destination.mkdir()
for name,lines in (('retained_owner_layout.h',cpp),('retained_owner_layout.slang',slang)):
    (destination/name).write_text('\n'.join(lines)+'\n')
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
result={'accepted':True,'source':seal(SOURCE),'protected_header':seal(HEADER),'constants':constants,
        'maps':{family+':'+stage:list(values) for (family,stage),values in maps.items()},
        'generated':{p.name:seal(p) for p in destination.iterdir()},
        'scope':'Generated immutable-prefix layout/maps only. Program references are copied from original arrays at future provisioning, not duplicated as C++ or Slang literal tables. No device allocation, compiled owner module, replay or qualification.'}
(WORK/'owner-layout.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'prefix_bytes':4*constants['PrefixWords'],'row_bytes':4*constants['RowWords'],'files':result['generated']},indent=2))
