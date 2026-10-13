"""Recheck saved producer seals and forty array byte spans in restored ELFs."""
from pathlib import Path
import hashlib,json,struct,subprocess

R = Path.cwd().resolve();Q=Path(__file__).resolve().parent;W=Q/'benefit'
p=json.loads((Q/'preparation.json').read_text());proof=json.loads((Q/'producer-acceptance.json').read_text())
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==p['source_revision']
assert not subprocess.check_output(['git','status','--porcelain'],text=True)
def identity(path):
    data=path.read_bytes();return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
records={}
for mode,item in p['executables'].items():
    path=R/item['path'];assert identity(path)=={key:item[key] for key in ['bytes','sha256']}
    data=path.read_bytes();h=struct.unpack_from('<16sHHIQQQIHHHHHH',data)
    assert h[0][:6]==b'\x7fELF\x02\x01' and h[2]==62
    sections=[struct.unpack_from('<IIQQQQIIQQ',data,h[6]+i*h[11]) for i in range(h[12])]
    found={};prefix='_ZN6sirius7backend16retained_program'
    expected=proof['actual_ELF_arrays'][mode]
    names={prefix+str(len(name))+name+'E':name for name in expected}
    for section in sections:
        if section[1]!=2:continue
        table=sections[section[6]];strings=data[table[4]:table[4]+table[5]]
        assert section[9]==24
        for off in range(section[4],section[4]+section[5],24):
            name,info,other,index,address,size=struct.unpack_from('<IBBHQQ',data,off)
            symbol=strings[name:strings.find(b'\0',name)].decode()
            if symbol not in names:continue
            label=names[symbol];assert label not in found
            sec=sections[index]
            assert sec[2]&2 and not(sec[2]&5) and sec[1]!=8
            assert address>=sec[3] and address+size<=sec[3]+sec[5]
            raw=data[sec[4]+address-sec[3]:sec[4]+address-sec[3]+size]
            pin={'bytes':size,'sha256':hashlib.sha256(raw).hexdigest()}
            assert pin=={key:expected[label][key] for key in ['bytes','sha256']}
            found[label]=pin
    assert len(found)==40
    records[mode]=found
for name,pin in p['current_byte_exact_dependency_joins'].items():
    assert identity(R/name)=={key:pin[key] for key in ['bytes','sha256']},name
assert sum(records['baseline'][name]!=records['candidate'][name] for name in records['baseline'])==5
result={'passed':True,'actual_ELF_arrays':records,'executables':p['executables'],
        'same_qualified_pair_as_original_controls':True,'current_dependency_joins':len(p['current_byte_exact_dependency_joins']),
        'scope':'Exact restored producer and array/data/current declared dependency evidence. Existing compiler/source/ABI proof reused; no recompile or new performance observation.'}
(Q/'producer-reuse.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'passed':True,'ELF_arrays_per_producer':40,'changed':5,'dependency_joins':result['current_dependency_joins']}))
