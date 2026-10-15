"""Recheck saved producer seals and forty array byte spans in restored ELFs."""
from pathlib import Path
import hashlib,json,struct,subprocess

R = Path.cwd().resolve();Q=Path(__file__).resolve().parent;W=Q/'benefit'
p=json.loads((Q/'preparation.json').read_text());proof=json.loads((Q/'producer-acceptance.json').read_text())
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==p['source_revision']
assert not subprocess.check_output(['git','status','--porcelain'],text=True)
def identity(path):
    data=path.read_bytes();return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
records={};owner_records={}
for mode,item in p['executables'].items():
    path=R/item['path'];assert identity(path)=={key:item[key] for key in ['bytes','sha256']}
    data=path.read_bytes();h=struct.unpack_from('<16sHHIQQQIHHHHHH',data)
    assert h[0][:6]==b'\x7fELF\x02\x01' and h[2]==62
    sections=[struct.unpack_from('<IIQQQQIIQQ',data,h[6]+i*h[11]) for i in range(h[12])]
    found={};prefix='_ZN6sirius7backend16retained_program'
    expected=dict(proof['actual_ELF_arrays'][mode])
    expected.update({'owner::'+name:pin for name,pin in proof['actual_ELF_owner_arrays'][mode].items()})
    names={prefix+str(len(name))+name+'E':name for name in proof['actual_ELF_arrays'][mode]}
    names.update({'_ZN6sirius7backend14retained_owner'+str(len(name))+name+'E':'owner::'+name for name in proof['actual_ELF_owner_arrays'][mode]})
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
    assert len(found)==42
    records[mode]={name:pin for name,pin in found.items() if not name.startswith('owner::')}
    owner_records[mode]={name.removeprefix('owner::'):pin for name,pin in found.items() if name.startswith('owner::')}
    assert len(records[mode])==40 and len(owner_records[mode])==2
for name,pin in p['current_byte_exact_dependency_joins'].items():
    assert identity(R/name)=={key:pin[key] for key in ['bytes','sha256']},name
assert sorted(name for name in proof['actual_ELF_arrays']['baseline'] if records['baseline'][name]!=records['candidate'][name])==p['changed_arrays']
assert len(p['changed_arrays'])==17
assert owner_records['baseline']==owner_records['candidate']
for name,pin in p['preserved_archive_joins'].items():
    assert identity(Path(name))==pin,name
assert p['prior_manifest_sha256']==json.loads((Q/'prerequisite-pin.json').read_text())['manifest_sha256']
assert identity(R/json.loads((Q/'prerequisite-pin.json').read_text())['archive']/'manifest.json')['sha256']==p['prior_manifest_sha256']
assert json.loads((Q/'producer-feature-acceptance.json').read_text())['passed']
result={'passed':True,'actual_ELF_owner_arrays':owner_records,'actual_ELF_arrays':records,'executables':p['executables'],
        'same_accepted_pair_as_original_controls':True,'changed_arrays':p['changed_arrays'],'owner_arrays_unchanged':True,'prerequisite_manifest_sha256':p['prior_manifest_sha256'],'current_dependency_joins':len(p['current_byte_exact_dependency_joins']),
        'scope':'Exact restored producer and array/data/current declared dependency evidence. Existing compiler/source/ABI proof reused; no recompile or new performance observation.'}
(Q/'producer-reuse.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'passed':True,'ELF_arrays_per_producer':40,'changed':17,'dependency_joins':result['current_dependency_joins']}))
