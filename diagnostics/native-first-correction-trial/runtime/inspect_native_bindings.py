"""Check complete built payloads in bounded read-only PE sections; no execution."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, re, struct, subprocess
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
def read(path):return json.loads(path.read_text(encoding='utf-8-sig'))
def seal(path):
    data=path.read_bytes()
    return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
def git(*args):return subprocess.check_output(['git',*args],text=True,timeout=30).strip()
def readonly_sections(data):
    assert len(data)>=64 and data[:2]==b'MZ'
    pe=struct.unpack_from('<I',data,60)[0]
    assert pe+24<=len(data) and data[pe:pe+4]==b'PE\0\0'
    machine,count,_,_,_,optional,_=struct.unpack_from('<HHIIIHH',data,pe+4)
    assert machine==0x8664 and 0<count<=96 and optional>=112
    assert pe+24+optional+40*count<=len(data)
    assert struct.unpack_from('<H',data,pe+24)[0]==0x20b
    selected=[]
    for index in range(count):
        offset=pe+24+optional+40*index
        name=data[offset:offset+8].split(b'\0',1)[0].decode('ascii')
        size,position=struct.unpack_from('<II',data,offset+16)
        flags=struct.unpack_from('<I',data,offset+36)[0]
        assert not size or position+size<=len(data),(name,position,size)
        if size and flags&0x40 and flags&0x40000000 and not flags&0xa0000000:
            selected.append({'name':name,'file_offset':position,'bytes':size,
                             'flags':flags,'payload':data[position:position+size]})
    assert selected
    return selected
assert ROOT.resolve()==Path('/home/astra/.project/Sirius')
assert WORK==ROOT/'out/native-division-first-correction-screen/runtime'
plan=read(WORK/'execution-plan.json');head=plan['source_revision']
assert git('rev-parse','HEAD')==head==git('rev-parse','@{u}') and not git('status','--porcelain')
review=read(WORK/'native-binding-reader-review.json')
assert review['accepted'] and review['reader_seal']==seal(Path(__file__))
compiled=read(WORK/'compiled-review.json');assert compiled['accepted'] and compiled['source_revision']==head
assert compiled['candidate_header']==seal(WORK/'candidate/retained_kernels.h')
paths={Path(__file__),WORK/'native-binding-reader-review.json',WORK/'execution-plan.json',
       WORK/'compiled-review.json',WORK/'product-inspection.json',WORK/'candidate/retained_kernels.h'}
expected={}
for ledger in (plan['input_seals'],compiled['input_seals']):
    for name,record in ledger.items():
        path=(Path(name) if Path(name).is_absolute() else ROOT/name).resolve()
        if path in expected:assert expected[path]==record,name
        expected[path]=record;paths.add(path)
for path,record in expected.items():assert seal(path)==record,path
producer=plan['producers']['candidate'];stage=ROOT/producer['build_root']
gate_path=stage/'generated/sirius/native_build_gate.json';gate=read(gate_path)
assert gate['source']=={'clean':True,'revision':head}
assert seal(gate_path)['sha256']==producer['gate_sha256']
paths.add(gate_path)
arrays={};header=(WORK/'candidate/retained_kernels.h').read_text()
for count,name,text in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',header,re.S):
    words=[int(word[:-1],0) for word in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',text)]
    assert len(words)==int(count) and name not in arrays
    arrays[name]=struct.pack('<'+'I'*len(words),*words)
assert len(arrays)==40
products=read(WORK/'product-inspection.json');assert products['accepted']
assert {name:{'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()} for name,data in arrays.items()}==products['arrays']
consumers={}
for name in ('sirius','sirius_backend_tests'):
    record=gate['tested_artifacts'][name]
    assert record['root']=='build'
    relative=Path(record['path']);assert not relative.is_absolute() and '..' not in relative.parts
    path=(stage/relative).resolve();assert path.is_relative_to(stage.resolve())
    assert seal(path)=={k:record[k] for k in ('bytes','sha256')}
    paths.add(path);consumers[name]={'path':str(path),'record':record}
inputs={str(path):seal(path) for path in sorted(paths)}
for name,consumer in consumers.items():
    data=Path(consumer['path']).read_bytes();sections=readonly_sections(data);bindings={}
    for array,payload in arrays.items():
        matches=[]
        for section in sections:
            offset=section['payload'].find(payload)
            if offset>=0:matches.append({'section':section['name'],'file_offset':section['file_offset']+offset,
                                         'bytes':len(payload),'section_flags':section['flags']})
        assert matches,(name,array)
        bindings[array]=matches
    consumer['whole_readonly_payload_bindings']=bindings
assert {str(path):seal(path) for path in sorted(paths)}==inputs
assert git('rev-parse','HEAD')==head==git('rev-parse','@{u}') and not git('status','--porcelain')
result={'accepted':True,'source_revision':head,'reader_seal':seal(Path(__file__)),'input_seals':inputs,
        'whole_generated_arrays':40,'consumer_count':2,'whole_readonly_bindings':80,'consumers':consumers,
        'scope':'Actual source-bound native PE app/backend read-only initialized non-executable payload membership; immutable plan/product/provider and pre/post input joins. No native execution, numerical result, runtime benefit, rendered frame or full scientific/release qualification.'}
target=WORK/'native-bindings.json';assert not target.exists();target.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'accepted':True,'consumers':2,'whole_readonly_bindings':80,'inputs':len(inputs)}))
