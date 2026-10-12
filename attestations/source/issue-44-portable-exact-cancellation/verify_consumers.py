"""Whole generated payload/ELF and private link authority checks."""
from pathlib import Path
import hashlib,json,re,struct
r=Path.cwd();w=Path(__file__).resolve().parent;b=r/'bin/linux-gcc/portable-exact-cancellation';records=[]
sha=lambda data:hashlib.sha256(data).hexdigest()
def symbols(data):
    assert data[4:6]==b'\x02\x01'
    shoff=struct.unpack_from('<Q',data,40)[0];entsize,num=struct.unpack_from('<HH',data,58)
    assert entsize==64
    sections=[struct.unpack_from('<IIQQQQIIQQ',data,shoff+i*entsize) for i in range(num)]
    found={}
    for section in sections:
        if section[1]!=2:continue
        assert section[9]==24
        strings=sections[section[6]];strings=data[strings[4]:strings[4]+strings[5]]
        for at in range(section[4],section[4]+section[5],24):
            name,info,other,index,value,size=struct.unpack_from('<IBBHQQ',data,at)
            text=strings[name:strings.index(b'\0',name)].decode()
            match=re.fullmatch(r'_ZN6sirius7backend16retained_program\d+(k\w+)E',text)
            if not match:continue
            target=sections[index];assert target[2]&2 and not target[2]&5;offset=target[4]+value-target[3]
            assert value>=target[3] and value+size<=target[3]+target[5]
            assert match[1] not in found
            found[match[1]]=(value,offset,size)
    return found
for mode in ('baseline','candidate'):
 header=r/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h' if mode=='baseline' else b/'products/retained_kernels.h'
 exe=b/('original_'+mode+'_controls');data=exe.read_bytes();assert data[:4]==b'\x7fELF';named=symbols(data);joins=[]
 for size,name,body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',header.read_text(),re.S):
  words=[int(v[:-1],0) for v in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',body)];assert len(words)==int(size);payload=struct.pack('<%dI'%len(words),*words);address,at,symbol_size=named[name];assert symbol_size==len(payload) and data[at:at+symbol_size]==payload,(mode,name);joins.append({'array':name,'bytes':len(payload),'sha256':sha(payload),'elf_offset':at,'symbol_address':address,'named_symbol':name})
 assert len(joins)==40
 mapping=(b/(mode+'.map')).read_text();members=('retained_compute.cpp.o','retained_integrator.cpp.o','retained_trace_executor.cpp.o','vulkan_device.cpp.o')
 if mode=='candidate':
  for m in members:assert not re.search(r'libsirius_backend\.a\('+re.escape(m)+r'\)',mapping),m
  for obj in ('compute','integrator','trace','vulkan'):assert str(b/(obj+'.o')) in mapping
 else:
  for m in members[:3]:assert re.search(r'libsirius_backend\.a\('+re.escape(m)+r'\)',mapping),m
  reuse=json.loads((w/'baseline-reuse.json').read_text());identity=reuse['baseline_provider_link_input_identity'];provider=(b/'vulkan.o').read_bytes();assert len(provider)==identity['bytes'] and sha(provider)==identity['sha256'];assert 'libsirius_backend.a(vulkan_device.cpp.o)' not in mapping and reuse['baseline_provider_map_path'] in mapping
 records.append({'producer':mode,'executable':str(exe.relative_to(r)),'executable_sha256':sha(data),'header_sha256':sha(header.read_bytes()),'full_payload_joins':joins,'link_authority_pass':True})
result={'passed':True,'full_payloads_checked':80,'same_current_vulkan_for_both':True,'records':records};(w/'consumer-check.json').write_text(json.dumps(result,indent=2)+'\n');print({'passed':True,'payloads':80})
