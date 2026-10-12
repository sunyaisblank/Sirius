"""Whole generated payload/ELF and private link authority checks."""
from pathlib import Path
import hashlib,json,re,struct
r=Path.cwd();w=Path(__file__).resolve().parent;b=r/'bin/linux-gcc/portable-integer-reciprocal';records=[]
sha=lambda data:hashlib.sha256(data).hexdigest()
for mode in ('baseline','candidate'):
 header=r/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h' if mode=='baseline' else b/'products/retained_kernels.h'
 exe=b/('original_'+mode+'_controls');data=exe.read_bytes();assert data[:4]==b'\x7fELF';joins=[]
 for size,name,body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',header.read_text(),re.S):
  words=[int(v[:-1],0) for v in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',body)];assert len(words)==int(size);payload=struct.pack('<%dI'%len(words),*words);at=data.find(payload);assert at>=0,(mode,name);assert data.find(payload,at+1)<0,(mode,name,'duplicated authority');joins.append({'array':name,'bytes':len(payload),'sha256':sha(payload),'elf_offset':at})
 assert len(joins)==40
 mapping=(b/(mode+'.map')).read_text();members=('retained_compute.cpp.o','retained_integrator.cpp.o','retained_trace_executor.cpp.o','vulkan_device.cpp.o')
 if mode=='candidate':
  for m in members:assert not re.search(r'libsirius_backend\.a\('+re.escape(m)+r'\)',mapping),m
  for obj in ('compute','integrator','trace','vulkan'):assert str(b/(obj+'.o')) in mapping
 else:
  for m in members[:3]:assert re.search(r'libsirius_backend\.a\('+re.escape(m)+r'\)',mapping),m
  assert 'libsirius_backend.a(vulkan_device.cpp.o)' not in mapping and str(b/'vulkan.o') in mapping
 records.append({'producer':mode,'executable':str(exe.relative_to(r)),'executable_sha256':sha(data),'header_sha256':sha(header.read_bytes()),'full_payload_joins':joins,'link_authority_pass':True})
result={'passed':True,'full_payloads_checked':80,'same_current_vulkan_for_both':True,'records':records};(w/'consumer-check.json').write_text(json.dumps(result,indent=2)+'\n');print({'passed':True,'payloads':80})
