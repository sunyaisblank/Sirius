"""Join saved ordinary driver artifacts; static footprint, never execution cost."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,re,struct,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent;STAGE=WORK/'placement'
BASELINE='2053ef79c90752d37b0075b8137cfb9b571b6820'
OLD=ROOT/'attestations/native-vulkan/3045067/ordinary-pipeline-generated-code'
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):
    b=p.read_bytes();return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
def git(*args):return subprocess.check_output(['git',*args],text=True,timeout=30).strip()
def statistics(data):
    assert len(data)==72
    u32=lambda p:struct.unpack_from('<I',data,p)[0]
    u64=lambda p:struct.unpack_from('<Q',data,p)[0]
    values={'shader_stage_mask':u32(0),'used_vgprs':u32(8),'used_sgprs':u32(12),
            'lds_bytes_per_workgroup_limit':u32(16),'lds_usage_bytes':u64(24),'scratch_usage_bytes':u64(32),
            'physical_vgprs':u32(40),'physical_sgprs':u32(44),'available_vgprs':u32(48),'available_sgprs':u32(52),
            'workgroup_x':u32(56),'workgroup_y':u32(60),'workgroup_z':u32(64)}
    assert values['shader_stage_mask']==32
    assert [values[n] for n in ('workgroup_x','workgroup_y','workgroup_z')]==[64,1,1]
    return values
def elf_sections(data):
    assert len(data)>=64 and data[:6]==b'\x7fELF\x02\x01'
    offset=struct.unpack_from('<Q',data,40)[0]
    stride,count,names=struct.unpack_from('<HHH',data,58)
    assert stride==64 and 0<names<count and offset+stride*count<=len(data)
    headers=[struct.unpack_from('<IIQQQQIIQQ',data,offset+i*stride) for i in range(count)]
    h=headers[names];assert h[1]==3 and h[4]+h[5]<=len(data)
    strings=data[h[4]:h[4]+h[5]];sections={};inventory=[]
    for h in headers:
        assert h[0]<len(strings)
        end=strings.index(0,h[0]);name=strings[h[0]:end].decode('ascii')
        assert name not in sections
        if h[1]!=8:assert h[4]+h[5]<=len(data)
        payload=data[h[4]:h[4]+h[5]] if h[1]!=8 else b''
        sections[name]=payload
        inventory.append({'name':name,'type':h[1],'flags':h[2],'offset':h[4],'bytes':h[5]})
    selected=next(x for x in inventory if x['name']=='.text')
    assert selected['type']==1 and selected['flags']==6
    assert '.AMDGPU.comment.llvmir' in sections and sections['.AMDGPU.comment.llvmir']
    return sections,inventory
def inventory(binary,assembly,stats):
    sections,metadata=elf_sections(binary);code=sections['.text'];ir=sections['.AMDGPU.comment.llvmir']
    assert assembly[-1]==0 and b'\0' not in assembly[:-1]
    lines=assembly[:-1].decode('utf-8').splitlines();printed=bytearray();instructions=[];mnemonics={}
    for number,line in enumerate(lines,1):
        match=re.search(r';\s*((?:[0-9a-f]{8}\s*)+)$',line)
        if not match:continue
        text=line.strip().split(';')[0].strip();assert text
        instructions.append({'line':number,'text_offset':len(printed),'instruction':text})
        for word in match[1].split():printed.extend(struct.pack('<I',int(word,16)))
        mnemonic=text.split()[0];mnemonics[mnemonic]=mnemonics.get(mnemonic,0)+1
    assert instructions and len(printed)<=len(code) and bytes(printed)==code[:len(printed)]
    assert not any('MSG_RTN_GET_REALTIME' in i['instruction'] for i in instructions)
    ir_text=ir.decode('utf-8');ir_lines=ir_text.splitlines()
    phis=[{'line':i+1,'text':line.strip()} for i,line in enumerate(ir_lines)
          if re.search(r'%q\..*sroa\.32.* = phi float',line)]
    return {'statistics':statistics(stats),'ELF_sections':metadata,'ELF_text_bytes':len(code),
            'printed_static_instruction_count':len(instructions),'printed_instruction_bytes':len(printed),
            'all_printed_bytes_match_ELF_text_prefix':True,'unprinted_ELF_tail_bytes':len(code)-len(printed),
            'mnemonics':mnemonics,'lexical_barrier_instruction_count':sum('s_barrier' in i['instruction'] for i in instructions),
            'IR_bytes':len(ir),'IR_sha256':hashlib.sha256(ir).hexdigest(),'IR_compiler_pass_provenance_verified':False,
            'q_error_phi_lexical_matches':phis,
            'phi_match_scope':'Lexical marker inventory only; absence/presence does not prove source-to-IR specialization or exclusive cost.'},code,ir

assert ROOT.resolve()==Path('/home/astra/.project/Sirius') and WORK==ROOT/'out/native-division-first-correction-screen'
assert git('rev-parse','HEAD')==BASELINE==git('rev-parse','@{u}') and not git('status','--porcelain')
reader_review=read(WORK/'reader-source-review.json')
assert reader_review['accepted'] and reader_review['reader_seal']==seal(Path(__file__))
plan=read(STAGE/'plan.json');bridge=read(STAGE/'placement-query-bridge.json');owner=read(STAGE/'run/owner.json')
receipt=read(STAGE/'clock-placement.json');retirement=read(STAGE/'retirement-run.json');census=read(STAGE/'native-retirement.json')
assert plan['source_revision']==bridge['source_revision']==BASELINE
assert bridge['accepted'] and not bridge['errors'] and not bridge['outer_stop']
assert owner['status']=='completed' and owner['returncode']==0 and not owner['outer_stop']
assert retirement['accepted'] and retirement['inputs_unchanged'] and census['pass'] and not census['matches']
assert receipt['completed'] and receipt['logical_devices_created']==receipt['pipelines_created']==1 and receipt['dispatches']==0
assert receipt['descriptor_bindings']==[0,1] and not receipt['shader_device_clock_enabled'] and not receipt['shader_subgroup_clock_enabled']
assert receipt['shader_float64_enabled'] and receipt['shader_fma32_enabled'] and not receipt['shader_fma16_enabled'] and not receipt['shader_fma64_enabled']
assert receipt['pipeline_cache_flags']==receipt['pipeline_cache_initial_bytes']==0 and not receipt['driver_cold_verified']
assert receipt['enabled_extensions']==['VK_KHR_shader_fma','VK_AMD_shader_info']
assert receipt['device']=='AMD Radeon 780M Graphics' and receipt['driver']=='26.8.1 (LLPC)'
providers={p['name'].lower():p for p in plan['providers']}
loaded={p['name'].lower():p for p in receipt['loaded_modules']}
assert len(receipt['loaded_modules'])==len(providers)==2 and set(loaded)==set(providers)
for name,p in loaded.items():
    expected=providers[name];assert p['bytes']==expected['bytes'] and p['sha256']==expected['sha256']
    assert p['path'].replace('\\\\','\\').lower()==expected['path'].replace('\\\\','\\').lower()
assert len(receipt['observations'])==1
observation=receipt['observations'][0]
assert observation['stage']=='Transport' and observation['pipeline_flags']==observation['stage_flags']==observation['pipeline_creation_result']==0
candidate=STAGE/'retained_transport.spv';assert seal(candidate)=={'bytes':observation['module_bytes'],'sha256':observation['module_sha256']}
assert candidate.read_bytes()==(WORK/'products/transport.spv').read_bytes()
queries={q['type']:q for q in observation['queries']};assert set(queries)=={0,1,2} and len(observation['queries'])==3
for kind,q in queries.items():
    assert q['available'] and q['captured'] and q['size_query_result']==q['complete_query_result']==0
    assert q['file']==f'transport.info{kind}' and 0<q['returned_bytes']<=q['capacity_bytes']<=8388608
    assert seal(STAGE/'query-payloads'/q['file'])=={'bytes':q['returned_bytes'],'sha256':q['sha256']}
old=read(OLD/'analysis.json');original=old['Transport']
for name in ('statistics','binary','disassembly'):
    row=original[name];assert seal(OLD/row['file'])=={k:row[k] for k in ('bytes','sha256')}
paths={Path(__file__),WORK/'reader-source-review.json',STAGE/'plan.json',STAGE/'source-review.json',STAGE/'semantics-review.json',
       STAGE/'placement-query-bridge.json',STAGE/'run/owner.json',STAGE/'owner-temp/bootstrap.json',STAGE/'clock-placement.json',
       STAGE/'run/driver-modules.json',STAGE/'retirement-run.json',STAGE/'native-retirement.json',
       WORK/'compiled-review.json',WORK/'products/transport.spv',OLD/'analysis.json',OLD/'source-code-join-bindings.json',
       OLD/'instruction-encoding-join.json',OLD/'generated-code/transport.llvmir'}
paths.update(STAGE/'query-payloads'/q['file'] for q in queries.values())
paths.update(OLD/original[name]['file'] for name in ('statistics','binary','disassembly'))
expected_seals={}
for ledger in (plan['input_seals'],bridge['input_seals'],bridge['outputs'],retirement['input_seals']):
    for name,expected in ledger.items():
        path=(Path(name) if Path(name).is_absolute() else ROOT/name).resolve()
        if path in expected_seals:assert expected_seals[path]==expected,name
        expected_seals[path]=expected
        paths.add(path)
for path,expected in expected_seals.items():assert seal(path)==expected,path
inputs={str(p):seal(p) for p in sorted(paths)}
current,code,ir=inventory(*[(STAGE/'query-payloads'/f'transport.info{k}').read_bytes() for k in (1,2,0)])
assert current['statistics']==queries[0]['statistics']
reference,oldcode,oldir=inventory(*[(OLD/'query-payloads'/f'transport.info{k}').read_bytes() for k in (1,2,0)])
assert oldir==(OLD/'generated-code/transport.llvmir').read_bytes()
assert reference['printed_static_instruction_count']==45524 and reference['printed_instruction_bytes']==237516 and reference['ELF_text_bytes']==237952
output=WORK/'generated-code';assert not output.exists();output.mkdir()
(output/'transport.text.bin').write_bytes(code);(output/'transport.llvmir').write_bytes(ir)
assert (output/'transport.llvmir').read_bytes()==ir and (output/'transport.text.bin').read_bytes()==code
assert {str(p):seal(p) for p in sorted(paths)}==inputs
assert git('rev-parse','HEAD')==BASELINE==git('rev-parse','@{u}') and not git('status','--porcelain')
result={'accepted':True,'baseline_revision':BASELINE,'input_seals':inputs,'module':seal(candidate),'candidate':current,'original':reference,
        'native_elapsed_seconds':owner['elapsed_seconds'],'native_peak_sampled_rss_bytes':owner['peak_sampled_rss_bytes'],
        'extracted_outputs':{str(p):seal(p) for p in output.iterdir()},
        'scope':'Saved one-pipeline/zero-dispatch static generated-code comparison with exact module/providers/feature-layout joins and whole printed ISA encoding/ELF prefix/embeddedIR bytes. Original Transport individual payloads complete; its historical aggregate query failed later on Endpoint and is not a whole-run pass. EmbeddedIR compiler-pass provenance, source-to-IR specialization, runtime words/cost/frequency/benefit and full scientific/release qualification remain unestablished. Lexical barrier absence does not remove source collectives.'}
target=WORK/'query-analysis.json';assert not target.exists();target.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({name:{k:r[k] for k in ('statistics','ELF_text_bytes','printed_static_instruction_count','printed_instruction_bytes','unprinted_ELF_tail_bytes','IR_bytes')} for name,r in [('candidate',current),('original',reference)]}))
