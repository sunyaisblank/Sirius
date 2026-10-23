"""Passive saved-code gate: exact full-zero substitution, no numerical trial."""
from pathlib import Path
import hashlib, json, re, zipfile
R = Path.cwd(); W = Path(__file__).resolve().parent
def doc(p): return json.loads(p.read_text())
def pin(p):
    with p.open('rb') as f: digest = hashlib.file_digest(f,'sha256').hexdigest()
    return {'bytes':p.stat().st_size,'sha256':digest}
def body(folder,v):
    p=folder/v['ordinary_ir_export']['path']
    assert pin(p)=={k:v['ordinary_ir_export'][k] for k in ['bytes','sha256']}
    return p.read_text()

def normalized_clear(ir,count,candidate):
    # Preserve all identifiers/instructions outside the one clear, and compare
    # attribute CONTENT rather than LLVM printer-assigned group numbers.
    extent=count*32
    assert ir.count(f'%reg = alloca [{count} x <8 x i32>], align 32')==1
    lines=ir.splitlines()
    if candidate:
        matches=[i for i,line in enumerate(lines) if re.fullmatch(
            rf'  (?:tail )?call void @llvm\.memset\.p0\.i64\(ptr ([^,]* )?%reg, i8 0, i64 {extent}, i1 false\)(?: #[0-9]+)?',line)]
        assert len(matches)==1, 'Require one exact full-size nonvolatile scratch memset'
        start,=matches
        assert re.search(r'\balign 32\b',lines[start])
        assert ir.count('call void @llvm.memset.p0.i64(')==1
        finish=start+1
    else:
        first='  store <8 x i32> zeroinitializer, ptr %reg, align 32'
        start,=[i for i,line in enumerate(lines) if line==first]
        expected=[first]
        for i in range(1,count):
            expected.extend([f'  %.fca.{i}.gep = getelementptr inbounds nuw i8, ptr %reg, i64 {32*i}',
                             f'  store <8 x i32> zeroinitializer, ptr %.fca.{i}.gep, align 32'])
        finish=start+len(expected); assert lines[start:finish]==expected
    lines[start:finish]=[f'  ; EQUIVALENT_FULL_ZERO %reg bytes={extent} align=32']
    # Keep the sentinel in the instruction stream; remove only printer comments.
    lines=[line for line in lines if not line.lstrip().startswith(';') or 'EQUIVALENT_FULL_ZERO' in line]
    if candidate:
        declarations=[i for i,line in enumerate(lines) if line.startswith('declare void @llvm.memset.p0.i64(')]
        assert len(declarations)==1
        del lines[declarations[0]]
    text='\n'.join(lines)
    groups=dict(re.findall(r'^attributes #(\d+) = (.*)$',text,re.M))
    text=re.sub(r'^attributes #\d+ = .*\n?', '', text, flags=re.M)
    text=re.sub(r'#(\d+)\b',lambda m:'#ATTR'+groups[m[1]],text)
    assert '@llvm.memset.' not in text
    # Empty lines and explanatory comments have no executable meaning; no SSA
    # renumbering, operation reordering, arithmetic/flag rewrite is permitted.
    return '\n'.join(line for line in text.splitlines() if line.strip())

compiler=doc(R/'out/native-cooperative-transport-intergroup/compiler-result.json')
provider=doc(W.parent/'provider-build-result.json'); assert provider['passed']
driver=R/provider['library']; assert pin(driver)==provider['library_identity']
archive=R/'attestations/qualification/33623035/lavapipe-owner-machine-codegen-finite-outcome'
manifest=doc(archive/'MANIFEST.json'); zpath=archive/'diagnostics-and-trial.zip'
assert pin(zpath)==manifest['archives'][zpath.name]
member='out/lavapipe-owner-machine-codegen/ordinary/ordinary-code-review.json'
with zipfile.ZipFile(zpath) as z: raw=z.read(member)
assert len(raw)==manifest['diagnostics'][member]['bytes'] and hashlib.sha256(raw).hexdigest()==manifest['diagnostics'][member]['sha256']
prior_transport=json.loads(raw)['rows']['Transport']
assert prior_transport['optimized_IR_unchanged'] and prior_transport['greedy_allocator_factories_in_selected_engine_emission']>0
assert prior_transport['fast_allocator_factories_in_selected_engine_emission']==0
rows={'Transport':{'structural_premise_observed':True,
    'complete_module':prior_transport['actual_complete_module'],
    'ordinary_Default_source_proof_reused':True,
    'actual_current_provider_Transport_capture_performed':False,
    'proof_archive':str(zpath.relative_to(R)), 'proof_archive_identity':pin(zpath),
    'proof_member':member, 'prior_actual_observation':prior_transport,
    'scope':'Historical cooperative T actual greedy/IR proof plus unchanged ordinary source/type-exclusion composition. Current E/DP provider captures do not retroactively make this a current-provider T code observation.'}}
for label,kind,count in [('endpoint','Endpoint',3060),('dopri_phase','DopriPhase',3145)]:
    folder=W/label; o=doc(folder/'owner.json')
    assert o['passed'] and o['source_unchanged'] and o['whole_inputs_unchanged']
    assert o['stop_reason'] is None and not o['cleanup_errors'] and o['remaining_owned_processes']==[]
    assert o['observed_births'] and all(x['absent'] for x in o['observed_births'])
    assert doc(folder/'inputs-before.json')==doc(folder/'inputs-after.json')
    v=doc(folder/'capture-verification.json'); b=doc(folder/'bindings.json'); p=doc(folder/'machine-policy.json')
    assert v['passed'] and v['actual_provider_sha256']==provider['library_identity']['sha256']
    assert v['actual_complete_module']['sha256']==compiler['modules'][kind]['sha256']
    entry=v['actual_entry']; assert entry['dispatch_id']==0 and entry['actual_mxcsr'] & 0xe040==0
    assert entry['entry_sysv_integer_arguments_2_3_4']==[64,1,1] and v['callback_restoration_observed']
    assert len(b['dispatches'])==1 and b['dispatches'][0]['completed'] and b['dispatches'][0]['groups_x']==1
    event,=[x for x in b['jit_events'] if x.get('file')==entry['file']]
    assert str(driver) in event['maps_at_first_call']
    assert pin(folder/event['file'])['sha256']==event['sha256']
    for section in event['live_sections_after_pointer_lookup']:
        assert pin(folder/section['path'])=={k:section[k] for k in ['bytes','sha256']}
    generation=v['ordinary_ir_export']['module_generation']
    lookups=[x for x in p['lookups'] if x['completed'] and x['pointer']==entry['actual_entry_pc'] and x['module_generation']==generation]
    assert lookups
    engine,={x['engine'] for x in lookups}
    creation,=[x for x in p['engines'] if x['completed'] and x['engine']==engine and x['module_generation']==generation]
    assert creation['initial_OptLevel']==2 and not p['bridges']
    emissions=[x for x in p['emissions'] if x['completed'] and x['engine']==engine and x['module_generation']==generation]
    assert emissions and all(any(l['sequence']<e['sequence']<e['return_sequence']<l['return_sequence'] and l['thread']==e['thread'] for l in lookups) for e in emissions)
    factories=[x for x in p['factories'] if x['emission_id'] in {e['id'] for e in emissions}]
    assert factories and all(x['allocator']=='greedy' and x['engine']==engine and x['module_generation']==generation for x in factories)
    for x in factories:
        e=p['emissions'][x['emission_id']]
        assert e['sequence']<x['sequence']<e['return_sequence'] and e['thread']==x['thread']
    clear,=[x for x in p['clears'] if x['module_generation']==generation and x['jit_file']==entry['file']]
    assert clear['completed'] and clear['complete_zero_bytes_observed'] and clear['bytes']==count*32 and clear['dispatch_id']==0
    assert clear['resulting_sha256']==hashlib.sha256(b'\0'*(count*32)).hexdigest()
    binding,=[x for x in p['clear_bindings'] if x['jit_file']==entry['file'] and x['address']==clear['address'] and x['symbol']==clear['symbol']]
    assert binding['relocation_kind']==1 and binding['addend']==0
    libc=doc(folder/'codegen-symbols.json')['libc']; assert pin(Path(libc['path']))['sha256']==libc['sha256']
    assert Path(clear['library']).resolve()==Path(libc['path']).resolve()
    prior=R/'out/lavapipe-intergroup-row-owner/ordinary'/label; pv=doc(prior/'capture-verification.json')
    current_norm=normalized_clear(body(folder,v),count,True)
    prior_norm=normalized_clear(body(prior,pv),count,False)
    assert current_norm==prior_norm, ('IR changed outside declared full-clear equivalence',kind)
    rows[kind]={'structural_premise_observed':True,
        'complete_bulk_clear_observed':True, 'remaining_optimized_IR_equivalent':True,
        'original_Default_machine_pipeline_observed':True, 'complete_module':v['actual_complete_module'],
        'actual_engine':engine, 'actual_initial_OptLevel':creation['initial_OptLevel'],
        'actual_greedy_factories':factories, 'actual_complete_clear':clear,
        'resolved_libc_binding':binding, 'libc':libc,
        'remaining_IR_normalized_sha256':hashlib.sha256(current_norm.encode()).hexdigest(),
        'IR_equivalence_rule':'One contiguous exact 32-byte-stride full-zero GEP/store span -> one exact complete nonvolatile scratch memset at same position; comments/empty lines and attribute-group numeric labels only. Remaining identifiers, operations, constants, flags, types and control flow exact.',
        'functions':event['functions'], 'total_called_function_bytes':sum(x['bytes'] for x in event['functions']),
        'capture_verification_identity':pin(folder/'capture-verification.json')}
out={'structural_premise_observed':True,'rows':rows,
    'complete_bulk_clear_observed':True,'remaining_optimized_IR_equivalent':True,
    'original_Default_machine_pipeline_observed':True,
    'actual_provider':{'path':str(driver.relative_to(R)),**pin(driver)},
    'scope':'Two 4B guard/canary E/DP captures. All scratch zeros retained and observed; no retained arithmetic acceptance, benefit, active-eight-row/frame/scientific/release claim. T proof explicitly reused, no new T capture.'}
target=W/'ordinary-code-review.json'; assert not target.exists()
target.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'passed':True,'E_DP_clear_bytes':{'Endpoint':97920,'DopriPhase':100640}}))
