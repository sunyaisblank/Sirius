"""Read actual completed code captures; no dispatch or numerical qualification."""
from pathlib import Path
import hashlib,json,re,zipfile
R=Path.cwd();W=Path(__file__).resolve().parent

def doc(p):return json.loads(p.read_text())
def pin(p):
 b=p.read_bytes();return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
def body(folder,v):
 p=folder/v['ordinary_ir_export']['path']
 assert pin(p)=={k:v['ordinary_ir_export'][k] for k in ['bytes','sha256']}
 return p.read_text()

def intergroup(ir, v, event, folder):
    assert v['actual_entry']['entry_sysv_integer_arguments_2_3_4'] == [64, 1, 1]
    bindings = json.loads((folder / 'bindings.json').read_text())
    assert bindings['dispatches'][0]['groups_x'] == 1
    assert len(event['functions']) == 1 and event['functions'][0]['name'] == 'cs_variant'
    lines = ir.splitlines()
    definitions = dict(re.findall(r'^  (%[A-Za-z0-9_.]+) = (.*)$', ir, re.M))
    groups = [name for name, body in definitions.items() if re.search(r'add <8 x i32> .*<i32 0, i32 1, i32 2, i32 3, i32 4, i32 5, i32 6, i32 7>', body)]
    assert len(groups) == 1
    group, = groups
    def path(start, target):
        queue = [start]; parent = {start: None}
        for name in queue:
            if name == target:
                result = []
                while name is not None:
                    result.append(name + ' = ' + definitions.get(name, '<function argument>'))
                    name = parent[name]
                return result[::-1]
            for ref in re.findall(r'%[A-Za-z0-9_.]+', definitions.get(name, '')):
                if ref not in parent: parent[ref] = name; queue.append(ref)
        return []
    assert path(group, '%5'), 'Lane group IDs must derive from the original first-group-X argument'
    assert definitions['%owner_remaining'] == 'sub i32 %8, %5'
    assert definitions['%owner_tail'] == 'icmp ult i32 %owner_remaining, 8'
    assert definitions['%owner_partials'] == 'select i1 %owner_tail, i32 %owner_remaining, i32 0'
    assert 'splat (i32 -1)' in definitions['%mask.0']
    assert 'insertelement <8 x i32> %mask.1, i32 0, i32 %loop_counter2.0' in ir
    assert ir.index('%mask.0 = phi <8 x i32>') < ir.index('@llvm.masked.gather')
    calls = []
    for line in lines:
        if ' call ' not in line: continue
        m = re.search(r'@llvm.masked.(gather|scatter)\.[^(]+\((.*)\)', line)
        if not m: continue
        masks = re.findall(r'<8 x i1> (%[A-Za-z0-9_.]+)', m[2])
        assert len(masks) == 1, line
        mask_path = path(masks[0], '%mask.0')
        assert mask_path, ('Actual memory mask must retain pre-prologue tail-mask dependence', line)
        pointers = re.findall(r'<8 x ptr> (%[A-Za-z0-9_.]+)', m[2]); assert len(pointers) == 1
        address_path = path(pointers[0], group)
        calls.append({'operation': m[1], 'line': line.strip(), 'tail_mask_path': mask_path,
                      'group_address_path': address_path})
    gathers = [x for x in calls if x['operation'] == 'gather' and x['group_address_path']]
    scatters = [x for x in calls if x['operation'] == 'scatter' and x['group_address_path']]
    assert gathers and scatters
    return {'structural_premise_observed': True, 'original_block_arguments': [64, 1, 1],
        'guard_original_dispatch_groups_x': 1, 'lane_group_id_definition': group + ' = ' + definitions[group],
        'original_first_group_X_def_use': path(group, '%5'),
        'actual_row_address_and_mask_examples': [gathers[0], scatters[0]],
        'masked_memory_operations_with_tail_dependence': len(calls),
        'row_dependent_gathers': len(gathers), 'row_dependent_scatters': len(scatters),
        'review_limit': 'Mask def-use checks bind emitted lowering to the independently reviewed source mask composition. Guard exits before retained arithmetic; no active eight-row or numerical result inferred.'}

compiler=doc(W.parent/'compiler-result.json')
provider=doc(W.parent/'provider/provider-build-result.json');assert provider['passed']
driver=R/provider['library'];assert pin(driver)==provider['library_identity']
archive=R/'attestations/qualification/337d0bae/lavapipe-owner-bulk-zero-finite-outcome/diagnostics-and-trial.zip'
assert pin(archive)['sha256']=='7959c6c6be2671de0059657ffee4ccfae1b2c0dc19a2e257b4f6f05ee605236b'
with zipfile.ZipFile(archive) as z:
 previous=json.loads(z.read('out/lavapipe-owner-bulk-zero/ordinary/ordinary-code-review.json'))
rows={}
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
    ir=body(folder,v)
    assert ir.count(f'%reg = alloca [{count} x <8 x i32>], align 32')==1
    assert len(re.findall(r'call void @llvm\.memset\.p0\.i64\(ptr[^\n]*%reg, i8 0, i64 '+str(count*32)+r', i1 false\)',ir))==1
    assert ir.count('call void @llvm.memset.p0.i64(')==1
    assert 'store <8 x i32> zeroinitializer, ptr %reg' not in ir
    structural=intergroup(ir,v,event,folder)
    old=previous['rows'][kind]['total_called_function_bytes']
    size=sum(x['bytes'] for x in event['functions'])
    assert size<old
    rows[kind]={**structural,'complete_bulk_clear_observed':True,
        'original_Default_machine_pipeline_observed':True,
        'actual_complete_module':v['actual_complete_module'],
        'functions':event['functions'],'total_called_function_bytes':size,
        'prior_best_called_text_bytes':old,'reduction_percent':100*(old-size)/old,
        'actual_optimized_IR':{'path':str((folder/v['ordinary_ir_export']['path']).relative_to(R)),**pin(folder/v['ordinary_ir_export']['path'])},
        'actual_capture_verification':pin(folder/'capture-verification.json'),
        'arithmetical_equivalence_basis':'Independent exact source control-flow composition, unchanged lower-level arithmetic and emitted controls; no whole-module inverse or machine-equivalence proof claimed.'}
out={'structural_premise_observed':True,'rows':rows,
 'complete_bulk_clear_observed':True,'original_Default_machine_pipeline_observed':True,
 'actual_provider':{'path':str(driver.relative_to(R)),**pin(driver)},
 'prior_best_archive':{'path':str(archive.relative_to(R)),**pin(archive)},
 'scope':'Two original4B guard/canary captures only. Actual called code/source structure and complete clear/engine/allocator/features/FP restore are bound. Guard exits before retained arithmetic; originalC8/all32 and unchanged22+Decimal remain required. No runtime-benefit, active8row, frame, science or release claim.'}
target=W/'ordinary-code-review.json';assert not target.exists();target.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'passed':True,'called_text':{k:v['total_called_function_bytes'] for k,v in rows.items()},'reduction_percent':{k:v['reduction_percent'] for k,v in rows.items()}}))
