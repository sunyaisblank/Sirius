"""Restore the accepted producer pair for one unchanged original benefit trial."""
from pathlib import Path
import ast, hashlib, json, re, shutil, struct, subprocess
import xml.etree.ElementTree as ET

R = Path.cwd().resolve()
Q = Path(__file__).resolve().parent
W = Q / 'benefit'
T = R / 'bin/linux-gcc/computed-carry-exit-benefit'
REVISION = 'a1d6e75299a2bf250e5ce7d44ad5965f8acb4299'
MANIFEST_SHA256 = '49a9f631abd228b60bd487a78e543fd2c121b3b55d4fe5f8d15c64c71742f7c2'
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip() == REVISION
assert not subprocess.check_output(['git','status','--porcelain'],text=True)
assert not T.exists() and not W.exists()
def identity(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for data in iter(lambda:stream.read(1024*1024),b''): h.update(data)
    return {'bytes':path.stat().st_size,'sha256':h.hexdigest()}
def sha(path): return identity(path)['sha256']
def write(path,value):path.write_text(json.dumps(value,indent=2)+'\n')
archive = R / 'attestations/software-vulkan/4bcd987/computed-carry-exit-original-prerequisites'
manifest_path = archive / 'manifest.json'
manifest = json.loads(manifest_path.read_text())
assert sha(manifest_path) == MANIFEST_SHA256
old_work = 'out/computed-carry-exit-original-prerequisites/'
old_stage = 'bin/linux-gcc/computed-carry-exit-original-prerequisites/controls/'
def preserved(relative):
    path = archive / 'files' / relative
    assert identity(path) == manifest['files']['files/' + relative], relative
    return path
disposition = json.loads(preserved(old_work+'disposition.json').read_text())
assert disposition['decision'] == 'PROCEED_TO_ORIGINAL_BENEFIT_GATES'
assert disposition['original_eight_controls_passed'] and not disposition['production_adopted']
producer = json.loads(preserved(old_work+'producer-acceptance.json').read_text())
assert producer['passed'] and producer['unchanged_array_count'] == 27 and len(producer['changed_arrays']) == 13
assert set(producer['actual_ELF_owner_arrays']) == {'baseline', 'candidate'}
# Read only the sealed prerequisites before restoring any executable or trial root.
source = json.loads(preserved(old_work+'operator-context.json').read_text())
def prerequisite_owner(label, provider=False):
    owner = json.loads(preserved(old_work+label+'/owner.json').read_text())
    assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0, label
    assert owner['source'] == {'revision':source['revision'], 'status':''}, label
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged'], label
    assert owner['stop_reason'] is None and owner['cleanup_errors'] == [] and owner['remaining_owned_processes'] == [], label
    assert owner['observed_births'] and all(b['absent'] for b in owner['observed_births']), label
    if provider:
        assert owner['provider_mappings'], label
    assert preserved(old_work+label+'/inputs-before.json').read_bytes() == preserved(old_work+label+'/inputs-after.json').read_bytes(), label
    return owner

controls = json.loads(preserved(old_work+'control-acceptance.json').read_text())
routes = json.loads(preserved(old_work+'endpoint-route-acceptance.json').read_text())
assert controls['passed'] and routes['passed'] and len(controls['original_controls']) == 8
schema = controls['timestamp_schema']
assert (schema['complete_readbacks'],schema['ordered_markers'],schema['actual_module_command_joins']) == (32,51,54)
assert schema['all_readbacks_exact_to_repaired_original_baseline']
assert routes['complete_readbacks'] == 56 and routes['complete_readbacks_exact_to_prior_original_control']
assert len(routes['routes']) == 4 and {(v['declined_RTE32'],v['wide']) for v in routes['routes']} == {(d,w) for d in [False,True] for w in [False,True]}
original_controls = json.loads((R/'attestations/software-vulkan/41fe5f2/portable-denominator-preparation/files/out/portable-denominator-preparation/original-controls/actions.json').read_text())
archived_actions = json.loads(preserved(old_work+'actions.json').read_text())
for index in range(1,9):
    label = 'control-'+str(index)
    owner = prerequisite_owner(label, True)
    assert owner['argv'] == archived_actions[label]['argv']
    assert owner['argv'][1] == original_controls[label]['argv'][1]
    assert owner['limits'] == original_controls[label]['limits']
    row = controls['original_controls'][index-1]
    assert row['label'] == label and row['owner'] == identity(preserved(old_work+label+'/owner.json'))
    assert row['xml'] == identity(preserved(old_work+label+'/gtest.xml'))
    xml = ET.parse(preserved(old_work+label+'/gtest.xml')).getroot()
    assert int(xml.get('tests','-1')) == 1
    assert not any(xml.findall('.//'+tag) for tag in ['failure','error','skipped'])
    for node in [xml,*xml.findall('.//testsuite')]:
        assert all(int(node.get(key,'0')) == 0 for key in ['failures','errors','skipped','disabled'])
    cases = xml.findall('.//testcase')
    assert len(cases) == 1 and cases[0].get('status') == 'run' and cases[0].get('result') == 'completed'
    assert owner['argv'][1] == '--gtest_filter='+cases[0].get('classname')+'.'+cases[0].get('name')
protected = json.loads(preserved(old_work+'protected-before.json').read_text())
assert len(protected) == controls['protected_payloads_unchanged'] == 114
for name,pin in protected.items():
    assert identity(R/name) == pin, name
products = json.loads(preserved(old_work+'products-acceptance.json').read_text())
assert products['passed'] and len(products['modules']) == 35 and sum(row['changed'] for row in products['modules']) == 15
for label in ['compile-products-stock','compile-products-candidate','verify-products','native-policy-selfcheck',
              'compile-compute-sealed','compile-observer','link-candidate','arithmetic/validate-arithmetic']:
    prerequisite_owner(label)
assert json.loads(preserved(old_work+'arithmetic/arithmetic-result.json').read_text())['passed']
for mode in ['stock','candidate']:
    label = 'arithmetic/execute-gpu-'+mode
    prerequisite_owner(label, True)
    assert preserved(old_work+label+'/actual.bin').read_bytes() == preserved(old_work+'arithmetic/expected.bin').read_bytes()
# Reuse the sealed original ELF reader's definitions, without its executable footer.
reader_tree = ast.parse(preserved(old_work+'verify_producer.py').read_text())
reader_defs = ast.Module(body=[n for n in reader_tree.body if isinstance(n,(ast.Import,ast.ImportFrom,ast.FunctionDef))],type_ignores=[])
reader = {}
exec(compile(reader_defs,'sealed_original_ELF_reader','exec'),reader)
for mode in ['baseline','candidate']:
    directory = 'stock' if mode == 'baseline' else 'candidate'
    generated = 'bin/linux-gcc/computed-carry-exit-original-prerequisites/'+directory+'/'
    executable = preserved(old_stage+'original_'+mode+'_controls')
    assert identity(executable) == producer[mode+'_executable']
    actual,_ = reader['elf_arrays'](executable,reader['arrays'](preserved(generated+'retained_kernels.h')))
    assert actual == producer['actual_ELF_arrays'][mode]
    owner_arrays = {}
    for count,name,body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',preserved(generated+'retained_owner_kernels.h').read_text(),re.S):
        words = [int(v) for v in re.findall(r'(\d+)u',body)]
        assert len(words) == int(count)
        owner_arrays[name] = struct.pack('<'+str(len(words))+'I',*words)
    assert set(owner_arrays) == {'kTransportShader','kEndpointShader'}
    actual,_ = reader['elf_arrays'](executable,owner_arrays,'retained_owner')
    assert actual == producer['actual_ELF_owner_arrays'][mode]
for label in ['compile-compute-sealed','compile-observer']:
    proof = producer['compile_provenance'][label]
    assert proof['owner'] == identity(preserved(old_work+label+'/owner.json'))
    assert proof['actual_dependencies'][old_stage+'headers/retained_kernels.h'] == identity(preserved(old_stage+'headers/retained_kernels.h'))
for name in ['retained_owner_kernels.h','retained_owner_layout.h']:
    assert producer['compile_provenance']['compile-compute-sealed']['actual_dependencies'][old_stage+'headers/'+name] == identity(preserved(old_stage+'headers/'+name))
assert producer['link_owner'] == identity(preserved(old_work+'link-candidate/owner.json'))
assert producer['link_map'] == identity(preserved(old_stage+'candidate.map')) and producer['old_archive_compute_member_excluded']
T.mkdir()
W.mkdir()
for name in ['producer-acceptance.json','control-acceptance.json','endpoint-route-acceptance.json',
             'audited-observer.cpp','protected-before.json']:
    shutil.copyfile(preserved(old_work+name),Q/name)
assert sha(Q/'audited-observer.cpp') == '399cd52cc6d9fc89e4b97ff63b06a54b06b4f80128f5a99c2c6d1278c62f3dc3'
for name,pin in json.loads((Q/'protected-before.json').read_text()).items():
    assert identity(R/name) == pin,name
executables = {}
for mode in ['baseline','candidate']:
    source = preserved(old_stage+'original_'+mode+'_controls')
    target = T/source.name
    shutil.copyfile(source,target)
    target.chmod(0o755)
    assert identity(target) == producer[mode+'_executable']
    executables[mode] = {'path':str(target.relative_to(R)),'source_revision':REVISION,
                         **identity(target),'observer_source_sha256':sha(Q/'audited-observer.cpp')}
resources = []
for archived,pin in manifest['files'].items():
    relative = archived.removeprefix('files/')
    if relative.startswith(old_stage+'resources/'):
        source = preserved(relative)
        target = T/'resources'/relative.removeprefix(old_stage+'resources/')
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(source,target)
        assert identity(target) == pin
        resources.append(str(target))
reference = R/'attestations/software-vulkan/8a27321/vulkan-descriptor-reuse/files/out/vulkan-descriptor-reuse/audited-benefit'
pins = {'software_contracts.py':'d6d47244f9a93c22d30e1a459f25d0207c9e64682a8814bead93eb455469cdb9',
        'compare_software.py':'2aa5e566bc1cb4a2767a73f1676b94cc894682f28d6b107ecec28548c026bac3',
        'run_sequence.py':'f8e754bcb137a43e610e3a794cb9065691971c227c557efb9c6e5a057c595342',
        'expected-readback-layout.json':'a40ce35b744ff05c987833de8337ec6f5445daea9553348cdb9209ad78314392'}
for name,pin in pins.items():
    assert sha(reference/name) == pin
    shutil.copyfile(reference/name,W/name)
for name in ['run_owned.py','owned_guard.py']:
    shutil.copyfile(preserved(old_work+name),W/name)
old_allowlist = "bin/linux-gcc/vulkan-descriptor-reuse/audited_"
wrapper = (reference/'run_limb_benefit_slot.py').read_text()
assert wrapper.count(old_allowlist) == 1
(W/'run_limb_benefit_slot.py').write_text(wrapper.replace(old_allowlist,'bin/linux-gcc/computed-carry-exit-benefit/original_'))
for mode,name in [('baseline','baseline.json'),('candidate','build-and-arrays.json')]:
    write(W/name,{'arrays':producer['actual_ELF_arrays'][mode]})
tracked = set(subprocess.check_output(['git','ls-files','-z']).decode().split('\0'))
prior_inputs = json.loads(preserved(old_work+'control-8/inputs-before.json').read_text())
current_joins = {}
shared = []
retired_roots = [(R/name).resolve() for name in ['out/computed-carry-exit-original-prerequisites', 'bin/linux-gcc/computed-carry-exit-original-prerequisites', 'out/computed-carry-exit-prerequisites', 'bin/linux-gcc/computed-carry-exit-prerequisites']]
for name,pin in prior_inputs.items():
    path = R/name
    if (name in tracked or name.startswith(('out/','attestations/'))
            or any(path.resolve().is_relative_to(owned) for owned in retired_roots)):
        continue
    assert identity(path) == {key:pin[key] for key in ['bytes','sha256']},name
    current_joins[name] = pin
    shared.append(str(path))
# The private candidate source is outside Git's tracked-stock source seal.
# Bind its preserved source, modules and both actual compile dependencies explicitly.
chain_names = [
    old_work+'candidate/src/sirius/kernels/retained_triple.slang',
    old_work+'candidate.patch', old_work+'source-freeze.json', old_work+'operator-context.json',
    old_work+'recovery.json', old_work+'prepared-phase.json', old_work+'products-acceptance.json',
    old_work+'arithmetic/arithmetic-result.json', old_work+'arithmetic/reference-validation.json',
    old_work+'arithmetic/fixture.json', old_work+'arithmetic/input.bin', old_work+'arithmetic/expected.bin',
    old_work+'preparation.json', old_work+'runtime-preparation.json',
    old_stage+'headers/retained_kernels.h', old_stage+'headers/retained_owner_kernels.h',
    old_stage+'headers/retained_owner_layout.h', old_stage+'compute-sealed.d',
    old_stage+'observer.d', old_stage+'candidate.map',
]
for label in ['compile-products-stock','compile-products-candidate','verify-products',
              'compile-compute-sealed','compile-observer','link-candidate', 'arithmetic/validate-arithmetic']:
    chain_names += [old_work+label+'/'+name for name in ['owner.json','inputs-before.json','inputs-after.json']]
for label in range(1,9):
    chain_names += [old_work+'control-'+str(label)+'/'+name for name in ['owner.json','gtest.xml']]
for archived in manifest['files']:
    relative = archived.removeprefix('files/')
    if relative.startswith(('bin/linux-gcc/computed-carry-exit-original-prerequisites/stock/',
                            'bin/linux-gcc/computed-carry-exit-original-prerequisites/candidate/')) and relative.endswith('.spv'):
        chain_names.append(relative)
candidate_chain = {relative: identity(preserved(relative)) for relative in chain_names}
shared += [str(preserved(relative)) for relative in chain_names]
assert json.loads((Q/'control-acceptance.json').read_text())['passed']
assert json.loads((Q/'endpoint-route-acceptance.json').read_text())['passed']
assert json.loads(preserved(old_work+'arithmetic/arithmetic-result.json').read_text())['passed']
source_freeze = json.loads(preserved(old_work+'source-freeze.json').read_text())
candidate_relative = old_work+'candidate/src/sirius/kernels/retained_triple.slang'
assert candidate_chain[candidate_relative] == source_freeze[str(R/candidate_relative)]
assert candidate_chain[candidate_relative] == {'bytes':15882,'sha256':'c9df98c1e05c129cd31b03c7ed354e7bac01d2568f4a9fc518e12433c95fda12'}
header_relative = old_stage+'headers/retained_kernels.h'
header_pin = candidate_chain[header_relative]
for label in ['compile-compute-sealed','compile-observer']:
    assert producer['compile_provenance'][label]['actual_dependencies'][header_relative] == header_pin
    assert identity(preserved(old_work+label+'/owner.json')) == producer['compile_provenance'][label]['owner']
assert identity(preserved(old_work+'link-candidate/owner.json')) == producer['link_owner']
for name in ['retained_owner_kernels.h','retained_owner_layout.h']:
    relative = old_stage+'headers/'+name
    assert producer['compile_provenance']['compile-compute-sealed']['actual_dependencies'][relative] == candidate_chain[relative]
# The already-preserved focused pass remains an explicit archive input.
focused = R/'attestations/software-vulkan/4bcd987/computed-carry-exit-focused'
assert sha(focused/'manifest.json') == 'ccadd9af9a99f37b997d43b7c7b7490d6e9cda851589977ffd1bc1364f16f7c9'
shared += [str(focused/'manifest.json'), str(focused/'files/out/computed-carry-exit-prerequisites/focused-result.json')]

shared += resources + [str(Q/name) for name in ['prepare.py','prove_reuse.py','finish_trial.py','producer-acceptance.json','control-acceptance.json',
                    'endpoint-route-acceptance.json','audited-observer.cpp','protected-before.json','preparation.json']]
shared += [str(manifest_path),str(Q/'producer-reuse.json')]
shared += [str(W/name) for name in ['software_contracts.py','compare_software.py','run_sequence.py','run_limb_benefit_slot.py',
                    'baseline.json','build-and-arrays.json','expected-readback-layout.json','retention-gate-frozen.json']]
assert not any(Path(name).resolve().is_relative_to(owned)
               for name in shared for owned in retired_roots)
old_actions = json.loads((reference/'actions.json').read_text())
actions = {}
for label in ['a1','b1','b2','a2']:
    mode = 'baseline' if label.startswith('a') else 'candidate'
    actions[label] = {'revision':REVISION,'kind':'test','expected_exit':0,
        'argv':['python3','-B',str(W/'run_limb_benefit_slot.py'),executables[mode]['path'],old_actions[label]['argv'][-1]],
        'limits':old_actions[label]['limits'],'software_provider':True,
        'frozen_inputs':list(dict.fromkeys(shared+[str(R/executables[mode]['path'])])),
        'scope':'One fixed original default A/B/B/A slot: original three repeats/family, complete numerical assertions,32readbacks51markers54actualcommand joins and20+2 gates. No pooling, retry or full scientific/rendering/release verdict.'}
write(W/'actions.json',actions)
freeze = json.loads((reference/'retention-gate-frozen.json').read_text())
freeze.update(candidate_revision=REVISION,baseline_revision=REVISION,executables=executables,
              actions_sha256=sha(W/'actions.json'),
              scope='One original fixed A1/B1/B2/A2 computed-carry exit trial. Previously accepted matched40+2-array producers and original8controls. Exact unchanged reader/comparator/launcher/gates; only wrapper executable allowlist adapted. No earlier rejected candidate observations reused as samples.')
assert freeze['sequence'] == ['a1','b1','b2','a2'] and freeze['within_slot_repeats'] == 3
assert freeze['complete_call_minimum_improvement_percent'] == 5
assert freeze['every_inclusive_stage_minimum_improvement_percent'] == -5
assert freeze['observed_instance_api_minimum_improvement_percent'] == 0 and freeze['all20_gates_required']
assert W.parents[2] == R
write(W/'retention-gate-frozen.json',freeze)
write(Q/'preparation.json',{'source_revision':REVISION,'prior_manifest_sha256':sha(manifest_path),
    'producer_prerequisites':{'actual40array_ELF_proof':True,'actual2owner_array_ELF_proof':True,'same_observer_CPP':True,'same_current_stock_adapter_ABI':True,
                           'original8controls_passed':True,'both_fresh_compiles_used_sealed_header':True},
    'executables':executables,'current_byte_exact_dependency_joins':current_joins,
    'candidate_source_and_compile_chain':candidate_chain,'unchanged_original_apparatus':pins,'resources':len(resources),'sequence':freeze['sequence'],
    'scope':'Exact producer/input reuse without recompilation, prior observations not performance samples. Frozendeclared inputs do not claim every implicit operating-system dependency.'})
print(json.dumps({'prepared':True,'sequence':freeze['sequence'],'reused_dependency_joins':len(current_joins),'resources':len(resources),'executables':executables}))
