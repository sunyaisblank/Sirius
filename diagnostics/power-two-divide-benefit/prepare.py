"""Restore the accepted producer pair for one unchanged original benefit trial."""
from pathlib import Path
import hashlib, json, shutil, subprocess

R = Path.cwd().resolve()
Q = Path(__file__).resolve().parent
W = Q / 'benefit'
T = R / 'bin/linux-gcc/power-two-divide-benefit'
REVISION = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip() == REVISION
assert not subprocess.check_output(['git','status','--porcelain'],text=True)
assert not T.exists() and not W.exists()
T.mkdir()
W.mkdir()
def identity(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for data in iter(lambda:stream.read(1024*1024),b''): h.update(data)
    return {'bytes':path.stat().st_size,'sha256':h.hexdigest()}
def sha(path): return identity(path)['sha256']
def write(path,value):path.write_text(json.dumps(value,indent=2)+'\n')
archive = R / 'attestations/software-vulkan/4bcd987/power-two-divide-prerequisites'
manifest_path = archive / 'whole-byte-manifest.json'
manifest = json.loads(manifest_path.read_text())
assert sha(manifest_path) == '36502f200a07fc4da0bba1df5baac0adc2fb8106342baa71a9943ce74ff8d5ef'
old_work = 'out/power-two-divide-controls/'
old_stage = 'bin/linux-gcc/power-two-divide-controls/'
def preserved(relative):
    path = archive / 'files' / relative
    assert identity(path) == manifest['files'][relative], relative
    return path
disposition = json.loads(preserved(old_work+'disposition.json').read_text())
assert disposition['decision'] == 'PROCEED_TO_ORIGINAL_BENEFIT_GATES'
assert disposition['original_eight_controls_passed'] and not disposition['production_adopted']
producer = json.loads(preserved(old_work+'producer-acceptance.json').read_text())
assert producer['passed'] and producer['unchanged_array_count'] == 35 and len(producer['changed_arrays']) == 5
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
for relative,pin in manifest['files'].items():
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
(W/'run_limb_benefit_slot.py').write_text(wrapper.replace(old_allowlist,'bin/linux-gcc/power-two-divide-benefit/original_'))
for mode,name in [('baseline','baseline.json'),('candidate','build-and-arrays.json')]:
    write(W/name,{'arrays':producer['actual_ELF_arrays'][mode]})
tracked = set(subprocess.check_output(['git','ls-files','-z']).decode().split('\0'))
prior_inputs = json.loads(preserved(old_work+'control-8/inputs-before.json').read_text())
current_joins = {}
shared = []
retired_roots = [(R/name).resolve() for name in manifest['roots']]
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
    'out/power-two-divide-screen/candidate/src/sirius/kernels/portable_binary32.h',
    'out/power-two-divide-screen/source-freeze.json',
    'out/power-two-divide-screen/module-bindings.json',
    'out/power-two-divide-screen/compiler-result.json',
    'out/power-two-divide-screen/compiler-screen-disposition.json',
    'out/power-two-divide-screen/independent-contract-review.json',
    'out/power-two-divide-arithmetic/arithmetic-result.json',
    old_work+'preparation.json', old_work+'runtime-preparation.json',
    old_stage+'headers/retained_kernels.h', old_stage+'compute-sealed.d',
    old_stage+'observer.d', old_stage+'candidate.map',
]
for label in ['compile-compute-sealed','compile-observer','link-candidate']:
    chain_names += [old_work+label+'/'+name for name in ['owner.json','inputs-before.json','inputs-after.json']]
for label in range(1,9):
    chain_names += [old_work+'control-'+str(label)+'/'+name for name in ['owner.json','gtest.xml']]
for relative in manifest['files']:
    if relative.startswith('bin/linux-gcc/power-two-divide-screen/modules/') and relative.endswith('.spv'):
        chain_names.append(relative)
candidate_chain = {relative: identity(preserved(relative)) for relative in chain_names}
shared += [str(preserved(relative)) for relative in chain_names]
assert json.loads((Q/'control-acceptance.json').read_text())['passed']
assert json.loads((Q/'endpoint-route-acceptance.json').read_text())['passed']
assert json.loads(preserved('out/power-two-divide-arithmetic/arithmetic-result.json').read_text())['passed']
source_freeze = json.loads(preserved('out/power-two-divide-screen/source-freeze.json').read_text())
assert candidate_chain['out/power-two-divide-screen/candidate/src/sirius/kernels/portable_binary32.h'] == source_freeze['candidate/src/sirius/kernels/portable_binary32.h']
header_relative = old_stage+'headers/retained_kernels.h'
header_pin = candidate_chain[header_relative]
for label in ['compile-compute-sealed','compile-observer']:
    assert producer['compile_provenance'][label]['actual_dependencies'][header_relative] == header_pin
    assert identity(preserved(old_work+label+'/owner.json')) == producer['compile_provenance'][label]['owner']
assert identity(preserved(old_work+'link-candidate/owner.json')) == producer['link_owner']
shared += resources + [str(Q/name) for name in ['prepare.py','prove_reuse.py','producer-acceptance.json','control-acceptance.json',
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
              scope='One original fixed A1/B1/B2/A2 exact-power-of-two division trial. Previously accepted matched40-array producers and original8controls. Exact unchanged reader/comparator/launcher/gates; only wrapper executable allowlist adapted. No earlier rejected candidate observations reused as samples.')
assert freeze['sequence'] == ['a1','b1','b2','a2'] and freeze['within_slot_repeats'] == 3
assert freeze['complete_call_minimum_improvement_percent'] == 5
assert freeze['every_inclusive_stage_minimum_improvement_percent'] == -5
assert freeze['observed_instance_api_minimum_improvement_percent'] == 0 and freeze['all20_gates_required']
assert W.parents[2] == R
write(W/'retention-gate-frozen.json',freeze)
write(Q/'preparation.json',{'source_revision':REVISION,'prior_manifest_sha256':sha(manifest_path),
    'producer_prerequisites':{'actual40array_ELF_proof':True,'same_observer_CPP':True,'same_current_stock_adapter_ABI':True,
                           'original8controls_passed':True,'both_fresh_compiles_used_sealed_header':True},
    'executables':executables,'current_byte_exact_dependency_joins':current_joins,
    'candidate_source_and_compile_chain':candidate_chain,'unchanged_original_apparatus':pins,'resources':len(resources),'sequence':freeze['sequence'],
    'scope':'Exact producer/input reuse without recompilation, prior observations not performance samples. Frozendeclared inputs do not claim every implicit operating-system dependency.'})
print(json.dumps({'prepared':True,'sequence':freeze['sequence'],'reused_dependency_joins':len(current_joins),'resources':len(resources),'executables':executables}))
