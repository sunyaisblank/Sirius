from pathlib import Path
import hashlib, json, re, shutil, struct, subprocess

ROOT = Path.cwd()
WORK = ROOT / 'out/portable-stage-seven-owner-construction'
PHASE = WORK / 'original-owner-benefit'
PRIOR = ROOT / 'attestations/software-vulkan/dc714b9/guarded-upmultiply-trial/diagnostic'
PHASE.mkdir()
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
edits = {}

def replace(text, old, new, ledger):
    assert text.count(old) == 1, (old, text.count(old))
    ledger.append({'old':old,'new':new})
    return text.replace(old,new)

reader = (PRIOR/'software_contracts.py').read_text()
reader_edits = []
reader = replace(reader, 'ROOT = Path(__file__).resolve().parent.parent.parent',
                 'ROOT = Path(__file__).resolve().parents[3]', reader_edits)
reader = replace(reader, "    result['expected_complete_readback_count'] = 32\n    return result\n", '''    result['expected_complete_readback_count'] = 32
    optional = document(WORK / 'owner-module-authority.json')
    assert {k for k in props if k.startswith('optional_module_')} == ({'optional_module_7', 'optional_module_8'} if producer == 'candidate' else set())
    if producer == 'candidate':
        for index, role in [(7, 'transport'), (8, 'endpoint')]:
            assert fields(props[f'optional_module_{index}']) == {'stage':role, 'bytes':str(optional[role]['bytes']), 'sha256':optional[role]['sha256']}
    expected_executed = set()
    executed = []
    handles = {}
    for marker_index, (_, stage_text) in enumerate(order):
        roles = ['endpoint', 'dense'] if stage_text == 'endpoint+dense' else [stage_text.split(':rows=')[0]]
        for command_index, role in enumerate(roles):
            key = f'executed_module_{marker_index}_{command_index}'
            expected_executed.add(key)
            module = fields(props[key])
            assert set(module) == {'stage','bytes','sha256','loaded_index','handle','bindings'}
            index = stages.index(role)
            if producer == 'candidate' and role == 'transport': index = 7
            if producer == 'candidate' and role == 'endpoint' and (':rows=' not in stage_text or stage_text.endswith(':rows=24')): index = 8
            authority = optional[role] if index >= 7 else {'bytes':selected[role]['bytes'], 'sha256':selected[role]['sha256']}
            assert module['stage'] == role and module['bytes'] == str(authority['bytes']) and module['sha256'] == authority['sha256']
            assert module['loaded_index'] == str(index)
            assert module['bindings'] == ('3' if index >= 7 else '2')
            handle = int(module['handle']); assert 0 <= handle < 2**32
            assert index not in handles or handles[index] == handle
            handles[index] = handle
            executed.append(module)
    assert {k for k in props if k.startswith('executed_module_')} == expected_executed
    assert len(executed) == 54 and len(set(handles.values())) == len(handles)
    result['executed_modules'] = executed
    create = number(props['observed_create_ms'])
    provisioning = number(props['owner_provisioning_ms'])
    assert create > 0 and provisioning <= create
    assert (provisioning > 0) == (producer == 'candidate')
    disabled_attempt = number(props['disabled_attempt_ms'])
    disabled_sample = number(props['disabled_sample_ms'])
    assert disabled_attempt > 0 and disabled_sample > 0
    result['instance_api_duration_ms'] = math.fsum([create, disabled_attempt, disabled_sample, *legacy_calls, *modern_calls])
    result['separate_provisioning_attribution_ms'] = provisioning
    return result
''', reader_edits)
(PHASE/'software_contracts.py').write_text(reader)

comparator = (PRIOR/'compare_software.py').read_text()
comparator_edits = []
comparator = replace(comparator,
    "    assert {stage for stage in old if new[stage]!=old[stage]}==(set() if producer=='baseline' else {'transport','endpoint'})",
    '    assert new == old  # Factory stock selections are unchanged; reader checks actual owner/stock commands separately.', comparator_edits)
comparator = replace(comparator, "result={'candidate_revision':freeze['candidate_revision']", '''assert freeze['observed_instance_api_minimum_improvement_percent'] == 0
accounting_checks = []
for direction, a_label, b_label in [('forward','a1','b1'),('reverse','a2','b2')]:
    a = records[a_label]['instance_api_duration_ms']; b = records[b_label]['instance_api_duration_ms']
    assert math.isfinite(a) and a > 0 and math.isfinite(b) and b > 0
    improvement = 100 * (a - b) / a
    accounting_checks.append({'direction':direction,'scope':'observed-instance Create and all11 original numerical API windows',
        'baseline_ms':a,'candidate_ms':b,'improvement_percent':improvement,'minimum_percent':0,'pass':improvement>=0})
result={'candidate_revision':freeze['candidate_revision']''', comparator_edits)
comparator = replace(comparator, "    'pass':all(check['pass'] for check in checks),'checks':checks,'failed_checks':[q for q in checks if not q['pass']],",
    "    'pass':all(q['pass'] for q in checks + accounting_checks),'original20_pass':all(q['pass'] for q in checks),\n"
    "    'provisioning_accounting_pass':all(q['pass'] for q in accounting_checks),'accounting_checks':accounting_checks,\n"
    "    'checks':checks,'failed_checks':[q for q in checks if not q['pass']],", comparator_edits)
(PHASE/'compare_software.py').write_text(comparator)
launcher = (PRIOR/'run_sequence.py').read_text()
launcher_edits = []
launcher = replace(launcher, 'W=Path(__file__).resolve().parent; R=W.parents[1]',
                   'W=Path(__file__).resolve().parent; R=W.parents[2]', launcher_edits)
(PHASE/'run_sequence.py').write_text(launcher)
for name in ['run_owned.py','owned_guard.py']:
    shutil.copyfile(WORK/'runtime-stream-v2'/name, PHASE/name)
shutil.copyfile(PRIOR/'expected-readback-layout.json',PHASE/'expected-readback-layout.json')

arrays = {}
header = ROOT/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
assert sha(header) == '3f7aca3372bb872934e978fb0270f2375ce8d2361555b35c8b891e42e07161cc'
for count, name, body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',header.read_text(),re.S):
    words = [int(x[:-1],0) for x in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',body)]
    assert len(words) == int(count)
    payload = struct.pack('<'+'I'*len(words),*words)
    arrays[name] = {'bytes':len(payload),'sha256':hashlib.sha256(payload).hexdigest()}
assert len(arrays) == 40
for name in ['baseline.json','build-and-arrays.json']:
    (PHASE/name).write_text(json.dumps({'arrays':arrays,'whole_header_sha256':sha(header)},indent=2)+'\n')
optional = {}
for role in ['transport','endpoint']:
    path = WORK/'full-builder/products'/('retained_'+role+'_portable_stage7_owner.spv')
    optional[role] = {'bytes':path.stat().st_size,'sha256':sha(path)}
(PHASE/'owner-module-authority.json').write_text(json.dumps(optional,indent=2)+'\n')
ledger = {'prior_reader_sha256':sha(PRIOR/'software_contracts.py'), 'prior_comparator_sha256':sha(PRIOR/'compare_software.py'),
          'prior_launcher_sha256':sha(PRIOR/'run_sequence.py'), 'reader_edits':reader_edits,
          'comparator_edits':comparator_edits,'launcher_edits':launcher_edits,
          'original20_numeric_loop_unchanged':True,'new_accounting_condition_before_data':'Both forward/reverse observed-instance API duration improvements >=0%'}
(PHASE/'edits.json').write_text(json.dumps(ledger,indent=2)+'\n')
executables = {}
for producer, revision in [('baseline','9fbd0e8b0d549a13476d066f2877dcd2cf9a94d9'),('candidate','87ae628f25ac99e8e4c412c4d21eb868753dce60')]:
    path = ROOT/'bin/linux-gcc/portable-stage-seven-owner'/('original_'+producer+'_controls')
    executables[producer] = {'path':str(path.relative_to(ROOT)),'source_revision':revision,'bytes':path.stat().st_size,'sha256':sha(path)}
frozen = [str(p.relative_to(ROOT)) for p in PHASE.iterdir() if p.is_file()]
frozen += [v['path'] for v in executables.values()]
frozen += ['bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h',
           '/usr/lib/x86_64-linux-gnu/libvulkan.so.1','/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so','/usr/share/vulkan/icd.d/lvp_icd.json',
           str((WORK/'prepare_benefit_controls.py').relative_to(ROOT))]
actions = {}
for label in ['a1','b1','b2','a2']:
    producer = 'baseline' if label.startswith('a') else 'candidate'
    # Preserve the original exact relative argv shape; cwd is the canonical checkout.
    actions[label] = {'revision':executables['candidate']['source_revision'],'kind':'test','expected_exit':0,
        'argv':[executables[producer]['path'],'--gtest_filter=RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults'],
        'limits':{'seconds':300,'rss_kib':4194304,'reserve_kib':2097152},'software_provider':True,
        'frozen_inputs':frozen+[str((PHASE/'retention-gate-frozen.json').relative_to(ROOT))],
        'scope':'One frozen original timestamp slot: unchanged32digests/51markers/three repetitions and20 numerical performance conditions, plus separately preregistered observed-instance API duration accounting; no frame/science/release claim.'}
(PHASE/'actions.json').write_text(json.dumps(actions,indent=2)+'\n')
freeze = {'candidate_revision':executables['candidate']['source_revision'],'baseline_revision':executables['baseline']['source_revision'],
    'executables':executables,'sequence':['a1','b1','b2','a2'],'within_slot_repeats':3,
    'complete_call_minimum_improvement_percent':5,'every_inclusive_stage_minimum_improvement_percent':-5,
    'all20_gates_required':True,'observed_instance_api_minimum_improvement_percent':0,
    'reader_sha256':sha(PHASE/'software_contracts.py'),'comparator_sha256':sha(PHASE/'compare_software.py'),
    'launcher_sha256':sha(PHASE/'run_sequence.py'),'actions_sha256':sha(PHASE/'actions.json'),
    'scope':'Frozen before fixedABBA invocation. Original20 numeric loop unchanged. Additional forward/reverse actual existing observed-instance API duration no-regression gates; no hidden provisioning subtraction or assumed reuse count.'}
(PHASE/'retention-gate-frozen.json').write_text(json.dumps(freeze,indent=2)+'\n')
print(json.dumps({'phase':str(PHASE.relative_to(ROOT)),'sequence':freeze['sequence'],'reader':freeze['reader_sha256']},indent=2))
