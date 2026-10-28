"""Read-only acceptance of one retained cleanup artifact; never repeat cleanup."""
from pathlib import Path
from collections import Counter
import hashlib
import json
import re
import shlex
import subprocess
import sys

R=Path.cwd();W=Path(__file__).resolve().parent
V=R/'out/native-cpu-strict-ir-linking';SDK=R/'bin/linux-gcc/native-cpu-strict-ir-linking';D=SDK/'toolchain/bin'
B=R/'bin/linux-gcc/native-cpu-canonical-closed-endpoint'
BC=B/'endpoint.closed.bc';LL=B/'endpoint.closed.ll';CPP=V/'audit_closed_endpoint_residual.cpp'
BC_PIN=dict(bytes=19497068,sha256='2c794bcf6a0171171366e62c45ea0011e0480045abc6dff4ff58ac5bbf746075')
LL_PIN=dict(bytes=119473814,sha256='46bcf17beabf9b74d882bc70eca31b99f70669dea90172788b674b665acecbd8')
FAILED_PIN=dict(bytes=1960,sha256='47c95bd43cbcca9174c6e9526700f1f893a59c98a4a472abae8e23f390829d40')
CONTEXT_PLAN_PIN=dict(bytes=17316,sha256='c4c82eb73c7e6947af6ab1d83d2c503991b69665c6b1bf1c54d7041f8e376813')
PLUGIN_OWNER=W/'cleanup-audit-plugin-endpoint/owner.json'
AUDIT_OWNER=W/'cleanup-audit-endpoint/owner.json'
NATIVE_OWNER=W/'native-audited-endpoint/owner.json'
PLAN=W/'endpoint-residual-plan.json';RESULT=W/'endpoint-cleanup-result.json';PROOF=W/'endpoint-residual-proof.json'
PASS='sirius-audit-closed-endpoint-residual'


def pin(p):
    with p.open('rb') as s:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(s,'sha256').hexdigest())

def doc(p):return json.loads(p.read_text())
def relative(p):return str(p.relative_to(R)) if p.is_relative_to(R) else str(p)
def row(p):return dict(path=relative(p),**pin(p))
def checked(p,identity):assert pin(p)=={k:identity[k] for k in ('bytes','sha256')},p

def write(p,value):
    with p.open('x') as s:json.dump(value,s,indent=2);s.write('\n')

def terminal(path,success=True):
    d=doc(path);assert d['status']=='terminal' and d['passed'] is success and d['child_exit']==(0 if success else 1)
    assert d['source_unchanged'] and d['whole_inputs_unchanged']
    assert not d['remaining_owned_processes'] and not d['cleanup_errors']
    assert d['observed_births'] and all(x['absent'] for x in d['observed_births'])
    assert path.with_name('inputs-before.json').read_bytes()==path.with_name('inputs-after.json').read_bytes()


def retained():
    checked(BC,BC_PIN);checked(LL,LL_PIN)
    failed=W/'cleanup-endpoint/owner.json';checked(failed,FAILED_PIN);terminal(failed,False)
    p=W/'endpoint-context-plan.json';checked(p,CONTEXT_PLAN_PIN);plan=doc(p)
    assert plan['commands']['cleanup']==[str(D/'opt'),'-passes='+plan['cleanup_pipeline'],'-verify-each','-mcpu=skylake-avx512',str(B/'endpoint.context.bc'),'-o',str(BC)]
    assert plan['outputs']['cleanup']=='endpoint.closed.bc' and plan['transactional_staging_words']==30632
    assert doc(failed)['argv']==['/usr/bin/python3','-B',str(W/'context.py'),'endpoint','cleanup']
    stderr=(failed.parent/'stderr.log').read_text()
    assert 'Unexpected Endpoint residual private storage' in stderr and 'AssertionError' in stderr
    assert 'programBase_3.i' in stderr and 'root_1.i' in stderr
    for name,identity in plan['inputs'].items():checked(R/name,identity)
    return plan


def freeze_history(inputs):
    for p in [W/'endpoint-context-plan.json',W/'endpoint-context-result.json',W/'context.py',W/'native_check.py',CPP,Path(__file__).resolve(),BC,LL]:
        inputs[relative(p)]=pin(p)
    for name in ('owner.json','inputs-before.json','inputs-after.json','stdout.log','stderr.log','samples.jsonl'):
        p=W/'cleanup-endpoint'/name;inputs[relative(p)]=pin(p)
    return inputs


def action(argv,inputs,scope):
    return dict(revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        kind='diagnostic',expected_exit=0,software_provider=False,
        limits=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152),argv=argv,
        frozen_inputs=sorted(inputs),scope=scope)


def append_actions(additions,prior_path):
    assert not subprocess.check_output(['git','status','--porcelain']).strip()
    path=W/'actions.json';before=doc(path)
    assert not set(additions)&set(before)
    with prior_path.open('xb') as s:s.write(path.read_bytes())
    after=dict(before);after.update(additions)
    assert all(after[k]==v for k,v in before.items())
    path.write_text(json.dumps(after,indent=2)+'\n')


def proof_checked(proof):
    for key in ('passed','read_only','full_module_fingerprint_unchanged','FP_call_inventory_unchanged',
        'strict_dynamic_exceptions','ordinary_entry_async_uwtable','math9_fixed_cc92_inaccessibleRW'):
        assert proof[key] is True,key
    assert proof['input_sha256']==BC_PIN['sha256'] and proof['text_sha256']==LL_PIN['sha256']
    assert (proof['definitions'],proof['owned_context_allocas'],proof['ordinary_FP_other_than_fneg'])==(1,0,0)
    assert (proof['staging_bytes'],proof['residual_bytes'],proof['residual_memory_operations'])==(122528,320,44)
    assert (proof['residual_loads'],proof['residual_stores'],proof['residual_masked_stores'],proof['residual_masked_loads'])==(17,11,16,0)
    assert proof['residual_root_zero_memsets']==1 and len(proof['residuals'])==6
    expected={'programBase_3.i':32,'count_2.i':32,'layers_2.i':32,'prefix_5.i':32,'component_4.i':32,'root_1.i':160}
    assert {x['name']:x['bytes'] for x in proof['residuals']}==expected
    assert all(x['alignment']==32 and x['full_initialization_dominates_every_memory_use'] and x['no_capture_or_dynamic_address'] for x in proof['residuals'])
    expected_math={'sirius_local_value_dense_multiply8':12,'sirius_local_value_dense_add8':11,
        'sirius_local_value_dense_divide8':6,'sirius_local_value_dense_sqrt8':1,
        'sirius_local_value_tagged_multiply8':1963,'sirius_local_value_tagged_divide8':154,
        'sirius_local_value_tagged_add8':1607,'sirius_local_value_tagged_sqrt8':3,'sirius_local_value_tagged_subtract8':197}
    expected_FP={'llvm.experimental.constrained.fadd.v8f32':96,'llvm.experimental.constrained.fsub.v8f32':13,
        'llvm.experimental.constrained.fcmp.v8f32':451,'llvm.experimental.constrained.fcmps.v8f32':28,
        'llvm.experimental.constrained.fmul.v8f32':1}
    assert proof['value_calls']==expected_math and proof['value_call_total']==sum(expected_math.values())==3954
    assert proof['constrained_calls']==expected_FP and proof['ordinary_fneg']==16 and proof['fabs_calls']==312
    assert proof['numerical_qualification_transferred'] is False


phase,=sys.argv[1:]
if phase=='prepare-audit':
    old=retained();assert not RESULT.exists() and not PROOF.exists() and not PLAN.exists()
    assert not (B/'endpoint.residual-audit.so').exists()
    config=shlex.split(subprocess.check_output([str(D/'llvm-config'),'--cxxflags'],text=True))
    commands={
        'plugin':[str(D/'clang++'),*config,'-O2','-fPIC','-shared',str(CPP),str(SDK/'sha256-only.a'),'-Wl,--exclude-libs,ALL','-o',str(B/'endpoint.residual-audit.so')],
        'audit':[str(D/'opt'),'-load-pass-plugin='+str(B/'endpoint.residual-audit.so'),'-passes='+PASS,'-verify-each','-disable-output',str(BC)]}
    inputs=freeze_history(dict(old['inputs']))
    for p in (V/'toolchain.json',SDK/'sha256-only.a',D/'llvm-config',D/'clang++',D/'opt'):
        inputs[relative(p)]=pin(p)
    tools=doc(V/'toolchain.json')['tools']
    for name in ('llvm-config','clang++','opt'):
        tool=tools[name];checked(R/tool['resolved'],tool);inputs[tool['resolved']]={k:tool[k] for k in ('bytes','sha256')}
    original_actions=W/'endpoint-residual-prior-actions.json'
    original_bytes=(W/'actions.json').read_bytes()
    original_snapshot=doc(W/'cleanup-endpoint/inputs-before.json')[relative(W/'actions.json')]
    assert dict(bytes=len(original_bytes),sha256=hashlib.sha256(original_bytes).hexdigest())=={k:original_snapshot[k] for k in ('bytes','sha256')}
    write(PLAN,dict(inputs=inputs,commands=commands,original_cleanup_argv=old['commands']['cleanup'],
        original_native_argv=old['commands']['native'],original_failed_owner=FAILED_PIN,original_context_plan=CONTEXT_PLAN_PIN,
        retained_BC=dict(path=relative(BC),**BC_PIN),retained_IR=dict(path=relative(LL),**LL_PIN),
        accepted_artifact_owner=relative(AUDIT_OWNER),cleanup_repeated=False,source_module_rewritten=False,
        scope='Correct only the unjustified zero-extra-allocas expectation. Exact320B bounded original metadata/projection locals are audited read-only; complete Context removal, staging, strict FP/math/masks/source/limits remain required. No numerical transfer.'))
    frozen=set(inputs)|{relative(PLAN),relative(original_actions)}
    additions={
        'cleanup-audit-plugin-endpoint':action(['/usr/bin/python3','-B',str(Path(__file__).resolve()),'audit-plugin'],frozen,
            'Compile one small read-only LLVM23 audit plugin only; do not rerun cleanup, source/IR conversion, native or workload.'),
        'cleanup-audit-endpoint':action(['/usr/bin/python3','-B',str(Path(__file__).resolve()),'audit'],
            frozen|{relative(B/'endpoint.residual-audit.so'),relative(PLUGIN_OWNER),relative(W/'endpoint-residual-plugin-result.json')},
            'Read-only disable-output audit of exact retained BC and independently parsed LL. Require only six named initialized bounded locals320B plus exact122528B staging, Context0, strict ABI/FP/math and unchanged full module. Preserve failed owner and artifacts; no compiler-pipeline rerun/numerical transfer.')}
    append_actions(additions,original_actions)
    print(json.dumps(dict(prepared=list(additions),cleanup_repeated=False,native_prepared=False)))
elif phase in ('audit-plugin','audit'):
    retained();plan=doc(PLAN)
    assert plan['original_failed_owner']==FAILED_PIN and plan['original_context_plan']==CONTEXT_PLAN_PIN
    assert not plan['cleanup_repeated'] and not plan['source_module_rewritten']
    for name,identity in plan['inputs'].items():checked(R/name,identity)
    if phase=='audit-plugin':
        assert not (B/'endpoint.residual-audit.so').exists()
        subprocess.run(plan['commands']['plugin'],check=True)
        write(W/'endpoint-residual-plugin-result.json',dict(passed=True,output=row(B/'endpoint.residual-audit.so'),argv=plan['commands']['plugin'],read_only_audit_plugin=True))
    else:
        terminal(PLUGIN_OWNER);plugin=doc(W/'endpoint-residual-plugin-result.json')
        assert plugin['passed'];checked(B/'endpoint.residual-audit.so',plugin['output'])
        assert not RESULT.exists() and not PROOF.exists()
        completed=subprocess.run(plan['commands']['audit'],check=True,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        if completed.stderr:print(completed.stderr,file=sys.stderr,end='')
        proof=json.loads(completed.stdout);proof_checked(proof)
        retained()  # Both compiler artifacts and complete failed history stay exact.
        text=LL.read_text();allocas=[x.strip() for x in re.findall(r'^\s+%[^=]+ = alloca .*$',text,re.M)]
        assert len(allocas)==7 and allocas.count('%staged_output = alloca [30632 x i32], align 64')==1
        constrained=dict(Counter(re.findall(r'\bcall\s+.*?@(llvm.experimental.constrained.[^ (]+)\(',text)))
        calls=dict(Counter(re.findall(r'\bcall\s+.*?@(sirius_local_value_[^ (]+)\(',text)))
        assert constrained==proof['constrained_calls'] and calls==proof['value_calls']
        write(PROOF,proof)
        write(RESULT,dict(passed=True,phase='cleanup',stage='endpoint',output=row(BC),ir=row(LL),
            argv=plan['original_cleanup_argv'],original_cleanup_argv=plan['original_cleanup_argv'],
            read_only_audit_argv=plan['commands']['audit'],accepted_artifact_owner=relative(AUDIT_OWNER),
            original_failed_owner=FAILED_PIN,original_context_plan=CONTEXT_PLAN_PIN,
            cleanup_repeated=False,source_module_rewritten=False,residual_audit=row(PROOF),
            definitions=1,constrained_calls=constrained,value_calls=calls,owned_context_allocas=0,
            memory_instruction_counts=dict(Counter(re.findall(r'^\s+(?:%[^=]+ = )?(load|store|alloca)\b',text,re.M))),
            masked_store_calls=text.count('call void @llvm.masked.store'),allocas=allocas,
            original_transactional_staging_bytes=122528,numerical_qualification_transferred=False,
            scope='Existing optimizer artifact accepted only through separate pinned read-only residual audit; original cleanup owner remains failed. Exact six320B original selector/projection locals are permitted after proof; Context0/staging/strict FP/math/source/mask/transaction contracts remain unchanged. No optimizer rerun or numerical acceptance transfer.'))
        print(json.dumps(dict(passed=True,accepted_output=row(BC),residual_bytes=320,cleanup_repeated=False)))
elif phase=='prepare-native-audited':
    old=retained();terminal(AUDIT_OWNER);accepted=doc(RESULT);proof=doc(PROOF);proof_checked(proof)
    assert accepted['accepted_artifact_owner']==relative(AUDIT_OWNER) and accepted['original_failed_owner']==FAILED_PIN
    assert accepted['output']==row(BC) and accepted['ir']==row(LL) and accepted['residual_audit']==row(PROOF)
    assert not accepted['cleanup_repeated'] and not (B/'endpoint.closed.o').exists()
    assert not (W/'native-endpoint').exists() and not (W/'native-audited-endpoint').exists()
    inputs=freeze_history(dict(doc(PLAN)['inputs']))
    for p in (PLAN,RESULT,PROOF,PLUGIN_OWNER,AUDIT_OWNER,W/'endpoint-residual-plugin-result.json',B/'endpoint.residual-audit.so'):
        inputs[relative(p)]=pin(p)
    for owner in (PLUGIN_OWNER,AUDIT_OWNER):
        for name in ('inputs-before.json','inputs-after.json'):inputs[relative(owner.with_name(name))]=pin(owner.with_name(name))
    native_plan=W/'endpoint-native-audited-plan.json'
    write(native_plan,dict(inputs=inputs,argv=old['commands']['native'],accepted_cleanup_owner=relative(AUDIT_OWNER),
        accepted_cleanup=row(RESULT),original_failed_cleanup_owner=FAILED_PIN,
        original_unrun_native_action=doc(W/'actions.json')['native-endpoint'],
        accepted_artifact_owner=relative(NATIVE_OWNER),cleanup_repeated=False,native_emitted_before=False))
    prior=W/'endpoint-native-audited-prior-actions.json'
    addition={'native-audited-endpoint':action(['/usr/bin/python3','-B',str(Path(__file__).resolve()),'native-audited'],
        set(inputs)|{relative(native_plan),relative(prior)},
        'First native emission only, exact original llc argv after separate read-only retained-cleanup acceptance. Original failed owner/compiler outputs and unrun native action preserved; no cleanup rerun or numerical transfer. Actual noFMA/C ABI/FDE checks required.')}
    append_actions(addition,prior);print(json.dumps(dict(prepared='native-audited-endpoint',native_argv=old['commands']['native'])))
elif phase=='native-audited':
    old=retained();terminal(AUDIT_OWNER);plan=doc(W/'endpoint-native-audited-plan.json')
    for name,identity in plan['inputs'].items():checked(R/name,identity)
    assert plan['argv']==old['commands']['native'] and plan['accepted_cleanup']==row(RESULT)
    assert plan['original_failed_cleanup_owner']==FAILED_PIN and plan['accepted_artifact_owner']==relative(NATIVE_OWNER)
    assert doc(W/'actions.json')['native-endpoint']==plan['original_unrun_native_action']
    assert not (W/'native-endpoint').exists() and not (B/'endpoint.closed.o').exists()
    assert not (W/'endpoint-native-result.json').exists()
    subprocess.run(plan['argv'],check=True)
    from native_check import check
    result=dict(phase='native',stage='endpoint',output=row(B/'endpoint.closed.o'),argv=plan['argv'],
        accepted_artifact_owner=relative(NATIVE_OWNER),accepted_cleanup_owner=relative(AUDIT_OWNER),
        original_failed_cleanup_owner=FAILED_PIN,cleanup_repeated=False,numerical_qualification_transferred=False)
    result.update(check(B/'endpoint.closed.o','sirius_closed_endpoint_packet'))
    result['cleaned_IR_value_calls']=doc(RESULT)['value_calls']
    write(W/'endpoint-native-result.json',result)
    print(json.dumps(dict(phase='native',stage='endpoint',output=result['output'],instructions=result['instructions'])))
else:raise ValueError(phase)
