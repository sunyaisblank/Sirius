from pathlib import Path
import hashlib,json,re,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/native-retained-ownership/providers'
cases={'control-8':'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults',
'control-1':'RetainedComputeTest.FactoredEndpointsPreserveIndependentRootsAndBoundaryAdmission',
'control-2':'RetainedComputeTest.MixedIndependentTransportLayersPreserveWordsAndRefusal',
'control-3':'RetainedComputeTest.ProjectedEndpointsKeepPhysicalColumnsAndRetainedContinuation',
'control-4':'RetainedDopriTest.CurvedTransportConnectsSamplerAndProjection',
'control-5':'RetainedDopriTest.ConnectedIntervalsPreserveTrialsAndIndependentBudgets',
'control-6':'RetainedDopriTest.SamplerPreservesPhysicalArrivalAndRejectsInconsistentRates',
'control-7':'RetainedComputeTest.CoupledIntervalsRequireEmbeddedAndIndependentDenseAgreement'}
inputs={str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}|{str(p.relative_to(R)) for p in B.rglob('*') if p.is_file()}
inputs|={'bin/linux-gcc/native-original-input-ownership/matched/audited_backend_tests','out/native-original-input-ownership/providers.json','out/native-original-input-ownership/ordinary/ordinary-code-review.json','out/native-original-input-ownership/reused/prior-proof/lavapipe-ieee32-validation/candidate-default/numerical-and-boundary-result.json','out/native-original-input-ownership/reused/prior-proof/lavapipe-ieee32-validation/candidate-poison/numerical-and-boundary-result.json','/usr/bin/python3','/etc/ld.so.cache'}
inputs|={str(p.relative_to(R)) for p in (R/'bin/linux-gcc/src/sirius/backend/retained').iterdir() if p.is_file()}
inputs|={str(p.relative_to(R)) for p in (W.parent/'readers').iterdir() if p.is_file()}
inputs|={str(p.relative_to(R)) for p in (W.parent/'provider-history').rglob('*') if p.is_file()}
inputs|={'out/native-original-input-ownership/previous-control-8-acceptance.json','out/native-original-input-ownership/provider-reuse.json','out/native-original-input-ownership/compiler-result.json','out/native-original-input-ownership/matched-result.json','out/native-original-input-ownership/headers/retained_kernels.h','bin/linux-gcc/native-original-input-ownership/modules/owned-endpoint.spv','bin/linux-gcc/native-original-input-ownership/modules/owned-transport.spv','bin/linux-gcc/native-retained-ownership/modules/owned-dopri_phase.spv'}
for p in [R/'bin/linux-gcc/native-original-input-ownership/matched/audited_backend_tests',B/'baseline/libvulkan_lvp.so',B/'candidate/libvulkan_lvp.so']:
    s=subprocess.check_output(['/usr/bin/ldd',str(p)],text=True);inputs|={str(Path(x).resolve(strict=True)) for x in re.findall(r'(/[^\s()]+)',s) if Path(x).is_file()}
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip();a={}
for label,arm in [('control-8','candidate'),('a1','baseline'),('b1','candidate'),('b2','candidate'),('a2','baseline')]+[(f'control-{i}','candidate') for i in range(1,8)]:
    prior=[] if label=='control-8' else (['control-8']+['a1','b1','b2','a2'][:['a1','b1','b2','a2'].index(label)] if label in ['a1','b1','b2','a2'] else ['control-8','a1','b1','b2','a2']+[f'control-{i}' for i in range(1,int(label[-1]))])
    a[label]={'revision':revision,'arm':arm,'case':cases.get(label,cases['control-8']),'prerequisites':prior,'kind':'test','expected_exit':0,'limits':{'seconds':300,'rss_kib':4194304,'reserve_kib':2097152},'argv':['/usr/bin/python3','-B',str(W/'run_case.py'),label,'{folder}'],'frozen_inputs':sorted(inputs),'software_provider':True,'runtime_icd':str(B/arm/'lvp_icd.json'),'runtime_provider':str(B/arm/'libvulkan_lvp.so'),'scope':'Original '+label+' with the same three-module ownership candidate audited public-factory ELF, preserving all35 original shipped math modules. Actual matched Mesa25.2.8/LLVM21 providers select existing portable vs native32 route; native32 Transport owns all46 normalized inputs, Endpoint preserves the exact ordered valid prefix through its first-invalid join and its single complete projection, and DopriPhase distributes all362 original checks. All original workgroups, programs, arithmetic and assertions retained. Same work, inputs, default product mode, counts, caches, governor, limits and assertions. C8 before sole ABBA22; remaining7 only after22 PASS. No performance/adoption/frame/science/release inferred from one control.'}
(W/'actions.json').write_text(json.dumps(a,indent=2)+'\n');print(json.dumps({'frozen_inputs':len(inputs),'new_CPP_compilations':0,'same_ELF_both_arms':True,'fixed_order':['a1','b1','b2','a2']}))
