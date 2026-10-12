"""Freeze original20 plus inclusive API-cost gates before the sole ABBA."""
from pathlib import Path
import hashlib,json,shutil,subprocess
r=Path.cwd();w=Path(__file__).resolve().parent;p=w/'benefit';p.mkdir();b=r/'bin/linux-gcc/portable-seed-scaling';old=r/'attestations/software-vulkan/dcd1368/sealed-stage-seven-owner/files/out/sealed-stage-seven-owner/benefit';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();edits=[]
def change(s,a,c):
 assert s.count(a)==1,(a,s.count(a));edits.append({'old':a,'new':c});return s.replace(a,c)
s=(old/'software_contracts.py').read_text()
start=s.index('    expected_tables = ');end=s.index('    expected_executed = set()',start)
s=change(s,s[start:end],'''    assert not any(key.startswith(('actual_table_verification_', 'optional_module_')) for key in props)
''')
s=change(s,"            if producer == 'candidate' and role == 'transport': index = 9\n            if producer == 'candidate' and role == 'endpoint' and (':rows=' not in stage_text or stage_text.endswith(':rows=24')): index = 10\n",'')
s=change(s,"            authority = optional[role] if index >= 7 else {'bytes':selected[role]['bytes'], 'sha256':selected[role]['sha256']}","            authority = {'bytes':selected[role]['bytes'], 'sha256':selected[role]['sha256']}")
s=change(s,"            assert module['bindings'] == ('3' if index >= 7 else '2')","            assert module['bindings'] == '2'")
s=change(s,"    assert (provisioning > 0) == (producer == 'candidate')","    assert provisioning == 0")
(p/'software_contracts.py').write_text(s)
s=(old/'compare_software.py').read_text();numeric=s[s.index('checks=[]'):s.index("result={'candidate_revision'")]
s=change(s,'    assert new == old  # Factory stock selections are unchanged; reader checks actual owner/stock commands separately.',"    assert {stage:new[stage]['array'] for stage in new} == {stage:old[stage]['array'] for stage in old}\n    assert new == records['a1' if producer == 'baseline' else 'b1']['selected_modules']  # Each producer joins its exact compiled-array authority in the reader.")
assert numeric==s[s.index('checks=[]'):s.index("result={'candidate_revision'")]
(p/'compare_software.py').write_text(s)
for n in ('run_sequence.py','expected-readback-layout.json'):shutil.copyfile(old/n,p/n)
for n in ('run_owned.py','owned_guard.py'):shutil.copyfile(w/n,p/n)
shutil.copyfile(w/'baseline-arrays.json',p/'baseline.json');shutil.copyfile(w/'candidate-arrays.json',p/'build-and-arrays.json')
(p/'edits.json').write_text(json.dumps({'original_numeric20_and2_accounting_loop_byte_exact':True,'original_comparator_sha256':sha(old/'compare_software.py'),'original_reader_sha256':sha(old/'software_contracts.py'),'edits':edits},indent=2)+'\n')
rev=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip();executables={}
for mode,revision in [('baseline','9fbd0e8b0d549a13476d066f2877dcd2cf9a94d9'),('candidate',rev)]:
 exe=b/('original_'+mode+'_controls');executables[mode]={'path':str(exe.relative_to(r)),'source_revision':revision,'bytes':exe.stat().st_size,'sha256':sha(exe)}
frozen=[str(x.relative_to(r)) for x in p.iterdir() if x.is_file()]+[x['path'] for x in executables.values()]+[str(b/'products/retained_kernels.h'),str(r/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'),'/usr/lib/x86_64-linux-gnu/libvulkan.so','/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so','/usr/share/vulkan/icd.d/lvp_icd.json',str(w/'prepare_benefit.py')]
act={}
for label in ('a1','b1','b2','a2'):
 mode='baseline' if label.startswith('a') else 'candidate';act[label]={'revision':rev,'kind':'test','expected_exit':0,'argv':[executables[mode]['path'],'--gtest_filter=RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults'],'limits':{'seconds':300,'rss_kib':4194304,'reserve_kib':2097152},'software_provider':True,'frozen_inputs':frozen+[str(p/'retention-gate-frozen.json')],'scope':'One fixed original timestamp slot, unchanged32completebuffers/51orderedmarkers/54actualcommands, original20benefit+2inclusiveCreate/all11APIwindows no-regression conditions. No cost subtraction/retry/pooled order.'}
(p/'actions.json').write_text(json.dumps(act,indent=2)+'\n')
freeze={'candidate_revision':rev,'baseline_revision':executables['baseline']['source_revision'],'executables':executables,'sequence':['a1','b1','b2','a2'],'within_slot_repeats':3,'complete_call_minimum_improvement_percent':5,'every_inclusive_stage_minimum_improvement_percent':-5,'all20_gates_required':True,'observed_instance_api_minimum_improvement_percent':0,'reader_sha256':sha(p/'software_contracts.py'),'comparator_sha256':sha(p/'compare_software.py'),'launcher_sha256':sha(p/'run_sequence.py'),'actions_sha256':sha(p/'actions.json'),'scope':'Frozen before soleABBA. Numeric20+2loop byte-exact. Actual compiled authority perproducer; stock7modules2bindings noowner/table readings; all original words/orders/costs required.'}
(p/'retention-gate-frozen.json').write_text(json.dumps(freeze,indent=2)+'\n');print(json.dumps(freeze))
