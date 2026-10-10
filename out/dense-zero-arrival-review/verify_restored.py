from pathlib import Path
import hashlib,json,struct,subprocess,xml.etree.ElementTree as ET,re
r=Path.cwd();w=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def git(*args):return subprocess.check_output(['git',*args],text=True).strip()
head=git('rev-parse','HEAD');assert head=='d92310c18ccb6141cf5c967f61be706f4603edd9' and not git('status','--porcelain')
s=json.loads((w/'source.json').read_text());assert s['source_revision']==head and s['trial_disposition']=='rejected_restored'
b=json.loads((w/'build-and-arrays.json').read_text());a=json.loads((w/'linux-readonly-payloads.json').read_text());assert b['source_revision']==a['source_revision']==head and b['all_40_arrays_exact'] and b['changed_array_names']==[] and b['returncode']==0 and a['pass_']
assert sha(r/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h')==a['whole_header_sha256']==b['whole_header_sha256']=='3f7aca3372bb872934e978fb0270f2375ce8d2361555b35c8b891e42e07161cc'
for name,seal in a['actual_consumers'].items():assert len(seal['arrays'])==40 and (r/name).stat().st_size==seal['bytes'] and sha(r/name)==seal['sha256']
for name in ['scripts/build-retained-kernels.py','src/sirius/kernels/retained_dense.slang','tests/backend/retained_dopri_test.cpp']:
 assert (r/name).read_bytes()==subprocess.check_output(['git','show','8d5b89d:'+name])
p=r/'tests/backend/retained_compute_test.cpp';old=subprocess.check_output(['git','show','8d5b89d:tests/backend/retained_compute_test.cpp'],text=True);current=p.read_text();name='TEST_F(RetainedComputeTest, DenseSegmentsPreserveSmallCovariantArrivalDerivatives) {'
def body(x):start=x.index(name);return re.sub(r'\s','',x[start:x.index('\nTEST_F(',start+len(name))])
assert body(old)==body(current)
q=json.loads((w/'restored-first-sequence.json').read_text());assert q['pass_'] and len(q['receipts'])==3
cases=['RetainedComputeTest.DenseSegmentsPreserveSmallCovariantArrivalDerivatives','RetainedComputeTest.ZeroFractionArrivalsPreserveProgramAuthorityAndCompleteRefusal','RetainedDopriTest.SamplerPreservesPhysicalArrivalAndRejectsInconsistentRates'];details=[];births=[]
for o,case in zip(q['receipts'],cases):
 assert o['case']==case and o['execution_source_revision']==head and o['mode']=='restored' and o['pass_'] and o['returncode']==0 and o['stop_reason'] is o['owner_error'] is None and not o['owned_group_remaining'] and o['unchanged_inputs'] and o['timeout_seconds']==180 and o['maximum_group_rss_mib']==4096
 x=o['xml'];assert x['tests']==1 and x['failures']==x['errors']==x['skipped']==0 and len(x['cases'])==1 and x['cases'][0]['name']==case
 for name,seal in o['whole_input_seals'].items():assert Path(name).stat().st_size==seal['bytes'] and sha(Path(name))==seal['sha256']
 for child in o['known_births']:
  p=Path('/proc')/str(child['pid'])/'stat';assert not p.exists() or p.read_text().rsplit(')',1)[1].split()[19]!=child['start_ticks'];births.append(child)
 details.append({k:o[k] for k in ['case','elapsed_seconds','peak_sampled_rss_kib','xml']})
xml=w/'controls/restored-first-zero-arrival/tests.xml';tree=ET.parse(xml);props={p.get('name'):p.get('value') for p in tree.findall('.//property')};assert len(props)==24
row_bytes=4279*4;expected_prefix=[8,8,3,3,3,3,3,4,1,1];streams={}
words=[0]*4279;words[0]=1
for field in range(40):words[4+5*field]=0x3f800000;words[8+5*field]=1
short_sha=hashlib.sha256(struct.pack('<'+'I'*len(words),*words)).hexdigest()
for mode in ['default','fp64']:
 prefix='zero_arrival_'+mode;assert props[prefix+'_device']=='llvmpipe (LLVM 20.1.2, 256 bits)' and props[prefix+'_driver'].startswith('llvmpipe: Mesa 25.2.8-0ubuntu0.24.04.2')
 commands=[dict(v.split('=',1) for v in props[prefix+'_command_'+str(i)].split(';')) for i in range(10)]
 for i,c in enumerate(commands):assert c['stage']==('dense-short' if i==9 else 'dense') and int(c['rows'])==expected_prefix[i] and int(c['bytes'])==expected_prefix[i]*row_bytes and re.fullmatch('[0-9a-f]{64}',c['sha256'])
 assert commands[0]==commands[1]
 assert all(commands[i]['sha256']==hashlib.sha256(bytes(3*row_bytes)).hexdigest() for i in [3,4,5])
 assert commands[9]['sha256']==short_sha
 streams[mode]=commands
assert streams['default']==streams['fp64']
report={'source_revision':head,'source_tree':git('rev-parse','HEAD^{tree}'),'pass_':True,'disposition':'Trial rejected and both production owners restored; useful independent regression retained','all_40_arrays_and_complete_header_exact':True,'all_160_full_readonly_bindings':True,'original_dense_case_token_exact':True,'original_sampler_source_exact':True,'controls':details,'actual_zero_arrival_command_prefix':expected_prefix,'complete_readbacks_per_product_mode':10,'actual_two_modes_same_complete_readback_hashes':True,'zeroed_malformed_readbacks_independently_rehashed':True,'short_table_expected_whole_row_hash':short_sha,'whole_input_seals_count_per_control':[len(o['whole_input_seals']) for o in q['receipts']],'owned_births_absent':births,'scope':'Three exact-restored-source finite software controls. Original physical/polynomial arrival and rollback assertions preserved. No candidate numerical pass, matched benefit, cold attribution, native/frame/fullscientific/release claim.'}
(w/'restored-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['source_revision','pass_','disposition','actual_zero_arrival_command_prefix','complete_readbacks_per_product_mode','actual_two_modes_same_complete_readback_hashes']}))
