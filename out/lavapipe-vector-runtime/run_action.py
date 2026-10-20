from pathlib import Path
import hashlib,json,os,subprocess,sys,xml.etree.ElementTree as ET
R=Path.cwd();W=Path(__file__).resolve().parent;case=sys.argv[1];folder=Path(sys.argv[2]).resolve(strict=True)
assert folder.parent==W
actions=json.loads((W/'actions.json').read_text());a=actions[case]
def sha(p):return hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
for dependency in a.get('prerequisites',[]):
 p=W/dependency/'owner.json';v=json.loads(p.read_text());assert v['passed'] and v['stop_reason'] is None and v['remaining_owned_processes']==[]
if 'corpus' in case or 'barrier16' in case:
 mode='corpus64' if 'corpus' in case else 'barrier16';local='64' if mode=='corpus64' else '16'
 B=R/'bin/linux-gcc/lavapipe-vector-validation';D=R/'out/lavapipe-vector-validation-build';actual=folder/'actual.bin'
 args=[str(B/'gpu_probe'),str(B/(mode+'.spv')),str(D/'input.bin'),str(D/'expected.bin'),str(actual),'--dispatch'+local]
 p=subprocess.run(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True);print(p.stdout,end='',flush=True);print(p.stderr,end='',file=sys.stderr,flush=True);assert p.returncode==0
 rows=[json.loads(x) for x in p.stdout.splitlines() if x.startswith('{')]
 first=[r for r in rows if 'resident_driver_path' in r];last=[r for r in rows if 'completed_dispatches' in r]
 assert len(first)==len(last)==1 and first[0]['resident_driver_path']==a['runtime_provider']
 assert last[0]['cases']==174268 and last[0]['mismatches']==last[0]['untouched_complement_words']==0
 assert last[0]['observed_words']==3*174268+48 and actual.read_bytes()==(D/'expected.bin').read_bytes()
 (folder/'corpus-acceptance.json').write_text(json.dumps({'passed':True,'mode':mode,'provider':first[0],
  'result':last[0],'actual_output_sha256':sha(actual),'complete_expected_output_sha256':sha(D/'expected.bin'),
  'scope':'All original174268 independent exact triples plus48 preserved tail words; barrier16 additionally checks cross-lane shared permutations and atomic active mask/count. No integrated or benefit qualification.'},indent=2)+'\n')
else:
 assert case in ['control-8','a1','b1','b2','a2']
 B=R/'bin/linux-gcc/lavapipe-vector-validation'
 subprocess.run([str(B/'audited_backend_tests'),'--gtest_filter=RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults'],check=True)
 v=ET.parse(folder/'gtest.xml').getroot();assert v.attrib['tests']=='1' and v.attrib['failures']==v.attrib['errors']==v.attrib.get('disabled','0')=='0'
 tests=v.findall('.//testcase');assert len(tests)==1 and tests[0].attrib['name']=='DeviceTimestampsPreserveOriginalIntervalResults'
 assert not tests[0].findall('failure') and not tests[0].findall('skipped')
