"""Explicitly adapt the preserved finite trial tools; launch nothing."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import ast, hashlib, json

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent/'runtime'
PRIOR = ROOT/'attestations/native-vulkan/99a4a6d/native-bound-selection-trial/files'
assert WORK == ROOT/'out/native-division-first-correction-screen/runtime'
def seal(path):
    data = path.read_bytes()
    return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
adaptations = []
def adapt(name, replacements):
    source = PRIOR/name;text = source.read_text();changes = []
    for old,new in replacements:
        assert text.count(old) == 1,(name,old,text.count(old))
        text = text.replace(old,new);changes.append({'old':old,'new':new})
    target = WORK/name;assert not target.exists();target.write_text(text)
    if target.suffix == '.py':ast.parse(text)
    adaptations.append({'file':name,'source':str(source),'source_seal':seal(source),
                        'output_seal':seal(target),'replacements':changes})
old_head = "HEAD='99a4a6d3a6bff69d9f20c995383a5f2e6b8cb901'"
new_head = "HEAD=json.loads((WORK/'trial-source.json').read_text())['source_revision']"
old_work = "assert WORK==ROOT/'out/native-bound-selection-trial'"
new_work = "assert WORK==ROOT/'out/native-division-first-correction-screen/runtime'"
adapt('build_owned.py',[(old_head,new_head),(old_work,new_work)])
old_helpers = '''  helpers={}
  for name,expected in (('RPUpAdd','OpFAdd'),('RPUpMultiply','OpFMul')):
   body=re.search(r'%'+name+r' = OpFunction\\b.*?OpFunctionEnd',after,re.S)[0]
   assert not re.search(r'OpBranchConditional|OpSelectionMerge|OpLoopMerge',body)
   fp=re.findall(r'(%\\S+) = (OpF(?:Add|Sub|Mul|Div)) ',body)
   assert len(fp)==1 and fp[0][1]==expected and fp[0][0] in decorated
   helpers[name]={'body':body,'fp_operations':1,'conditional_branches':0,'NoContraction':True,'uint_selects':body.count('OpSelect %uint')}
'''
new_helpers = '''  helpers={}
  division_bodies=re.findall(r'%RTDivide(?:_\\d+)? = OpFunction\\b.*?OpFunctionEnd',after,re.S)
  old_division_bodies=re.findall(r'%RTDivide(?:_\\d+)? = OpFunction\\b.*?OpFunctionEnd',before,re.S)
  assert len(division_bodies)==len(old_division_bodies)==2
  for index,(old_body,body) in enumerate(zip(old_division_bodies,division_bodies)):
   assert body.count('OpLoopMerge')==old_body.count('OpLoopMerge')==1
   assert re.search(r'OpStore %i(?:_\\d+)? %uint_1\\b',body)
   assert re.search(r'OpStore %i(?:_\\d+)? %uint_0\\b',old_body)
   assert 'DontUnroll' in body and 'DontUnroll' in old_body
   for helper,old_count,new_count in (('RTMultiply',2,3),('RTSubtract',2,3),('RTAdd',1,2),('RTCenter',1,2),('RPDivideSeed',2,3)):
    pattern=r'OpFunctionCall %\\S+ %'+helper+r'(?:_\\d+)?\\b'
    assert len(re.findall(pattern,old_body))==old_count
    assert len(re.findall(pattern,body))==new_count
   helpers[str(index)]={'baseline_body':old_body,'candidate_body':body,
                        'remaining_loop_starts_at_one':True,'static_product_sites':3,
                        'dynamic_corrections':3,'final_product_remains_separate':True}
  if stem=='retained_transport':
   assert seal(path)==seal(WORK.parent/'products/transport.spv')
'''
adapt('inspect_products.py',[(old_head,new_head),(old_work,new_work),
    ("assert changed_functions==['RPUpAdd','RPUpMultiply'] and old.keys()==new.keys(),changed_functions",
     "assert changed_functions==['RTDivide'] and old.keys()==new.keys(),changed_functions"),
    (old_helpers,new_helpers)])
adapt('prepare_native_stage.py',[
    ("assert HEAD == '99a4a6d3a6bff69d9f20c995383a5f2e6b8cb901'",
     "assert HEAD == json.loads((WORK/'trial-source.json').read_text())['source_revision']"),
    ("STAGE = ROOT / 'bin/windows-msvc/native-99a4a6d-bound-select-trial'",
     "STAGE = ROOT / ('bin/windows-msvc/native-'+HEAD[:7]+'-division-first-correction-trial')"),
    ("assert ci['databaseId'] == 38159862169 and ci['headSha'] == HEAD",
     "assert ci['databaseId'] == json.loads((WORK/'ci-admission.json').read_text())['databaseId'] and ci['headSha'] == HEAD"),
    ("'prepare_native_stage.py','native-preparation-review.json'",
     "'trial-source.json','ci-admission.json','prepare_native_stage.py','native-preparation-review.json'")])
for name in ('run_native.py','run_native_control.ps1','emergency_native_cleanup.ps1',
             'run_native_retirement.py','verify_native_retirement.ps1'):
    adapt(name,[])
adapt('compare_matched.py',[
    ('ROOT = WORKSPACE.parent.parent','ROOT = WORKSPACE.parent.parent.parent'),
    ('    owner_births, test_births = set(), set()\n',
     '    owner_births, test_births = set(), set()\n    validated_seals = {}\n'),
    ('            validate_execution(label+"-"+case)\n',
     '''            for path,expected in validate_execution(label+"-"+case).items():
                if path in validated_seals:require(validated_seals[path]==expected, "Input seal changed during comparison: "+path)
                validated_seals[path]=expected
'''),
    ('    payload = json.dumps(report, indent=2, allow_nan=False) + "\\n"\n',
     '''    for path,expected in validated_seals.items():
        require(whole_seal(Path(path) if Path(path).is_absolute() else ROOT/path)==expected, "Input changed before comparison completion: "+path)
    report['complete_input_seals']=validated_seals
    payload = json.dumps(report, indent=2, allow_nan=False) + "\\n"
''')])
adapt('analyze_controls_v2.py',[
    ("WORK/'controls-analysis.stdout',WORK/'controls-analysis.stderr',",''),
    ("report={'accepted':True",'''for name,expected in input_seals.items():assert m.whole_seal(ROOT/name)==expected,name
report={'accepted':True''')])
# Consumers use the current corrected reader name; preserve its actual source adaptation.
target = WORK/'analyze_controls_v2.py';destination = WORK/'analyze_controls.py'
assert not destination.exists();destination.write_bytes(target.read_bytes())
adapt('evaluate_retention.py',[
    ('r=w.parent.parent','r=w.parent.parent.parent'),
    ("assert ci['databaseId']==38159862169 and ci['headSha']==head and ci['status']=='completed' and ci['conclusion']=='success'",
     "assert ci['databaseId']==doc(w/'ci-admission.json')['databaseId'] and ci['headSha']==head and ci['status']=='completed' and ci['conclusion']=='success'"),
    ("'kind':'Root joined original native bound-selector retention decision'",
     "'kind':'Root joined original native first-correction peeling retention decision'")])
assert not (WORK/'helper-adaptations.json').exists()
(WORK/'helper-adaptations.json').write_text(json.dumps({'accepted':True,'preparer_seal':seal(Path(__file__)),
    'adaptations':adaptations,'reader_name_copy':{'source':'analyze_controls_v2.py','target':'analyze_controls.py','seal':seal(destination)},
    'limits':'Source preparation only. Independent source review and immutable runtime source/product/provider plan required before any owned launch.'},indent=2)+'\n')
print(json.dumps({'accepted':True,'helpers':len(adaptations),'launches':0}))
