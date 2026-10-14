from pathlib import Path
import copy,hashlib,json,shutil
root=Path.cwd();prior=root/'out/forensic-kerr-calibration';work=root/'out/forensic-kerr-ported';assert not work.exists();work.mkdir()
def ident(p):
 b=p.read_bytes();return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
m=json.loads((prior/'preparation.json').read_text());pins=json.loads((prior/'pins.json').read_text())['files']
for name,row in pins.items():assert ident(Path(name))==row,name
first=prior/'legacy-kerr-original-size';owner=json.loads((first/'owner.json').read_text());assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged'] and not owner['remaining_owned_processes'] and not owner['cleanup_errors']
legacy=Path(m['output']);assert legacy.is_file()
output=root/'renders/forensic-reconstruction-b739f92/ported/kerr-png.png';assert not output.exists();output.parent.mkdir(parents=True,exist_ok=True)
for name in ('run_owned.py','owned_guard.py','run_action.py'):shutil.copy2(prior/name,work/name)
shutil.copy2(Path(__file__),work/'prepare.py')
m['role']='ported';m['output']=str(output);m['legacy_reference']={'path':str(legacy),**ident(legacy)};m['budget_scope']='One full-size ported Kerr PNG, unchanged original scene/sampling and300s diagnostic cap; no scaled substitute, original tape identity or performance qualification.';m['prior_legacy_owner']=str(first/'owner.json')
(work/'preparation.json').write_text(json.dumps(m,indent=2,sort_keys=True)+'\n')
a=copy.deepcopy(json.loads((prior/'actions.json').read_text())['legacy-kerr-original-size'])
oldcmd=m['case']['future_render_argv_review_only']['legacy'];newcmd=m['case']['future_render_argv_review_only']['ported'];assert a['argv'][-len(oldcmd):]==oldcmd;a['argv']=a['argv'][:-len(oldcmd)]+newcmd;a['scope']=m['budget_scope']
files=sorted(set(list(pins)+[str(prior/'pins.json'),str(prior/'actions.json'),str(first/'owner.json'),str(first/'live-bindings.json'),str(first/'live-maps.txt'),str(legacy)]+[str(p) for p in work.iterdir() if p.is_file()]+[str(work/'actions.json'),str(work/'pins.json')]))
a['frozen_inputs']=files;(work/'actions.json').write_text(json.dumps({'ported-kerr-original-size':a},indent=2)+'\n')
sealed={name:ident(Path(name)) for name in files if name!=str(work/'pins.json')};(work/'pins.json').write_text(json.dumps({'files':sealed,'file_count':len(sealed),'file_records_sha256':hashlib.sha256(json.dumps(sealed,sort_keys=True,separators=(',',':')).encode()).hexdigest()},indent=2,sort_keys=True)+'\n')
print(json.dumps({'prepared':True,'pins':len(sealed),'case':'ported-kerr-original-size','output':str(output),'legacy_reference':m['legacy_reference']}))
