from pathlib import Path
import hashlib,json,shutil,subprocess
r=Path.cwd();w=Path(__file__).resolve().parent;a=r/'attestations/software-vulkan/d92310c/dense-zero-arrival-rejected'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()=='d92310c18ccb6141cf5c967f61be706f4603edd9' and not subprocess.check_output(['git','status','--porcelain'])
assert json.loads((w/'restored-verification.json').read_text())['pass_'] and json.loads((w/'source-preservation.json').read_text())['remote_exact']
files={};excluded=[]
a.mkdir(parents=True,exist_ok=False)
for p in sorted(w.rglob('*')):
 if not p.is_file():continue
 rel=p.relative_to(w)
 if str(rel)=='baseline-8d5b89d/sirius_backend_tests':excluded.append({'path':str(p),'reason':'Prepared unused copy of preserved8d ELF; baseline benefit/Sampler was neverstarted. Original protected8d archive remains.'});continue
 if p.is_symlink():raise AssertionError(str(p))
 t=a/rel;t.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,t);assert p.stat().st_size==t.stat().st_size and sha(p)==sha(t);files[str(rel)]={'bytes':t.stat().st_size,'sha256':sha(t),'origin':str(p)}
inputs={}
seq=json.loads((w/'restored-first-sequence.json').read_text())
for o in seq['receipts']:
 for name,seal in o['whole_input_seals'].items():
  p=Path(name);assert p.stat().st_size==seal['bytes'] and sha(p)==seal['sha256'];inputs[name]=seal
b=json.loads((w/'linux-readonly-payloads.json').read_text())
for name,seal in b['actual_consumers'].items():inputs[str((r/name).resolve())]={'bytes':seal['bytes'],'sha256':seal['sha256']}
for name,seal in inputs.items():
 p=Path(name);rel=Path('restored-executed-inputs')/str(p.resolve()).lstrip('/');t=a/rel;t.parent.mkdir(parents=True,exist_ok=True);assert p.stat().st_size==seal['bytes'] and sha(p)==seal['sha256'];shutil.copy2(p,t);assert sha(t)==seal['sha256'];files[str(rel)]={'bytes':seal['bytes'],'sha256':seal['sha256'],'origin':name}
manifest={'source_revision':'d92310c18ccb6141cf5c967f61be706f4603edd9','files':files,'excluded_workspace_files':excluded,'workspace_payloads':len(files)-len(inputs),'restored_inputs':len(inputs),'total_payload_bytes':sum(v['bytes'] for v in files.values()),'scope':'Rejected Dense trial and exact restored-source finite controls. Local whole-byte preservation of task scratch and sealed inputs, including both incomplete guards and all restored modes/readbacks; no raw-evidence remote backup or full scientific/native/frame/release claim.'}
(a/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');assert {str(p.relative_to(a)) for p in a.rglob('*') if p.is_file()}==set(files)|{'manifest.json'}
for name,seal in files.items():p=a/name;assert p.stat().st_size==seal['bytes'] and sha(p)==seal['sha256']
print(json.dumps({'archive':str(a),'payloads':len(files),'bytes':manifest['total_payload_bytes'],'manifest_sha256':sha(a/'manifest.json'),'scope':manifest['scope']}))
