"""Remove only the two empty, owned frame temporary directories after recovery push."""
from pathlib import Path
import hashlib,json,subprocess
R=Path.cwd();W=R/'out/native-cpu-local-identities-cohort/cold-frame';J=lambda p:json.loads(p.read_text())
assert not subprocess.check_output(['git','status','--porcelain']).strip()
assert subprocess.check_output(['git','rev-parse','HEAD'])==subprocess.check_output(['git','rev-parse','@{upstream}'])
o=J(W/'frame/owner.json');assert o['status']=='terminal' and not o['remaining_owned_processes'] and not o['cleanup_errors']
assert o['source_unchanged'] and o['whole_inputs_unchanged'] and all(x['absent'] for x in o['observed_births'])
bundles=list((R/'attestations/qualification'/o['source']['revision'][:8]).glob('native-cpu-local-identities-*/MANIFEST.json'));assert len(bundles)==1
p=bundles[0];m=J(p)
for path in [p,*(p.parent/a['file'] for a in m['archives'])]:
 assert subprocess.check_output(['git','show','@{upstream}:'+str(path.relative_to(R))])==path.read_bytes()
action=J(W/'actions.json')['frame'];temp=Path(action['build_environment']['TMPDIR']);parent=temp.parent
assert parent==R/('renders/native-cpu-local-identities-cohort-'+o['source']['revision'][:8])
assert temp.name=='cold-frame-tmp' and temp.is_dir() and not temp.is_symlink() and list(temp.iterdir())==[]
assert parent.is_dir() and not parent.is_symlink() and list(parent.iterdir())==[temp]
record=W/'temporary-retirement.json';assert not record.exists()
report=dict(status='started',source=o['source'],verified_pushed_recovery=str(p.relative_to(R)),removed_empty_directories=[],scope='Only the two verified empty task-owned render temporary directories. No render payload or unique evidence removed; original resource cache preserved.')
def save():record.write_text(json.dumps(report,indent=2)+'\n')
save()
for path in [temp,parent]:
 path.rmdir();report['removed_empty_directories'].append(str(path.relative_to(R)));save()
report['status']='complete';save();print(json.dumps(report))
