"""Copy the retired trial exactly; retain live source until recovery is reviewed."""
from pathlib import Path
import hashlib
import json
import resource
import shutil
import stat
import subprocess
import sys

sys.dont_write_bytecode=True
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
resource.setrlimit(resource.RLIMIT_FSIZE,(256*1024*1024,256*1024*1024))
root=Path.cwd()
work=Path(__file__).resolve().parent
build=root/'bin/linux-gcc/portable-indexed-limb-accumulation'
archive=root/'attestations/software-vulkan/6014f91/portable-indexed-limb-accumulation'
assert sys.argv[1:]==['--archive'] and not archive.exists()
assert subprocess.check_output(['git','rev-parse','HEAD']).decode().strip()=='5b11d31b813b47a55117f5bccff5119683e9a7a0'
assert not subprocess.check_output(['git','status','--porcelain'])
assert not subprocess.check_output(['git','diff','c8199bd','HEAD'])

def identity(path):
    value=hashlib.sha256()
    with path.open('rb') as stream:
        for data in iter(lambda:stream.read(1024*1024),b''): value.update(data)
    return {'bytes':path.stat().st_size,'sha256':value.hexdigest()}

receipts=[]
for path in sorted(work.rglob('owner.json')):
    owner=json.loads(path.read_text())
    assert owner['status']=='terminal' and owner['remaining_owned_processes']==[]
    assert not owner['cleanup_errors'] and all(v['absent'] for v in owner['observed_births'])
    receipts.append(str(path.relative_to(root)))
assert len(receipts)==38
for name,pin in json.loads((work/'protected-before.json').read_text()).items():
    assert identity(root/name)==pin,name
assert json.loads((work/'benefit/software-comparison.json').read_text())['pass'] is False
disposition=json.loads((work/'disposition.json').read_text())
assert disposition['disposition']=='REJECT_FOR_ADOPTION' and disposition['all_owned_execution_terminal']
archive.mkdir(parents=True)
manifest={}
modes={}
for owned in (work,build):
    assert owned.is_dir() and not owned.is_symlink()
    for path in sorted(owned.rglob('*')):
        mode=path.lstat().st_mode
        assert not stat.S_ISLNK(mode) and (stat.S_ISREG(mode) or stat.S_ISDIR(mode)),path
        if not stat.S_ISREG(mode): continue
        name=str(path.relative_to(root))
        target=archive/'files'/name
        target.parent.mkdir(parents=True,exist_ok=True)
        expected=identity(path)
        shutil.copy2(path,target)
        assert identity(path)==expected and identity(target)==expected,name
        manifest[name]=expected
        modes[name]=oct(mode&0o777)
        assert oct(target.stat().st_mode&0o777)==modes[name]
record={'source_revision':'6014f91d28ddeadbc7b916f213be7f89a4d1f16c',
    'restored_revision':'5b11d31b813b47a55117f5bccff5119683e9a7a0',
    'roots':[str(p.relative_to(root)) for p in (work,build)],
    'files':manifest,'preserved_file_modes':modes,'file_count':len(manifest),
    'total_bytes':sum(v['bytes'] for v in manifest.values()),'terminal_owner_receipts':receipts,
    'canonical_protected_files':114,
    'meaning':'Exact local raw recovery for a rejected trial. No qualification, adoption or remote binary backup claim.'}
(archive/'whole-byte-manifest.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'archive':str(archive.relative_to(root)),'files':len(manifest),
    'bytes':record['total_bytes'],'terminal_owner_receipts':len(receipts)}))
