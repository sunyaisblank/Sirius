"""Recover every consumed input version using exact archived bytes or Git blobs."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import resource
import subprocess
import sys

sys.dont_write_bytecode=True
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
resource.setrlimit(resource.RLIMIT_FSIZE,(256*1024*1024,256*1024*1024))
root=Path.cwd()
archive=root/'attestations/software-vulkan/661e157/portable-limb-accumulation'
bundle=json.loads((archive/'whole-byte-manifest.json').read_text())

def identity(path):
    value=hashlib.sha256()
    with path.open('rb') as stream:
        for data in iter(lambda:stream.read(1024*1024),b''): value.update(data)
    return {'bytes':path.stat().st_size,'sha256':value.hexdigest()}

aliases={}
for name,value in bundle['files'].items():
    aliases.setdefault((value['bytes'],value['sha256']),[]).append('files/'+name)
records={}
receipt_count=0
for path in sorted((archive/'files/out/portable-execution-attribution').rglob('inputs-before.json')):
    before=json.loads(path.read_text())
    after=json.loads(path.with_name('inputs-after.json').read_text())
    assert before==after,path
    owner=json.loads(path.with_name('owner.json').read_text())
    revision=owner['source']['revision']
    assert revision in ('95c4293f598c1dbf4cd13ece7b06ffff5383b81a','661e157b53f9f9fbbba0347cfda32863a1684c75')
    receipt_count+=1
    for name,value in before.items():
        key=(name,value['bytes'],value['sha256'])
        records.setdefault(key,{'identity':value,'source_revisions':set()})['source_revisions'].add(revision)
assert receipt_count==55
tracked={}
for revision in {rev for item in records.values() for rev in item['source_revisions']}:
    tracked[revision]={p.decode() for p in subprocess.check_output(['git','ls-tree','-r','--name-only','-z',revision]).split(b'\0') if p}
snapshots=archive/'retained-input-snapshots'
snapshots.mkdir()
counts=Counter()
proof=[]
for (name,size,digest),item in sorted(records.items()):
    expected={'bytes':size,'sha256':digest}
    path=Path(name) if Path(name).is_absolute() else root/name
    relative=str(path.relative_to(root)) if path.is_relative_to(root) else None
    key=(size,digest)
    if key in aliases:
        binding=aliases[key]
        assert identity(archive/binding[0])==expected
        kind='archive'
    else:
        matched=None
        for revision in sorted(item['source_revisions']):
            if relative not in tracked[revision]: continue
            data=subprocess.check_output(['git','show',revision+':'+relative])
            if len(data)==size and hashlib.sha256(data).hexdigest()==digest:
                matched={'source_revision':revision,'path':relative,
                    'blob':subprocess.check_output(['git','rev-parse',revision+':'+relative]).decode().strip()}
                break
        if matched is not None:
            kind='git'
            binding=matched
        else:
            assert relative is None or not relative.startswith(tuple(v+'/' for v in bundle['roots'])),(name,'owned input lacks archive recovery')
            assert path.is_file() and identity(path)==expected,(name,'retained input mismatch')
            kind='retained'
            binding={'path':str(path),'resolved_path':str(path.resolve()),**expected}
            if relative is None:
                target=snapshots/digest
                if not target.exists():
                    with path.open('rb') as source,target.open('wb') as output:
                        for data in iter(lambda:source.read(1024*1024),b''): output.write(data)
                assert identity(target)==expected
                binding['archive_copy']=str(target.relative_to(archive))
    counts[kind]+=1
    proof.append({'recorded_name':name,**expected,'recovery_kind':kind,'binding':binding,
        'consumed_source_revisions':sorted(item['source_revisions'])})
result={'passed':True,'recorded_input_versions':len(proof),'receipt_pairs':receipt_count,
    'recovery_counts':dict(counts),'records':proof,
    'scope':'Every recorded admitted input version has exact local archive, Git or retained-byte recovery. '
        'Unrecorded system dependencies and remote raw backup are not asserted.'}
(archive/'recorded-input-recovery.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({key:result[key] for key in ('passed','recorded_input_versions','receipt_pairs','recovery_counts')}))
