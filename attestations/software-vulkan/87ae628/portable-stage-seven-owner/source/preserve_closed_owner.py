from pathlib import Path
import hashlib, json, os, shutil, subprocess

ROOT = Path.cwd()
WORK = ROOT / 'out/portable-stage-seven-owner-construction'
BIN = ROOT / 'bin/linux-gcc/portable-stage-seven-owner'
ARCHIVE = ROOT / 'attestations/software-vulkan/87ae628/portable-stage-seven-owner'
assert not ARCHIVE.exists()
ARCHIVE.mkdir(parents=True)
sha = lambda data: hashlib.sha256(data).hexdigest()

def query(args, **kwargs):
    return subprocess.run(args, check=True, capture_output=True, **kwargs).stdout

assert query(['git','status','--porcelain']).strip() == b''
manifest = {}
for folder in [WORK, BIN]:
    for path in sorted(folder.rglob('*')):
        if not path.is_file(): continue
        assert not path.is_symlink(), path
        relative = str(path.relative_to(ROOT))
        destination = ARCHIVE/'files'/relative
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,destination)
        data = path.read_bytes()
        assert data == destination.read_bytes()
        manifest[relative] = {'bytes':len(data),'sha256':sha(data)}
(ARCHIVE/'whole-byte-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

source_paths = []
for path in WORK.rglob('*.py'):
    if path.name.startswith('baseline_'): continue
    source_paths.append(path)
source_paths += [WORK/p for p in [
    'cxx-screen/owner_host_tests.cpp','runtime-inventory/inventory.cpp',
    'runtime-stream/stream_probe.cpp','runtime-stream/stream_probe_v2.cpp',
    'runtime_sensitivity.h','runtime_module_sha.h']]
source_paths += [path for path in WORK.rglob('*.json') if path.name in {
    'actions.json','actions-generated.json','additional-frozen-inputs.json',
    'retention-gate-frozen.json','preregistered-control.json'}]
index = WORK/'preservation-git-index'
assert not index.exists()
env = dict(os.environ, GIT_INDEX_FILE=str(index))
query(['git','read-tree','87ae628f25ac99e8e4c412c4d21eb868753dce60'],env=env)
source_bindings = {}
seen = {}
for path in sorted(set(source_paths)):
    data = path.read_bytes(); digest = sha(data)
    if digest not in seen:
        target = 'attestations/software-vulkan/87ae628/portable-stage-seven-owner/source/'+str(path.relative_to(WORK))
        blob = query(['git','hash-object','-w','--stdin'],input=data).decode().strip()
        query(['git','update-index','--add','--cacheinfo','100644',blob,target],env=env)
        seen[digest] = target
    source_bindings[str(path.relative_to(ROOT))] = {'sha256':digest,'git_path':seen[digest]}
tree = query(['git','write-tree'],env=env).decode().strip()
commit = query(['git','commit-tree',tree,'-p','87ae628f25ac99e8e4c412c4d21eb868753dce60'],
    input=b'Preserve rejected portable stage-seven owner executable evidence\n\nFinite word replay and all eight numerical controls pass. Fixed latency gate passes12/20; both observed-instance accounting gates fail. Automatic renderer adoption is withdrawn separately. No frame or full qualification claim.\n').decode().strip()
tag = 'evidence/issue-44-portable-stage7-owner-87ae628'
query(['git','tag','-a',tag,commit,'-m','Immutable executable evidence for the rejected portable partial stage-seven owner'])
index.unlink()
(ARCHIVE/'source-git-bindings.json').write_text(json.dumps({'tag':tag,'commit':commit,'sources':source_bindings},indent=2)+'\n')

# Resolve every historical consumed origin/version, preserving exact old source
# epochs through Git and exact private outputs through the archive.
by_digest = {}
for relative, identity in manifest.items():
    by_digest.setdefault((identity['bytes'],identity['sha256']),[]).append(relative)
resolved = {}; unresolved = []; cache = {}
for ledger in sorted(WORK.rglob('inputs-before.json')):
    owner_path = ledger.parent/'owner.json'
    owner = json.loads(owner_path.read_text())
    revision = owner['source']['revision']
    for origin, identity in json.loads(ledger.read_text()).items():
        key = (origin,identity['bytes'],identity['sha256'])
        if key in resolved: continue
        wanted = (identity['bytes'],identity['sha256'])
        path = Path(origin) if Path(origin).is_absolute() else ROOT/origin
        recovery = None
        if path.is_file():
            data = path.read_bytes()
            if (len(data),sha(data)) == wanted:
                recovery = {'kind':'archive' if origin in manifest else 'retained', 'path':origin}
        if recovery is None and not Path(origin).is_absolute():
            git_key = (revision,origin)
            if git_key not in cache:
                result = subprocess.run(['git','show',revision+':'+origin],capture_output=True)
                cache[git_key] = result.stdout if result.returncode == 0 else None
            data = cache[git_key]
            if data is not None and (len(data),sha(data)) == wanted:
                recovery = {'kind':'git','revision':revision,'path':origin}
        if recovery is None and wanted in by_digest:
            recovery = {'kind':'archive','path':by_digest[wanted][0]}
        entry = {'origin':origin,'bytes':wanted[0],'sha256':wanted[1], 'recovery':recovery}
        if recovery is None: unresolved.append(entry)
        else: resolved[key] = entry
(ARCHIVE/'historical-input-recovery.json').write_text(json.dumps({'resolved':list(resolved.values()),'unresolved':unresolved},indent=2)+'\n')
print(json.dumps({'archive':str(ARCHIVE.relative_to(ROOT)),'files':len(manifest), 'bytes':sum(x['bytes'] for x in manifest.values()),
    'source_commit':commit,'tag':tag,'unique_source_files':len(seen),'resolved_origin_versions':len(resolved),'unresolved':unresolved},indent=2))
assert not unresolved, 'No cleanup until every consumed input version resolves exactly'
