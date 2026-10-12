from pathlib import Path
import hashlib,json,shutil
root=Path.cwd();work=Path(__file__).resolve().parent
bundle=root/'attestations/software-vulkan/dcd1368/sealed-stage-seven-owner'
bundle.mkdir(parents=True,exist_ok=False)
for source in [work,root/'bin/linux-gcc/sealed-stage-seven-owner']:
    destination=bundle/'files'/source.relative_to(root)
    shutil.copytree(source,destination)
identity=lambda f:{'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
manifest={str(f.relative_to(bundle)):identity(f) for f in sorted((bundle/'files').rglob('*')) if f.is_file()}
(bundle/'whole-byte-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
for name,expected in manifest.items():
    f=bundle/name
    assert identity(f)==expected
    original=root/str(Path(name).relative_to('files'))
    assert identity(original)==expected
print(json.dumps({'bundle':str(bundle.relative_to(root)),'files':len(manifest),'bytes':sum(v['bytes'] for v in manifest.values()),'manifest_sha256':identity(bundle/'whole-byte-manifest.json')['sha256'],'scope':'Raw experimental source/build/runtime evidence only. Failed original20+2 gates; no release/renderer adoption.'},indent=2))
