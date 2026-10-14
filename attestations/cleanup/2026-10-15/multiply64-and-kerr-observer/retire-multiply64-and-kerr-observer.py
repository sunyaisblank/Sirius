from pathlib import Path
import importlib.util,json,os,shutil,subprocess
root=Path.cwd().resolve()
helper=root/'attestations/cleanup/2026-10-15/terminal-isolated-parity-preflight-hints/cleanup-terminal-parity-preflight-hints.py'
spec=importlib.util.spec_from_file_location('fixed_cleanup',helper);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
m.ROOT=root;m.SOURCE='46a9599a0cb30b0667e09838eac7a87690e982e4'
rows=[
('evidence/issue-44-multiply64-original-controls','bf23b994d658e2a70e9931d6f9a4c1e239044973','attestations/software-vulkan/4bcd987/multiply64-original-controls','fc2f642794b6d17be5b614d25d0d5387fd4e5946eaa67b5b715c998bf4544bba',('out/portable-multiply64-controls-o0-lines','bin/linux-gcc/portable-multiply64-controls-o0-lines')),
('evidence/issue-44-multiply64-original-benefit','27cbda645f4602dc008b67ccd296ffb35b9cf59e','attestations/software-vulkan/4bcd987/multiply64-original-benefit','465d81593a707f5835274e8b831f22b6ffe39c4bda5c4c6722e491c28a5262b3',('out/portable-multiply64-benefit','bin/linux-gcc/portable-multiply64-benefit')),
('evidence/issue-53-p1-kerr-minimum-step-observer-rss-stop','11d17cc42761f5c40e925d7ec414178e8697991e','attestations/scientific-core/4bcd987/p1-kerr-minimum-step-observer-rss-stop','a13e169bbc619866fb8344487a7f27845a1eaed21a124c019873ea9757e9125f',('out/p1-kerr-minimum-step-observer','bin/linux-gcc/p1-kerr-minimum-step-observer'))]
m.TARGETS=tuple(root/name for row in rows for name in row[4]);m.prefix_target=lambda path,targets=None:any(path==target or path.is_relative_to(target) for target in (m.TARGETS if targets is None else targets))
m.BUNDLE=root/'attestations/cleanup/2026-10-15/multiply64-and-kerr-observer';m.TAG='evidence/issue-44-multiply64-and-kerr-observer-cleanup'
m.clean_source();assert not m.BUNDLE.exists()
assert subprocess.run(['git','show-ref','--verify','--quiet','refs/tags/'+m.TAG]).returncode==1
preserved=[m.preserved(row) for row in rows]
# Pending source changes consume canonical source and immutable archives only.
freeze_path=root/'out/minimum-step-floor-repair/source-freeze.json'
freeze=json.loads(freeze_path.read_text());dependencies=[]
for name,pin in freeze['pins'].items():
 path=root/name;assert m.identity(path)==pin
 assert not m.prefix_target(path) and not m.prefix_target(path.resolve())
 dependencies.append(name)
contracts=json.loads((root/'out/minimum-step-floor-repair/source-contracts.json').read_text())
lineage=contracts['finding_lineage'];finding=root/lineage['source_finding_path'];assert m.identity(finding)==lineage['source_finding']
assert not m.prefix_target(finding) and m.git('show',lineage['commit']+':'+str(finding.relative_to(root)))==finding.read_bytes()
for relative in ['out/macos-failure-diagnostics/ci.yml','out/macos-failure-diagnostics/changes.patch']:
 path=root/relative;assert path.is_file() and not m.prefix_target(path)
 dependencies.append(relative)
metadata=[]
for name in ['bin/linux-gcc/build.ninja','bin/linux-gcc/compile_commands.json']:
 path=root/name;raw=path.read_bytes()
 for target in m.TARGETS:
  # CMake-generated build source paths cannot consume any private experiment.
  assert os.fsencode(target) not in raw and os.fsencode(str(target.relative_to(root))) not in raw,name
 metadata.append({'path':name,**m.identity(path)})
owners=m.retired_owners();census=m.consumer_census();m.clean_source();final=m.consumer_census()
m.BUNDLE.mkdir(parents=True)
shutil.copy2(Path(__file__),m.BUNDLE/Path(__file__).name)
receipt={'source':m.SOURCE,'status':'verified-not-removed','preservation':preserved,'terminal_owners':owners,'future_source_inputs_verified':dependencies,'future_build_metadata':metadata,'future_observer_lineage_uses_preserved_archive':True,'initial_consumer_census':census,'final_consumer_census':final,'removed':[],'qualification_claimed':False,'global_consumer_absence_claimed':False}
m.write(m.BUNDLE/'receipt.json',receipt)
for target in m.TARGETS:
 assert target.is_dir() and not target.is_symlink() and target.resolve()==target
 shutil.rmtree(target);assert not target.exists()
 receipt['removed'].append(str(target.relative_to(root)));m.write(m.BUNDLE/'receipt.json',receipt)
m.clean_source();receipt['status']='preserved-terminal-roots-removed';m.write(m.BUNDLE/'receipt.json',receipt)
m.write(m.BUNDLE/'manifest.json',{'files':{path.name:m.identity(path) for path in sorted(m.BUNDLE.iterdir()) if path.is_file()},'qualification_claimed':False})
record=m.record_tag();print(json.dumps({**record,'removed':receipt['removed'],'terminal_owners_rechecked':len(owners),'accessible_local_consumers':len(census['accessible_local_consumers']),'permission_gaps':len(census['permission_gaps'])}))
