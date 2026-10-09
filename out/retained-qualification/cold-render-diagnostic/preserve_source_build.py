"""Capture exact strict source/build authority and portable module deltas."""
import hashlib,json,os,pathlib,shutil,subprocess,sys

root=pathlib.Path.cwd().resolve()
folder=pathlib.Path(__file__).resolve().parent
revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()
short=revision[:7]
assert not subprocess.check_output(['git','status','--porcelain'],cwd=root,text=True)
assert len(sys.argv) in [2,3] and pathlib.Path(sys.argv[1]).name==sys.argv[1]
stages=['camera','ray_camera']
if len(sys.argv)==3:
    assert sys.argv[2]=='--all-portable'
    stages=['camera','dense','endpoint','initialize','ray_camera','transport']
baseline_path=root/'out/retained-qualification/portable-integration'/sys.argv[1]/'identity.json'
baseline=json.loads(baseline_path.read_text())
assert baseline['status']=='terminal' and baseline['passed']
build=root/'bin/linux-gcc'
generated=build/'generated/sirius'
gate=json.loads((generated/'native_build_gate.json').read_text())
assert gate['source']=={'revision':revision,'clean':True}
assert gate['status']=='passed' and gate['ctest']['executed']==9
assert all(gate['ctest'][k]==0 for k in ['failures','errors','skipped'])
reader=folder/f'authority-reader-{short}.log'
with reader.open('wb') as output:
    subprocess.run([sys.executable,str(root/'scripts/verify-build-gate.py'),
      '--action','check-native-build','--stamp',str(generated/'native_build_gate.json'),
      '--alignment-mode','qualification','--source-root',str(root),'--build-dir',str(build),
      '--source-revision',revision,'--ctest','/usr/bin/ctest','--config','Release'],
      cwd=root,stdout=output,stderr=subprocess.STDOUT,check=True,timeout=60)
runtime=folder/f'ray-camera-reuse{short}-runtime-identity.json'
environment=os.environ.copy()
environment.update(VK_DRIVER_FILES='/usr/share/vulkan/icd.d/lvp_icd.json',
                   SIRIUS_VULKAN_DEVICE='0',MESA_SHADER_CACHE_DISABLE='true')
subprocess.run([sys.executable,str(root/'scripts/runtime_identity.py'),
  '--source-root',str(root),'--source-revision',revision,'--source-tree-clean','true',
  '--executable',str(build/'src/sirius/app/sirius'),'--require-vulkan','true',
  '--output',str(runtime)],cwd=root,env=environment,check=True,timeout=60)
sys.path.insert(0,str(root/'scripts'))
from runtime_identity import validate_record,identity
validate_record(json.loads(runtime.read_text()),revision,
                identity(build/'src/sirius/app/sirius'),required=True,clean=True)
old={pathlib.Path(k).name:r['sha256'] for k,r in baseline['artifacts'].items() if k.endswith('.spv')}
modules=sorted((build/'src/sirius/backend/retained').glob('retained_*.spv'))
assert len(modules)==len(old)==24
current={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in modules}
changed=sorted(name for name,value in current.items() if value!=old[name])
assert set(changed)=={f'retained_{stage}_{mode}.spv' for stage in stages
                     for mode in ['portable','portable_fp64']}
for stage in stages:
    assert current[f'retained_{stage}_portable.spv']==current[f'retained_{stage}_portable_fp64.spv']
comparison={'source_revision':revision,'baseline_revision':baseline['source']['revision'],
 'baseline_identity_sha256':hashlib.sha256(baseline_path.read_bytes()).hexdigest(),
 'changed':changed,'unchanged':sorted(set(current)-set(changed)),
 'all24_module_sha256':current,'two_public_portable_camera_mode_bytes_equal':True,
 'portable_mode_byte_equivalence_stages':stages}
comparison_path=folder/f'modules-{short}.json'
comparison_path.write_text(json.dumps(comparison,indent=2)+'\n')
destination=root/'attestations/software-vulkan'/short/'source-build-authority'
destination.mkdir(parents=True,exist_ok=False)
files={'configure.log':folder/f'configure-{short}.log',
 'build-and-authority.log':folder/f'build-and-authority-{short}.log',
 'native_build_gate.json':generated/'native_build_gate.json',
 'native_build_gate_junit.xml':generated/'native_build_gate_junit.xml',
 'native_build_gate_ctest.log':generated/'native_build_gate_ctest.log',
 'alignment_receipt.json':generated/'alignment_receipt.json',
 'authority-reader.log':reader,'runtime-identity.json':runtime,
 'module-comparison.json':comparison_path,'retained_pair.slang':root/'src/sirius/kernels/retained_pair.slang',
 'build-retained-kernels.py':root/'scripts/build-retained-kernels.py'}
for extra in ['configure-optimizer-check', 'configure-missing-optimizer-verified', 'optimizer-build-identity']:
    suffix = 'log' if extra == 'configure-missing-optimizer-verified' else 'json'
    evidence = folder / f'{extra}-{short}.{suffix}'
    if evidence.is_file():
        files[f'{extra}.{suffix}'] = evidence
if len(sys.argv)==3:
    for stage in ['triple','endpoint','transport','dense','initialize']:
        files[f'retained_{stage}.slang']=root/f'src/sirius/kernels/retained_{stage}.slang'
for name,path in files.items():
    shutil.copyfile(path,destination/name)
    assert (destination/name).read_bytes()==path.read_bytes()
classification={'scope':'Strict source/build and nine non-render authority controls only; no full Mandatory/scientific/native-runtime/release qualification',
 'source_revision':revision,'original_topology':{k:len(gate[k]) for k in ['tested_artifacts','product_artifacts','test_input_artifacts']},
 'ctest':gate['ctest'],'actual_runtime_identity_reader':'PASS',
 'receipt_and_log_bytes_preserved_before_further_build':True,'all24_module_deltas_checked':True,'qualification_claimed':False}
(destination/'classification.json').write_text(json.dumps(classification,indent=2)+'\n')
manifest={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
          for p in sorted(destination.iterdir()) if p.is_file()}
(destination/'sha256-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()==revision
assert not subprocess.check_output(['git','status','--porcelain'],cwd=root,text=True)
print(json.dumps({'bundle':str(destination.relative_to(root)),'payload_count':len(manifest),
 'manifest_sha256':hashlib.sha256((destination/'sha256-manifest.json').read_bytes()).hexdigest(),
 'changed_modules':changed,'unchanged_modules':len(current)-len(changed)}))
