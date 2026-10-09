"""Read-only static comparison of frozen original artifacts; no optimizer/device calls."""
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
sys.dont_write_bytecode = True
ROOT = Path.cwd()
OWN = Path(__file__).resolve().parent
BASE = ROOT / 'out/retained-qualification/cold-render-diagnostic/shader-ssa/o1-structure-review'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
binding = lambda p: {'bytes': p.stat().st_size, 'sha256': sha(p)}
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result

def main():
    revision = '8424fd72b16e599877a8057667cf07314b19752d'
    assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == revision
    assert not subprocess.check_output(['git', 'status', '--porcelain'], text=True)
    helpers = {'structure': BASE/'reproduce.py', 'constraints': BASE.parent/'scalar-replacement-review/review.py', 'surface': BASE/'transport-clean-compact-review.py'}
    before_helpers = {n: binding(p) for n,p in helpers.items()}
    stats = load('existing_structure', helpers['structure'])
    checks = load('existing_constraints', helpers['constraints'])
    semantics = load('existing_surface', helpers['surface'])
    oldpath = ROOT/'attestations/software-vulkan/125a890/ray-camera-first-and-repeat/verification.json'
    newpath = ROOT/'attestations/software-vulkan/8424fd7/ray-camera-first-and-repeat/verification.json'
    old = json.loads(oldpath.read_text())['all_recorded_artifacts_rehashed']
    new = json.loads(newpath.read_text())['all_recorded_artifacts_rehashed']
    modules = [p for p in new if '/retained/' in p and p.endswith('.spv')]
    assert len(modules) == 24
    for p in modules: assert binding(ROOT/p) == new[p]
    native = [p for p in modules if '_portable' not in p]
    assert len(native) == 12 and all(new[p] == old[p] for p in native)
    dis = Path('/usr/bin/spirv-dis')
    tool = binding(dis)
    prior = json.loads((ROOT/'out/retained-qualification/cold-render-diagnostic/leading-bit-normalization-review/current125-structure/report.json').read_text())
    assert tool == {k:prior['disassembler'][k] for k in ['bytes','sha256']}
    archive = ROOT/'attestations/software-vulkan/125a890/ci-integration-non-render/windows-native-build.zip'
    assert sha(archive) == 'eb20401b1e9ffa17a4f2a2f079b5bc1983974cf35559f1935cfef3a34c97a332'
    with zipfile.ZipFile(archive) as z:
        record = json.loads(z.read('windows-build.json'))
        backend = z.read('native-build-tested-sirius_backend_tests')
        assert len(backend) == record['artifacts']['native-build-tested-sirius_backend_tests']['bytes']
        assert hashlib.sha256(backend).hexdigest() == record['artifacts']['native-build-tested-sirius_backend_tests']['sha256']
    offsets = [m.start() for m in re.finditer(re.escape(b'\x03\x02\x23\x07'), backend)]
    report = {'scope':'Offline source/module structural comparison only; no compiler, optimizer, GPU, timing or allocation-cause verdict', 'baseline_revision':'125a8906639df0f09e461e6b0ee55bd76e366dae','current_revision':revision,'baseline_captured_verification':binding(oldpath),'current_captured_verification':binding(newpath),'disassembler':{'path':str(dis),**tool},'existing_helpers':{n:{'path':str(p),**before_helpers[n]} for n,p in helpers.items()},'baseline_provider_archive':binding(archive),'baseline_embedded_backend':record['artifacts']['native-build-tested-sirius_backend_tests'],'all12native_byte_identical':True,'all12portable_changed':all(old[p]!=new[p] for p in modules if '_portable' in p),'stages':[]}
    surface = lambda text: semantics.normalized_surface(text[:re.search(r'^\s*%\S+\s*=\s*OpFunction ',text,re.M).start()])
    with tempfile.TemporaryDirectory(prefix='readonly-spv-',dir=OWN) as temporary:
        tmp=Path(temporary)
        for stage in ['ray_camera','transport','dense']:
            name='bin/linux-gcc/src/sirius/backend/retained/retained_'+stage+'_portable.spv'
            size=old[name]['bytes'];matches=[p for p in offsets if hashlib.sha256(backend[p:p+size]).hexdigest()==old[name]['sha256']]
            assert matches, (stage,matches)  # Identical portable-mode arrays may each be embedded.
            previous=tmp/(stage+'-125.spv');previous.write_bytes(backend[matches[0]:matches[0]+size]);assert binding(previous)==old[name]
            current=ROOT/name
            structures=[]
            for label,module in [('baseline',previous),('current',current)]:
                assembly=tmp/(stage+'-'+label+'.spvasm')
                subprocess.run([str(dis),str(module),'-o',str(assembly)],check=True,timeout=30)
                text=assembly.read_text();checks.constraints(text);structure=stats.structure(assembly)
                structure['bytes']=module.stat().st_size
                structure['static_opcode_instructions']=sum(structure['opcode_histogram'].values())
                structures.append((structure,surface(text)))
            assert structures[0][1]==structures[1][1]
            assert old[name]==old[name.replace('_portable.spv','_portable_fp64.spv')]
            assert new[name]==new[name.replace('_portable.spv','_portable_fp64.spv')]
            report['stages'].append({'stage':stage,'baseline_module':old[name],'current_module':new[name],'baseline_bytes_recovered_from_exact_hashed_backend_offsets':matches,'baseline_structure':structures[0][0],'current_structure':structures[1][0],'integer_only_LocalSize1_no_barriers':True,'exact_normalized_interface_layout_globals':True,'both_portable_product_mode_modules_byte_identical':True})
    change=subprocess.check_output(['git','diff','125a890..8424fd7','--','src/sirius/kernels/portable_binary32.h'])
    (OWN/'source-change.patch').write_bytes(change)
    report['source_change_patch']=binding(OWN/'source-change.patch')
    report['source_change']='Only PB32ProductResidual changed: bounded finite nonzero branch computes one exact wide product and directly reuses it for the same high projection instead of PB32Multiply plus a second unpack/wide product; zero branch gives the same XOR sign and +0 residual. Domain/representation/rounding/caller sequence remains as separately reviewed.'
    report['limits']='Static instruction/size counts cannot identify driver compiler memory, optimization phase, runtime speed or benefit. Completed numerical/timing receipts are separate; remaining7 was active during this static review.'
    for p in modules: assert binding(ROOT/p)==new[p]
    assert binding(dis)==tool and {n:binding(p) for n,p in helpers.items()}==before_helpers
    assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==revision and not subprocess.check_output(['git','status','--porcelain'],text=True)
    report['bindings_unchanged_after']=True
    output=OWN/'report.json';assert not output.exists();output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'report':binding(output),'stages':[{'stage':r['stage'],'bytes_before':r['baseline_module']['bytes'],'bytes_after':r['current_module']['bytes'],'static_instructions_before':r['baseline_structure']['static_opcode_instructions'],'static_instructions_after':r['current_structure']['static_opcode_instructions'],'IMul_before':r['baseline_structure']['opcode_histogram'].get('OpIMul',0),'IMul_after':r['current_structure']['opcode_histogram'].get('OpIMul',0),'functions_before':r['baseline_structure']['functions'],'functions_after':r['current_structure']['functions'],'calls_before':r['baseline_structure']['function_calls'],'calls_after':r['current_structure']['function_calls'],'loops_before':r['baseline_structure']['loops'],'loops_after':r['current_structure']['loops']} for r in report['stages']]},indent=2))
if __name__=='__main__':main()
