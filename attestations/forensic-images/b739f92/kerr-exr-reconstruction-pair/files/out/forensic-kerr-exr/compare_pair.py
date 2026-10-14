from pathlib import Path
import hashlib,json,subprocess,time
root=Path.cwd();work=Path(__file__).resolve().parent
m=json.loads((work/'preparation.json').read_text());actions=json.loads((work/'actions.json').read_text());results={}
for name,action in actions.items():
 r=json.loads((work/name/'result.json').read_text());assert r['classification']=='COMPLETE_ORIGINAL_SIZE_SINGLE_RENDER' and r['complete_original_size_render'] and r['self_decode']['passed'] and not r['proof_gaps'];results[action['role']]=r
assert set(results)=={'legacy','ported'}
left=Path(results['legacy']['output']['path']);right=Path(results['ported']['output']['path']);tool=work/'image_compare'
def identity(p):
 b=p.read_bytes();return {'path':str(p),'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
inputs=[identity(p) for p in (left,right,tool)]
for n,role in enumerate(('legacy','ported')):assert {k:inputs[n][k] for k in ('bytes','sha256')}=={k:results[role]['output'][k] for k in ('bytes','sha256')}
assert {k:inputs[2][k] for k in ('bytes','sha256')}==results['legacy']['self_decode']['binary']==results['ported']['self_decode']['binary']
argv=[str(tool),'compare','exr','512','512',str(left),str(right)];start=time.monotonic();r=subprocess.run(argv,cwd=root,env={'PATH':'/usr/bin:/bin','LC_ALL':'C','TZ':'UTC'},capture_output=True,text=True,timeout=30);assert r.returncode in (0,1) and not r.stderr;assert inputs==[identity(p) for p in (left,right,tool)]
d=json.loads(r.stdout);assert d['compared_pixels']==262144 and d['compared_samples']==1048576 and d['sample_domain']=='decoded-RGBA-IEEE754-binary32-words-from-HALF' and d['equal']==(r.returncode==0)
result={'evidence_kind':'NEW_RECONSTRUCTION','case':'kerr-exr','source_export_revision':m['source_export_revision'],'legacy_runtime_origin':m['legacy_runtime_origin'],'runner_revision':m['runner_revision'],'source_scaffold_overlay':'one115-byte GLFW CMake config template; runtime/asset bytes unchanged','requested_spp':32,'effective_spp':25,'rays_per_renderer':6553600,'argv':argv,'inputs':inputs,'render_results':{k:{'owner_wall_seconds':v['ownership']['wall_seconds'],'owner_peak_sampled_rss_kib':v['ownership']['peak_sampled_rss_kib'],'workers':v['runtime_join']['workers']} for k,v in results.items()},'comparison':d,'comparison_wall_seconds':time.monotonic()-start,'historical_tape_identity_claimed':False,'current_physics_or_scientific_or_release_qualification_claimed':False,'reconstruction_pair_complete':True}
(work/'pair-result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
