"""Freeze original producer/input/provider and complete historical readback authority."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, subprocess, xml.etree.ElementTree as ET
root=Path.cwd();work=Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
source=json.loads((work/'source.json').read_text())
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==source['source_revision']
assert not subprocess.check_output(['git','status','--porcelain'])
archive=root/'attestations/native-vulkan/21e55ce/guarded-native-default-fma-rejected/manifest.json'
assert sha(archive)=='aea15fed97a6af8d384fb6b628bf6357c5b3c9403a1b50f3ee85ac5ccb70ad53'
manifest=json.loads(archive.read_text());refs={}
for ending in ('native/matched/a1/timestamps/gtest.xml','native/matched/a1/timestamps/driver-modules.json'):
    rows=[v for v in manifest['diagnostics'] if v['origin'].endswith('/'+ending)]
    assert len(rows)==1;v=rows[0];p=root/v['recovery']['path']
    assert p.stat().st_size==v['bytes'] and sha(p)==v['sha256']
    refs[ending]={'path':v['recovery']['path'],'bytes':v['bytes'],'sha256':v['sha256']}
old=ET.fromstring((root/refs['native/matched/a1/timestamps/gtest.xml']['path']).read_bytes())
assert old.get('tests')=='1' and old.get('failures')==old.get('errors')=='0'
case=old.findall('.//testcase');assert len(case)==1
assert case[0].get('classname')+'.'+case[0].get('name')=='RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults'
assert case[0].get('status')=='run' and case[0].get('result')=='completed'
props={}
for item in old.findall('.//property'):
    assert item.get('name') not in props;props[item.get('name')]=item.get('value')
readbacks={k:v for k,v in props.items() if '_readback_' in k};modules={k:v for k,v in props.items() if k.startswith('selected_module_')}
assert len(readbacks)==32 and len(modules)==7 and props['device_observations']=='51'
provider=json.loads((root/refs['native/matched/a1/timestamps/driver-modules.json']['path']).read_text(encoding='utf-8-sig'))
assert len(provider)==2 and {v['name'].lower() for v in provider}=={'vulkan-1.dll','amdvlk64.dll'}
staged=json.loads((work/'baseline-reconstruction.json').read_text());stage=root/staged['stage']
assert staged['whole_source_trees_exact'] and staged['native_producer_revision']==source['baseline_revision'] and staged['live_source_revision']==source['source_revision']
files={}
def seal(p):
    p=p.resolve();name=str(p.relative_to(root)) if p.is_relative_to(root) else str(p)
    b=p.read_bytes();value={'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
    assert name not in files or files[name]==value;files[name]=value
for p in stage.rglob('*'):
    if p.is_file(): seal(p)
for filename in ('tools-extraction.json','tools-gui-extraction.json'):
    inventory=json.loads((work/filename).read_text());seal(work/filename)
    for item in inventory['files']:
        p=root/item['path'];seal(p);assert files[item['path']]=={k:item[k] for k in ('bytes','sha256')}
for item in provider:
    p=Path('/mnt/'+item['path'][0].lower()+item['path'][2:].replace('\\','/'))
    seal(p);assert files[str(p.resolve())]=={k:item[k] for k in ('bytes','sha256')}
for item in refs.values():seal(root/item['path'])
for n in ('source.json','baseline-reconstruction.json','capture-range-source-review.json','service-info-loopback-result.json','run_original_control.ps1'):
    seal(work/n)
reference={'historical_producer':source['baseline_revision'],'historical_raw_A1':refs,'readbacks':readbacks,'selected_modules':modules,'provider_properties':{k:props[k] for k in ('selected_device','selected_driver','retained_scalar_route','product_mode','timestamp_valid_bits','timestamp_period_ns','coupled_rows','coupled_projection_budget','fma_fp32_enabled','preserves_fp32_signed_zero_inf_nan')},'ordered_marker_labels':[dict(part.split('=',1) for part in props['device_observation_'+str(i)].split(';')) for i in range(51)],'driver_modules':provider,'scope':'Original complete same-provider readbacks and module/stream authority; historical timings are not a performance comparison for a profiled run.'}
reference_path=work/'capture-reference.json';assert not reference_path.exists();reference_path.write_text(json.dumps(reference,indent=2)+'\n');seal(reference_path)
plan={'source_revision':source['source_revision'],'native_producer_revision':source['baseline_revision'],'stage':staged['stage'],'gate_sha256':staged['gate_sha256'],'tools':'bin/windows-msvc/qualification-tools/rdts-2026-10-02-1872/tools','filter':'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults','dispatch_spec':'dispatch:43:2','expected_dispatch_dimensions':[[12,1,1],[24,1,1]],'original_control_seconds':90,'original_control_rss_bytes':4294967296,'tool_rss_bytes':1073741824,'bridge_seconds':300,'clock_or_driver_experiments':[],'shader_instrumentation':False,'files':files,'scope':'One actual instruction/counter trace attempt of the unchanged original default timestamp control. Trace must independently confirm dispatch identity/dimensions and sample availability. No benefit/frame/full-science/release qualification.'}
p=work/'capture-plan.json';assert not p.exists();p.write_text(json.dumps(plan,indent=2)+'\n')
print(json.dumps({'frozen_files':len(files),'complete_readbacks':32,'ordered_markers':51,'original_selected_modules':7,'provider_modules_exact':True}))
