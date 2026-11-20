"""Replay the registered cold Kerr command with the ordinary qualification products."""

import json
import re
from pathlib import Path
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

ROOT = Path.cwd()
sys.path.insert(0, str(ROOT/'out/science-replay'))
from run_science import GUARD, ICD, MODEL, PROVIDER, elf_dependencies, pin, read, require, write

BUILD = ROOT/'bin/linux-ci'
WORK = ROOT/'out/strict-backend-science'/sys.argv[1]/'original-cold'
RENDERS = ROOT/'renders/strict-backend-science'/sys.argv[1]/'original-cold'
EXPLORATORY = sys.argv[1] == 'certified'
REGISTERED = 'OperationalVulkan.RequiredRuntimeDispatch'
DIRECT = 'RenderCommandParse.ExplicitGpuRequestRunsVulkanWhenDevicePresent'
AUTHORITY = ROOT/'attestations/qualification/34e40355/native-cpu-certified-original-gates/MANIFEST.json'
SCENE_AUTHORITY = ROOT/'attestations/qualification/2ec9eea2/native-cpu-simd16-original-gates/MANIFEST.json'
PROFILE = ROOT/'tests/operational/vulkan_runtime_profile.cmake'
TEST_SOURCE = ROOT/'tests/render/app_render_execution_test.cpp'


def run():
    require(shutil.disk_usage('/').free>=100_000_000_000, 'Less than 100 GB free on /')
    revision = subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    require(not subprocess.check_output(['git','status','--porcelain']), 'Cold frame must start clean')
    receipt = BUILD/'generated/sirius/alignment_receipt.json'
    alignment = read(receipt)
    require(alignment['source_revision']==revision and alignment['source_tree_clean']
            and alignment['alignment_mode']=='qualification','Use the exact strict build')
    for path in (PROFILE,TEST_SOURCE):
        original=subprocess.check_output(['git','show',f'2ec9eea2:{path.relative_to(ROOT)}'])
        require(path.read_bytes()==original,f'Original runtime workload differs at {path}')
    ctest=['/usr/bin/ctest','--test-dir',str(BUILD),'-C','Release','--show-only=json-v1','-R',
           '^'+re.escape(REGISTERED)+'$']
    listing = json.loads(subprocess.check_output(ctest))
    registered,=listing['tests']
    props={row['name']:row['value'] for row in registered['properties']}
    require(registered['name']==REGISTERED and props['TIMEOUT']==300 and
            'Mandatory' in props['LABELS'] and props['RUN_SERIAL'], 'Original registration changed')
    runtime=BUILD/'src/sirius/app'
    executable=runtime/'sirius_render_tests'
    argv=['/usr/bin/cmake',f'-DSIRIUS_RENDER_TESTS={executable}','-P',str(PROFILE)]
    require(registered['command']==argv, 'Cold command differs from the registered one')
    extra=[]
    if EXPLORATORY:
        candidate=WORK.parent.parent/'certified-link'
        executable=candidate/'sirius_render_tests'
        link=read(candidate/'link.json')
        require(link['checkout']==revision and pin(executable)==link['executable'],
                'The exploratory executable must bind the current original inputs')
        argv[1]=f'-DSIRIUS_RENDER_TESTS={executable}'
        extra=[candidate/name for name in ('link.json','vulkan_renderer.cpp','vulkan_renderer.o')]
    prior=read(AUTHORITY)['original_gates']['cold_frame']
    limits=prior['action']['limits']
    require(limits=={'reserve_kib':2097152,'rss_kib':4194304,'seconds':300},
            'Original cold-frame governors differ')
    resources=executable.parent/'resources'
    require((resources/'assets/Starfield.png').samefile(ROOT/'assets/Starfield.png'),
            'Build must use the existing starfield inode')
    frozen=[executable,runtime/'sirius',PROVIDER,ICD,PROFILE,TEST_SOURCE,AUTHORITY,SCENE_AUTHORITY,
            MODEL,receipt,WORK/'registered.json',Path(__file__).resolve(),Path('/usr/bin/cmake')]
    frozen+=sorted(path for path in resources.rglob('*') if path.is_file())
    frozen+=sorted((BUILD/'src/sirius/backend/retained').glob('*.spv'))
    frozen+=elf_dependencies({executable,PROVIDER})+extra
    WORK.mkdir(parents=True)
    RENDERS.mkdir(parents=True)
    write(WORK/'registered.json',registered)
    shutil.copyfile(GUARD/'run_owned.py',WORK/'run_owned.py')
    (WORK/'owned_guard.py').symlink_to(GUARD/'owned_guard.py')
    action={'revision':revision,'kind':'test','expected_exit':0,'argv':argv,'limits':limits,
        'software_provider':True,'disable_software_shader_cache':True,'runtime_icd':str(ICD),
        'runtime_provider':str(PROVIDER),'frozen_inputs':list(map(str,dict.fromkeys(frozen))),
        'build_environment':{'TMPDIR':str(RENDERS)},
        'scope':('Exploratory stock provider and certified continuation only' if EXPLORATORY else
            'Original cold 128x128 renderer on the ordinary stock Vulkan product')}
    write(WORK/'actions.json',{'original-cold':action})
    finished=subprocess.run([sys.executable,'-B',str(WORK/'run_owned.py'),'original-cold'])
    capture=WORK/'original-cold'
    owner=read(capture/'owner.json')
    require(owner['status']=='terminal' and owner['source_unchanged']
        and owner['whole_inputs_unchanged'] and not owner['remaining_owned_processes']
        and not owner['cleanup_errors'], 'Cold frame ownership or input verification failed')
    out=capture/'stdout.log'
    err=capture/'stderr.log'
    text=out.read_text()+'\n'+err.read_text()
    expected=read(SCENE_AUTHORITY)['original_cold_frame']['scene']
    scenes=[json.loads(line.split('Scene evidence: ',1)[1]) for line in text.splitlines()
            if line.startswith('[Session] Scene evidence: ')]
    require(not any(record!=expected for record in scenes), 'Observed cold scene differs')
    render=[json.loads(line.split('Render evidence: ',1)[1]) for line in text.splitlines()
            if line.startswith('[Vulkan] Render evidence: ')]
    xml=capture/'gtest.xml'
    cases=[]
    if xml.exists():
        cases=list(ET.parse(xml).iter('testcase'))
    confirmed= (not finished.returncode and owner['passed'] and owner['wall_seconds']<=300
                and len(cases)==1 and cases[0].get('classname')+'.'+cases[0].get('name')==DIRECT
                and all(cases[0].find(tag) is None for tag in ('failure','skipped','error'))
                and len(render)==1 and bool(scenes))
    if confirmed:
        require(render[0]['route']=='retained' and render[0]['accepted_intervals']>0
                and (render[0]['width'],render[0]['height'])==(128,128),
                'A different route cannot qualify the cold frame')
    progress=[line for line in text.splitlines() if line.startswith('[Vulkan] Retained progress: ')]
    details=[line for line in text.splitlines() if line.startswith((
        '[Vulkan] initialising','[Vulkan] software retained','[Vulkan] Certified spherical escapes:'))]
    result={'revision':revision,'accepted':bool(confirmed),'registered_command':registered,'exploratory':EXPLORATORY,
        'profile':pin(PROFILE),'model':pin(MODEL),'executable':pin(executable),'alignment':pin(receipt),
        'provider':pin(PROVIDER),'icd':pin(ICD),'action':pin(WORK/'actions.json'),'input_snapshot':pin(capture/'inputs-before.json'),
        'owner':{key:owner[key] for key in ('status','stop_reason','child_exit','wall_seconds','peak_sampled_rss_kib',
                    'provider_mappings','source_unchanged','whole_inputs_unchanged','remaining_owned_processes')},
        'stdout':pin(out),'stderr':pin(err),'xml':pin(xml) if xml.exists() else None,
        'scene':scenes[0] if scenes else None,'render_evidence':render,'preparation_and_route':details,
        'last_progress':progress[-3:],'diagnostic_tail':text.splitlines()[-45:]}
    write(WORK/'result.json',result)
    print(json.dumps({key:result[key] for key in ('accepted','owner','last_progress')},indent=2))
    return confirmed


if __name__=='__main__':
    sys.exit(0 if run() else 1)
