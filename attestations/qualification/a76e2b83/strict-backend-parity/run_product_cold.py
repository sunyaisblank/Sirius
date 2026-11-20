"""Replay the registered cold Kerr command with the pinned private provider."""

import json
import os
import re
from pathlib import Path
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

ROOT = Path.cwd()
sys.path.insert(0, str(ROOT/'out/science-replay'))
from run_science import GUARD, MODEL, elf_dependencies, pin, read, require, write

BUILD = ROOT/'bin/linux-ci'
TASK = ROOT/'out/strict-backend-science'
WORK = TASK/'private/original-cold'
RENDERS = ROOT/'renders/strict-backend-science/private/original-cold'
NATIVE = TASK/'native-provider'
ADMISSION = ROOT/'out/simd16-provider-admission'
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
    candidate=NATIVE/'matched-runtime'
    executable=candidate/'sirius_render_tests'
    link=read(candidate/'link.json')
    require(link['revision']==revision and pin(executable)==link['executable'],
            'The privately linked executable must bind the current original inputs')
    argv[1]=f'-DSIRIUS_RENDER_TESTS={executable}'
    authority=read(AUTHORITY)
    prior=authority['original_gates']['cold_frame']
    provider=ROOT/prior['provider']['path']
    icd=ADMISSION/'lvp_icd.json'
    c8=read(ADMISSION/'c8-acceptance.json')
    require(authority['adoption_chain_passed'] and prior['accepted'] and
            pin(provider)==prior['provider']==c8['provider'] and c8['source']==revision and
            len(c8['record']['readbacks'])==32 and
            read(icd)['ICD']['library_path']==str(provider),
            'Require the original provider, historical controls and current C8')
    require('SIRIUS_CERTIFIED_SPHERICAL_ESCAPE' not in os.environ,
            'The private strict candidate must use its compiled continuation policy')
    extra=[candidate/'link.json',NATIVE/'selector/owned-frame-modules.json',
           ADMISSION/'linked-provider.json',ADMISSION/'c8-acceptance.json',
           ADMISSION/'c8-owner/control-8/owner.json',TASK/'certified-link/link.json']
    limits=prior['action']['limits']
    require(limits=={'reserve_kib':2097152,'rss_kib':4194304,'seconds':300},
            'Original cold-frame governors differ')
    resources=executable.parent/'resources'
    require((resources/'assets/Starfield.png').samefile(ROOT/'assets/Starfield.png'),
            'Build must use the existing starfield inode')
    frozen=[executable,runtime/'sirius',provider,icd,PROFILE,TEST_SOURCE,AUTHORITY,SCENE_AUTHORITY,
            MODEL,receipt,WORK/'registered.json',Path(__file__).resolve(),Path('/usr/bin/cmake')]
    frozen+=sorted(path for path in resources.rglob('*') if path.is_file())
    frozen+=sorted((BUILD/'src/sirius/backend/retained').glob('*.spv'))
    frozen+=elf_dependencies({executable,provider})+extra
    WORK.mkdir(parents=True)
    RENDERS.mkdir(parents=True)
    write(WORK/'registered.json',registered)
    shutil.copyfile(GUARD/'run_owned.py',WORK/'run_owned.py')
    (WORK/'owned_guard.py').symlink_to(GUARD/'owned_guard.py')
    action={'revision':revision,'kind':'test','expected_exit':0,'argv':argv,'limits':limits,
        'software_provider':True,'disable_software_shader_cache':True,'runtime_icd':str(icd),
        'runtime_provider':str(provider),'frozen_inputs':list(map(str,dict.fromkeys(frozen))),
        'build_environment':{'TMPDIR':str(RENDERS)},
        'scope':'Private exact provider and strict relink diagnostic under concurrent host tests'}
    write(WORK/'actions.json',{'original-cold':action})
    finished=subprocess.run([sys.executable,'-B',str(WORK/'run_owned.py'),'original-cold'])
    capture=WORK/'original-cold'
    owner=read(capture/'owner.json')
    require(owner['status']=='terminal' and owner['source_unchanged']
        and owner['whole_inputs_unchanged'] and not owner['remaining_owned_processes']
        and not owner['cleanup_errors'] and owner['provider_mappings'],
        'Cold frame ownership, provider mapping or input verification failed')
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
    certified=[line for line in details if line.startswith('[Vulkan] Certified spherical escapes: ')]
    if confirmed:
        fields,=certified
        events=re.fullmatch(r'\[Vulkan\] Certified spherical escapes: ([1-9]\d*), '
                            r'exterior turns: (\d+), quadrature panels: ([1-9]\d*)',fields)
        require(events and int(events[2])<=int(events[1]) and
                render[0]['work_items']==16384 and render[0]['precision']=='fp32',
                'Require actual retained and certified work for the original private frame')
    result={'revision':revision,'accepted':bool(confirmed),'registered_command':registered,'exploratory':True,
        'profile':pin(PROFILE),'model':pin(MODEL),'executable':pin(executable),'alignment':pin(receipt),
        'provider':pin(provider),'icd':pin(icd),'action':pin(WORK/'actions.json'),'input_snapshot':pin(capture/'inputs-before.json'),
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
