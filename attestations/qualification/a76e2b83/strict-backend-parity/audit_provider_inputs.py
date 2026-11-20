"""Verify recovery inputs for the successful native-owned provider without extracting them."""

import hashlib
import json
from pathlib import Path
import subprocess
from zipfile import ZipFile

ROOT=Path.cwd()
QUALIFICATION=ROOT/'attestations/qualification'
CONTRACT=ROOT/'out/native-cpu-closed-all-provider/provider/source-contract.json'
SIMD16=QUALIFICATION/'2ec9eea2/native-cpu-simd16-original-gates/MANIFEST.json'
MODULES=QUALIFICATION/'bfa095ba/native-cpu-row-ownership-finite-outcome/new-actual-products-01.zip'
MESA_REVISION='43389d91ca54667d12a41ee0a33d35e596fc42bd'
MESA_MEMBER='attestations/qualification/8c2be67d/lavapipe-vector-workgroup-rejected/mesa-25.2.8.tar.xz'
MESA_PIN=(43813260,'097842f3e49d996868b38688db87b006f7d4541e93ce86d2f341d8b3e7be7c93')
OBJECT_ARCHIVES=(
    ('a47e219b/native-cpu-dispatch-admission-token-finite-compatibility/component.zip', range(3), ''),
    ('5ae04321/native-cpu-endpoint-header-frame-timeout/actual-products.zip', (3,), ''),
    ('86843614/native-cpu-local-identities-finite-compatibility/actual-products.zip', range(4,7), 'files/'),
    ('363aa0c2/native-cpu-closed-all-original-c8/products.zip', (7,), ''),
)


def identity(stream):
    size=0
    digest=hashlib.sha256()
    while chunk:=stream.read(1024*1024):
        size+=len(chunk)
        digest.update(chunk)
    return size,digest.hexdigest()


def verify_member(archive,member,expected):
    with archive.open(member) as stream:
        found=identity(stream)
    assert found==(int(expected['bytes']),expected['sha256']), member
    return {'member':member,'bytes':found[0],'sha256':found[1]}


def version(*argv):
    return subprocess.check_output(argv,text=True).strip()


def main():
    with subprocess.Popen(['git','show',MESA_REVISION+':'+MESA_MEMBER],stdout=subprocess.PIPE) as process:
        assert identity(process.stdout)==MESA_PIN
    assert process.returncode==0
    contract=json.loads(CONTRACT.read_text())
    objects=[]
    for path,indexes,prefix in OBJECT_ARCHIVES:
        with ZipFile(QUALIFICATION/path) as archive:
            for index in indexes:
                original=contract['native_objects'][index]
                record=verify_member(archive,prefix+original['path'],original)
                record['archive']=str((QUALIFICATION/path).relative_to(ROOT))
                objects.append(record)
    assert len(objects)==8
    selected=json.loads(SIMD16.read_text())['original_c8']['selected_modules']
    modules={}
    with ZipFile(MODULES) as archive:
        for name in ('transport','endpoint','dopri_phase'):
            original=selected[name]
            modules[name]=verify_member(archive,'sha256/'+original['sha256'],original)
    tools={'gcc':version('g++-14','-dumpfullversion'),
           'llvm':version('/usr/lib/llvm-21/bin/llvm-config','--version'),
           'installed_meson':version('meson','--version'),
           'sealed_build_meson':'1.9.1'}
    report={'scope':'Read-only verification of the recoverable archive inputs; no provider rebuild or runtime admission',
            'upstream_mesa':{'git_revision':MESA_REVISION,'member':MESA_MEMBER,'bytes':MESA_PIN[0],
                             'sha256':MESA_PIN[1]},
            'native_object_members':objects,'native_module_members':modules,
            'native_module_archive':str(MODULES.relative_to(ROOT)),'tool_versions':tools}
    path=ROOT/'out/strict-backend-science/provider-input-audit.json'
    path.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print('Verified the Mesa source, eight pinned native objects and three native-owned modules')


if __name__=='__main__':
    main()
