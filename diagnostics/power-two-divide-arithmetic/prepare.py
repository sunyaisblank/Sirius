"""Prepare exact-corpus diagnostics after the independently reviewed compiler screen."""
from pathlib import Path
import hashlib
import importlib.util
import json
import random
import re
import shlex
import shutil
import struct
import subprocess
import sys
sys.dont_write_bytecode = True
R = Path.cwd().resolve()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc'
T = B / 'power-two-divide-arithmetic'
S = R / 'out/power-two-divide-screen'
REV = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == REV
assert not subprocess.check_output(['git', 'status', '--porcelain=v1'])
assert not T.exists()
T.mkdir()
def record(path):
    with Path(path).open('rb') as f:
        digest = hashlib.file_digest(f, 'sha256').hexdigest()
    return {'bytes': Path(path).stat().st_size, 'sha256': digest}
def write(name, value):
    (W / name).write_text(json.dumps(value, indent=2) + '\n')
for name in ['run_owned.py', 'owned_guard.py']:
    shutil.copyfile(S / name, W / name)
reference_path = R / 'tests/support/portable_binary32/reference.py'
spec = importlib.util.spec_from_file_location('issue99_exact_reference', reference_path)
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)
ref.selfcheck()
fixture_path = B / 'tests/backend/portable_binary32_reference.bin'
raw = fixture_path.read_bytes()
assert record(fixture_path)['sha256'] == '3fb1433f9cd14fbc9fd2b99a38a4893a4a571150df10f2fbe58172318fc7691e'
assert struct.unpack_from('<III', raw) == (ref.MAGIC, ref.VERSION, ref.RECORD_COUNT)
assert len(raw) == 12 + 24 * ref.RECORD_COUNT
records = list(struct.iter_unpack('<6I', raw[12:]))
numerators = [0,1,2,3,4,5,7,8,0x3fffff,0x400000,0x400001,0x7ffffe,0x7fffff,
    0x800000,0x800001,0xffffff,0x1000000,0x3f000000,0x3f000001,0x3f7fffff,
    0x3f800000,0x3f800001,0x3fc00000,0x3fffffff,0x40000000,0x40000001,
    0x7f000000,0x7f7ffffe,0x7f7fffff,0x7f800000,0x7f800001,0x7fc00000,0x7fffffff]
focus = []
for e in range(-149,128):
    denominator = (1 << (e+149)) if e < -126 else ((e+127) << 23)
    for a in numerators:
        for sa in [0,0x80000000]:
            for sb in [0,0x80000000]:
                aa,bb = a | sa, denominator | sb
                focus.append((3,aa,bb,ref.scalar_expected(3,aa,bb),0,1))
rng = random.Random(990014)
for _ in range(4096):
    a,b = rng.getrandbits(32),rng.getrandbits(32)
    focus.append((3,a,b,ref.scalar_expected(3,a,b),0,1))
records.extend(focus)
inputs = [len(records)]
expected = []
for op,a,b,high,low,valid in records:
    inputs.extend([op,a,b,0]);expected.extend([high,low,valid,0])
(W / 'input.bin').write_bytes(struct.pack('<%dI' % len(inputs), *inputs))
(W / 'expected.bin').write_bytes(struct.pack('<%dI' % len(expected), *expected))
resource = T / 'resources/tests/backend/portable_binary32_reference.bin'
resource.parent.mkdir(parents=True)
shutil.copyfile(fixture_path, resource)
write('fixture.json', {'source_revision':REV,'original_records':ref.RECORD_COUNT,
    'focused_division_records':len(focus),'denominator_exponents':[-149,127],
    'denominator_exponent_count':277,'sign_combinations':4,'numerator_boundary_words':numerators,
    'additional_seeded_fallback_pairs':4096,'complete_records':len(records),
    'complete_output_words':len(expected),'requested_device_bytes':4*(len(inputs)+len(expected)),
    'reference':record(reference_path),'original_fixture':record(fixture_path),
    'input':record(W/'input.bin'),'expected':record(W/'expected.bin'),
    'scope':'Original ordered exact corpus plus separate finite Fraction-derived division boundaries and fallback pairs; not all-input or full retained-path proof.'})
rows=json.loads((B/'compile_commands.json').read_text())
row=next(x for x in rows if x['file'].endswith('/portable_binary32_test.cpp'))
host=shlex.split(row['command']);index=host.index('-o');host[index+1]=str(T/'host.o')
host[1:1]=['-I'+str(S/'candidate/src')]
host += ['-MMD','-MF',str(T/'host.d')]
gpu=['/usr/bin/g++-14','-DSIRIUS_CONTRACT_MODE=2','-DSIRIUS_HAS_RETAINED_COMPUTE=1',
    '-DSIRIUS_RELEASE_RESOURCE_LOCKED=1','-I'+str(R/'src'),'-O3','-DNDEBUG','-std=c++2c',
    '-Wall','-Wextra','-Wpedantic','-Wno-missing-field-initializers','-Werror',
    '-MMD','-MF',str(T/'gpu.d'),'-c',str(W/'gpu.cpp'),'-o',str(T/'gpu.o')]
adapter=B/'src/sirius/backend/CMakeFiles/sirius_backend.dir/vulkan/vulkan_device.cpp.o'
base=B/'src/sirius/libsirius_base.a'
common=[str(p.resolve()) for p in W.glob('*') if p.is_file()]
common += [str(p.resolve()) for p in S.rglob('*') if p.is_file() and p.suffix in ['.h','.slang','.py']]
common += [str(p) for p in [reference_path,fixture_path,resource,R/'tests/backend/portable_binary32_test.cpp',adapter,base]]
common += ['/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so','/usr/lib/x86_64-linux-gnu/libLLVM.so.20.1','/usr/share/vulkan/icd.d/lvp_icd.json']
limits={'seconds':300,'rss_kib':4*1024**2,'reserve_kib':2*1024**2}
actions={}
def action(name,argv,extra=(),software=False,test=False,scope='Private finite arithmetic evidence, no integrated control, performance or qualification verdict.'):
    actions[name]={'revision':REV,'kind':'test' if test else 'diagnostic','expected_exit':0,
        'argv':list(map(str,argv)),'limits':limits,'frozen_inputs':list(dict.fromkeys(common+list(map(str,extra)))),
        'software_provider':software,'scope':scope}
action('compile-host',host)
action('link-host',['/usr/bin/g++-14',T/'host.o',B/'lib/libgtest_main.a',B/'lib/libgtest.a',base,'-ldl','-pthread','-o',T/'host'],[T/'host.o',B/'lib/libgtest_main.a',B/'lib/libgtest.a'])
action('compile-gpu',gpu)
action('link-gpu',['/usr/bin/g++-14',T/'gpu.o',adapter,base,'-lvulkan','-ldl','-pthread','-o',T/'gpu'],[T/'gpu.o'])
cache=(B/'CMakeCache.txt').read_text()
tool=lambda n:re.search(r'^'+n+r':FILEPATH=(.*)$',cache,re.M)[1]
compiler,optimizer,validator,disassembler=[tool(n) for n in ['SIRIUS_SLANGC','SIRIUS_SPIRV_OPT','SIRIUS_SPIRV_VAL','SIRIUS_SPIRV_DIS']]
for producer in ['stock','candidate']:
    action('compile-shader-'+producer,[compiler,W/'probe.slang','-I',S/producer/'src/sirius/kernels','-O0','-target','spirv','-profile','spirv_1_5','-entry','ComputeMain','-stage','compute','-o',T/(producer+'.raw.spv')],[compiler])
    action('optimize-shader-'+producer,[optimizer,'--target-env=vulkan1.2','--ssa-rewrite','--eliminate-dead-code-aggressive','--preserve-bindings','--preserve-interface',T/(producer+'.raw.spv'),'-o',T/(producer+'.spv')],[optimizer,T/(producer+'.raw.spv')])
    action('validate-shader-'+producer,[validator,'--target-env','vulkan1.2',T/(producer+'.spv')],[validator,T/(producer+'.spv')])
    action('disassemble-shader-'+producer,[disassembler,T/(producer+'.spv'),'-o',T/(producer+'.spvasm')],[disassembler,T/(producer+'.spv')])
    action('execute-gpu-'+producer,[T/'gpu',T/(producer+'.spv'),W/'input.bin',W/'expected.bin','{folder}/actual.bin'],[T/'gpu',T/(producer+'.spv')],software=True)
action('execute-host',[T/'host','--gtest_filter=PortableBinary32.MatchesIndependentExactCorpus','--gtest_color=no'],[T/'host'],test=True)
write('actions.json',actions)
write('preparation.json',{'source_revision':REV,'actions':list(actions),'host_test_source_unchanged':True,
    'candidate_header':record(S/'candidate/src/sirius/kernels/portable_binary32.h'),
    'scope':'Private original host exact corpus and raw integer shader arithmetic; no production source adoption, full controls or benefit verdict.'})
print(json.dumps({'prepared':True,'records':len(records),'actions':len(actions)}))
