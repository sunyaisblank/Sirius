"""Prepare a private, exact-body first-correction peeling compiler screen."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, re, shutil, struct, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
BASELINE = '2053ef79c90752d37b0075b8137cfb9b571b6820'
PRIOR = ROOT/'attestations/native-vulkan/451e2cb/native-division-phase-observer/files'
ORDINARY = ROOT/'attestations/native-vulkan/3045067/ordinary-pipeline-generated-code'

def seal(path):
    data = path.read_bytes()
    return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
def read(path):return json.loads(path.read_text(encoding='utf-8-sig'))
def git(*args):return subprocess.check_output(['git',*args],text=True,timeout=30).strip()

assert ROOT.resolve() == Path('/home/astra/.project/Sirius')
assert WORK == ROOT/'out/native-division-first-correction-screen'
assert git('rev-parse','HEAD') == BASELINE == git('rev-parse','@{u}')
assert not git('status','--porcelain')
baseline = read(PRIOR/'baseline-source.json')
assert baseline['revision'] == BASELINE
assert set(git('ls-files').splitlines()) == set(baseline['tracked_source_seals'])
for name,expected in baseline['tracked_source_seals'].items():assert seal(ROOT/name) == expected,name
header = ROOT/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
assert seal(header) == baseline['production_header']
entry = ORDINARY/'query-payloads/transport.spv'
assert seal(entry) == {'bytes':85364,'sha256':'b52e6373f6b30aac6c7fe877f6e07f3c47c7ac9eafa1671840eacd19960d7cb1'}
match = re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> kTransportShader\{\{(.*?)\}\};',header.read_text(),re.S)
assert len(match) == 1
count,array = match[0]
words = [int(w[:-1],0) for w in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',array)]
assert len(words) == int(count) and struct.pack('<'+'I'*len(words),*words) == entry.read_bytes()

original = (ROOT/'src/sirius/kernels/retained_triple.slang').read_text()
assert seal(ROOT/'src/sirius/kernels/retained_triple.slang')['sha256'] == '5a4a948a2527526644ec7b82a6e18409577debe997ae47be899cdecb8f4654bd'
function = re.search(r'RetainedTriple RTDivide\(RetainedTriple a,RetainedTriple b\) \{\n#ifdef SIRIUS_RETAINED_PORTABLE\n.*?#else\n(.*?)#endif\n\}',original,re.S)
assert function
native = function[1]
body = '''        RetainedTriple residual=RTSubtract(a,RTMultiply(b,q));
        if(residual.valid==0||RPScalarGreater(RPRatioUpper(RPAbs(residual.hi),RPAbs(b.hi)),RPFromBits(0x5d000000u)))
            return RTInvalid();
        q=RTCenter(RTAdd(q,RTConstant(RPDivideSeed(residual.hi,b.hi))));
        if(q.valid==0) return q;
'''
old = '    [loop] for(uint i=0;i<3;++i) {\n'+body+'    }\n'
new = '    {\n'+body+'    }\n    [loop] for(uint i=1;i<3;++i) {\n'+body+'    }\n'
assert native.count(old) == original.count(old) == 1
assert not re.search(r'\bi\b',body)
assert native.count('RetainedTriple q=RTConstant(RPDivideSeed(a.hi,b.hi));') == 1
candidate = original.replace(old,new)
assert candidate.replace(new,old) == original
assert candidate.count(new) == 1
assert candidate.split('#else\n',1)[0] == original.split('#else\n',1)[0]
destination = WORK/'candidate'
assert not destination.exists() and not (WORK/'preparation.json').exists()
destination.mkdir()
sources = sorted((ROOT/'src/sirius/kernels').glob('retained_*.slang'))
sources.append(ROOT/'src/sirius/kernels/portable_binary32.h')
for path in sources:shutil.copyfile(path,destination/path.name)
(destination/'retained_triple.slang').write_text(candidate)
for path in sources:
    if path.name != 'retained_triple.slang':assert seal(path) == seal(destination/path.name)
paths = [Path(__file__),PRIOR/'baseline-source.json',header,entry,
         ORDINARY/'source-code-join-bindings.json',ORDINARY/'analysis.json',ORDINARY/'manifest.json',
         ORDINARY/'instruction-encoding-join.json',ORDINARY/'generated-code/transport.llvmir',
         ORDINARY/'query-payloads/transport.info1',ORDINARY/'query-payloads/transport.info2',
         ROOT/'scripts/build-retained-kernels.py',*sources]
record = {'accepted':True,'baseline_revision':BASELINE,'input_seals':{str(p):seal(p) for p in paths},
    'output_seals':{str(destination/p.name):seal(destination/p.name) for p in sources},
    'source_paths_exact':len(baseline['tracked_source_seals']),'production_header':seal(header),
    'original_transport_spv':seal(entry),'source_inverse_byte_exact':True,
    'loop_body_byte_exact':True,'body_variable_i_unused':True,
    'correction_count':3,'peeled_iterations':[0],'remaining_iterations':[1,2],
    'constructor_calls_unchanged':True,'guards_returns_final_residual_unchanged':True,
    'scope':'Private compiler screen only. First correction body copied verbatim, remaining two corrections unchanged. Initial RTConstant/RTPack and low limbs preserved; subsequent RTCenter supplies literal zero error. No source adoption, runtime dispatch, exclusive cost or ordinary performance/qualification claim.'}
assert git('rev-parse','HEAD') == BASELINE and not git('status','--porcelain')
(WORK/'preparation.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:record[k] for k in ('accepted','baseline_revision','source_paths_exact','source_inverse_byte_exact','loop_body_byte_exact','scope')}))
