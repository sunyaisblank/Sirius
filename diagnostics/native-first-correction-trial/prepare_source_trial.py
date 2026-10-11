"""Prepare a guarded native-default source trial without changing production."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, re, shutil, struct, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
TRIAL = WORK/'runtime'
BASELINE = '2053ef79c90752d37b0075b8137cfb9b571b6820'
PRIOR = ROOT/'attestations/native-vulkan/99a4a6d/native-bound-selection-trial/files'

def seal(path):
    data = path.read_bytes()
    return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
def read(path):return json.loads(path.read_text(encoding='utf-8-sig'))
def git(*args):return subprocess.check_output(['git',*args],text=True,timeout=30).strip()

assert ROOT.resolve() == Path('/home/astra/.project/Sirius')
assert WORK == ROOT/'out/native-division-first-correction-screen'
assert git('rev-parse','HEAD') == BASELINE == git('rev-parse','@{u}')
assert not git('status','--porcelain') and not TRIAL.exists()
analysis = read(WORK/'query-analysis.json')
assert analysis['accepted'] and read(WORK/'compiled-review.json')['accepted']
for name, expected in analysis['input_seals'].items():assert seal(Path(name)) == expected,name
for name, expected in analysis['extracted_outputs'].items():assert seal(Path(name)) == expected,name
original = (ROOT/'src/sirius/kernels/retained_triple.slang').read_text()
candidate = (WORK/'candidate/retained_triple.slang').read_text()
body = '''        RetainedTriple residual=RTSubtract(a,RTMultiply(b,q));
        if(residual.valid==0||RPScalarGreater(RPRatioUpper(RPAbs(residual.hi),RPAbs(b.hi)),RPFromBits(0x5d000000u)))
            return RTInvalid();
        q=RTCenter(RTAdd(q,RTConstant(RPDivideSeed(residual.hi,b.hi))));
        if(q.valid==0) return q;
'''
old = '    [loop] for(uint i=0;i<3;++i) {\n'+body+'    }\n'
peeled = '    {\n'+body+'    }\n    [loop] for(uint i=1;i<3;++i) {\n'+body+'    }\n'
guard = '#if defined(SIRIUS_RETAINED_DIVISION_FIRST_CORRECTION32) && !defined(SIRIUS_RETAINED_PORTABLE) && !defined(SIRIUS_RETAINED_FP64) && !defined(SIRIUS_RETAINED_FMA32)\n'
guarded = guard+peeled+'#else\n'+old+'#endif\n'
assert original.count(old) == 1 and candidate == original.replace(old,peeled)
assert not re.search(r'\bi\b',body)
trial_source = original.replace(old,guarded)
assert trial_source.replace(guarded,old) == original
assert trial_source.replace(guarded,peeled) == candidate
builder = (ROOT/'scripts/build-retained-kernels.py').read_text()
point = '    definitions = ["-DSIRIUS_RETAINED_FP64=1"] if fp64 else []\n'
addition = '''    # The first correction establishes centered metadata for the remaining loop.
    # Keep every optional and portable product on its original division body.
    if not portable and not fp64 and not fma and source.stem in ("retained_transport", "retained_endpoint"):
        definitions.append("-DSIRIUS_RETAINED_DIVISION_FIRST_CORRECTION32=1")
'''
assert builder.count(point) == 1
trial_builder = builder.replace(point,point+addition)
assert trial_builder.replace(addition,'') == builder
TRIAL.mkdir()
inputs = {str(ROOT/name):seal(ROOT/name) for name in git('ls-files').splitlines()}
paths = [Path(__file__),WORK/'query-analysis.json',WORK/'compiled-review.json',
         WORK/'source-review.json',WORK/'candidate/retained_triple.slang',
         PRIOR/'baseline-arrays.json',PRIOR/'baseline-native-binding.json',
         ROOT/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h',
         ROOT/'bin/linux-gcc/src/sirius/backend/retained/retained_transport.spv',
         ROOT/'bin/linux-gcc/src/sirius/backend/retained/retained_endpoint.spv']
inputs.update({str(p):seal(p) for p in paths})
for name,text in [('src/sirius/kernels/retained_triple.slang',trial_source),
                  ('scripts/build-retained-kernels.py',trial_builder)]:
    target = TRIAL/'source'/name;target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(text)
    baseline = TRIAL/'baseline'/name;baseline.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(ROOT/name,baseline)
baseline_header = TRIAL/'baseline/bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
baseline_header.parent.mkdir(parents=True,exist_ok=True)
shutil.copyfile(ROOT/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h',baseline_header)
for name in ('retained_transport.spv','retained_endpoint.spv'):
    shutil.copyfile(ROOT/'bin/linux-gcc/src/sirius/backend/retained'/name,TRIAL/'baseline'/name)
for name in ('baseline-arrays.json','baseline-native-binding.json','native-cases.json'):
    shutil.copyfile(PRIOR/name,TRIAL/name)
arrays = {}
for count,name,text in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',baseline_header.read_text(),re.S):
    words = [int(word[:-1],0) for word in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',text)]
    assert len(words) == int(count)
    data = struct.pack('<'+'I'*len(words),*words)
    arrays[name] = {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
assert len(arrays) == 40 and arrays == read(TRIAL/'baseline-arrays.json')['arrays']
assert {name:seal(Path(name)) for name in inputs} == inputs
outputs = {str(p):seal(p) for p in TRIAL.rglob('*') if p.is_file()}
receipt = {'accepted':True,'baseline_revision':BASELINE,'input_seals':inputs,'output_seals':outputs,
           'allowed_changed_arrays':['kTransportShader','kEndpointShader'],'required_unchanged_arrays':38,
           'native_default_guard_reproduces_screen_source':True,'whole_inverse_recovers_original_two_files':True,
           'scope':'Prepared source for reviewed finite native-default trial only. No production writes, numerical execution, benefit, source adoption or full qualification.'}
(TRIAL/'preparation.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'accepted':True,'prepared_source_files':2,'baseline_arrays':40,'input_files':len(inputs),'outputs':len(outputs)}))
