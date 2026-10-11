"""Emit one unadopted native-narrow product branch and bounded builder changes."""
from pathlib import Path
import hashlib
import json
import shutil
import sys
sys.dont_write_bytecode = True
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent

def seal(path):
    raw = path.read_bytes()
    return {'path': str(path.relative_to(ROOT)), 'bytes': len(raw),
            'sha256': hashlib.sha256(raw).hexdigest()}

model = json.loads((WORK/'model-result.json').read_text())
assert model['pass'] and model['device_workloads'] == 0
assert json.loads((WORK/'model-independent-review.json').read_text())['accepted']
source = ROOT/'src/sirius/kernels'
kernels = WORK/'kernels'
kernels.mkdir()
original_inputs = []
for path in sorted(source.iterdir()):
    if path.is_file():
        assert not path.is_symlink()
        original_inputs.append(seal(path))
        shutil.copyfile(path, kernels/path.name)

pair = kernels/'retained_pair.slang'
original = pair.read_text()
anchor = '''#else
 RPScalar ah,al,bh,bl;Split32(a,ah,al);Split32(b,bh,bl);hi=OrderedMultiply32(a,b);'''
addition = '''#else
#ifdef SIRIUS_RETAINED_FMA32
 // Only nonzero exact residuals on the bounded normal lattice substitute for
 // the default split product. Zero signs and all other inputs retain its body.
 uint am=RPBits(a)&0x7fffffffu,bm=RPBits(b)&0x7fffffffu;
 if(am>=0x2b800000u&&am<=0x5d800000u&&bm>=0x2b800000u&&bm<=0x5d800000u) {
  RPScalar high=OrderedMultiply32(a,b);
  RPScalar negativeHigh=RPScalarNegate(high);
  RPScalar residual=spirv_asm {
   OpFmaKHR $$float %v $a $b $negativeHigh;
   OpDecorate %v NoContraction;
   OpCopyObject $$float result %v;
  };
  if((RPBits(residual)&0x7fffffffu)!=0u) {hi=high;lo=residual;return;}
 }
#endif
 RPScalar ah,al,bh,bl;Split32(a,ah,al);Split32(b,bh,bl);hi=OrderedMultiply32(a,b);'''
assert original.count(anchor) == 1
changed = original.replace(anchor, addition)
assert changed.replace(addition, anchor) == original
pair.write_text(changed)
assert all((kernels/p.name).read_bytes()==p.read_bytes() for p in source.iterdir()
           if p.is_file() and p.name!='retained_pair.slang')

builder_source = ROOT/'scripts/build-retained-kernels.py'
original_builder = builder_source.read_text()
changes = [
 ('def native_controls(text, fma=False):',
  'def native_controls(text, fma=False, signed_zero=True):'),
 ('        capabilities += "\\n               OpCapability FMAKHR\\n               OpCapability SignedZeroInfNanPreserve"',
  '        capabilities += "\\n               OpCapability FMAKHR"\n    if fma and signed_zero:\n        capabilities += "\\n               OpCapability SignedZeroInfNanPreserve"'),
 ('    if fma:\n        modes += "\\n               OpExecutionMode " + entry + " SignedZeroInfNanPreserve 32"',
  '    if fma and signed_zero:\n        modes += "\\n               OpExecutionMode " + entry + " SignedZeroInfNanPreserve 32"'),
 ('        if not fp64 or portable or source.stem not in ("retained_transport", "retained_endpoint"):',
  '        if portable or source.stem not in ("retained_transport", "retained_endpoint"):'),
 ('            raise ValueError("FMA32 is qualified only for native-wide retained Transport and Endpoint")',
  '            raise ValueError("diagnostic FMA32 belongs only to native retained Transport and Endpoint")'),
 ('    text = native_controls(text, fma=fma)',
  '    text = native_controls(text, fma=fma, signed_zero=fp64)'),
]
builder = original_builder
for old, new in changes:
    assert builder.count(old) == 1, old
    builder = builder.replace(old, new)
reverse = builder
for old, new in reversed(changes):
    assert reverse.count(new) == 1, new
    reverse = reverse.replace(new, old)
assert reverse == original_builder
(WORK/'builder.py').write_text(builder)
result = {'pass':True,'production_adopted':False,
          'source_revision':model['source_revision'],'original_inputs':original_inputs,
          'original_builder':seal(builder_source),'emitter':seal(Path(__file__).resolve()),
          'candidate_pair':seal(pair),'candidate_builder':seal(WORK/'builder.py'),
          'unchanged_imports_exact':True,'source_reverse_exact':True,'builder_reverse_exact':True,
          'builder_replacements':changes,
          'scope':'Only native default FMA32 RetainedProduct32 gains guarded nonzero substitution; full old fallback and wide/portable sources unchanged. Builder adds only optional FMA capability/extension to default modes, no global SignedZero mode.'}
(WORK/'emission.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k: result[k] for k in ('pass','candidate_pair','candidate_builder','scope')},indent=2))
