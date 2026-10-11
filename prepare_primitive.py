"""Compile the actual upward-product owner and independent exact expectations."""
from collections import Counter
from pathlib import Path
import hashlib, importlib.util, itertools, json, random, re, shlex, struct, subprocess, sys
sys.dont_write_bytecode = True
W = Path(__file__).resolve().parent
R = W.parents[1]
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module
ref = load('fraction_reference', R/'tests/support/portable_binary32/reference.py')
builder = load('product_builder', R/'scripts/build-retained-kernels.py')
cache = {}
for line in (R/'bin/linux-gcc/CMakeCache.txt').read_text().splitlines():
    match = re.match(r'([^#/:][^:]*):[^=]*=(.*)$', line)
    if match: cache[match[1]] = match[2]
tools = {kind: cache[key] for kind, key in [('compiler','SIRIUS_SLANGC'),
    ('assembler','SIRIUS_SPIRV_AS'),('disassembler','SIRIUS_SPIRV_DIS'),
    ('validator','SIRIUS_SPIRV_VAL'),('optimizer','SIRIUS_SPIRV_OPT')]}
commands = []
def run(args):
    args = list(map(str,args)); commands.append(args)
    subprocess.run(args,check=True,cwd=R/'bin/linux-gcc')
# A separate immutable baseline import uses the exact old pair owner. Every
# other dependency remains the current, unchanged production source.
base = W/'baseline-kernels'; base.mkdir()
for name in ['retained_scalar.slang','portable_binary32.h']:
    (base/name).write_bytes((R/'src/sirius/kernels'/name).read_bytes())
(base/'retained_pair.slang').write_bytes((W/'baseline-retained_pair.slang').read_bytes())
for label, source in [('baseline',base),('candidate',R/'src/sirius/kernels')]:
    raw=W/(label+'-upward.compiler.spv'); opt=W/(label+'-upward.optimized.spv')
    assembly=W/(label+'-upward.spvasm'); final=W/(label+'-upward.spv')
    run([tools['compiler'],W/'upward-probe.slang','-DSIRIUS_RETAINED_PORTABLE=1',
        '-DSIRIUS_RETAINED_NORMAL_SUM32=1','-DSIRIUS_RETAINED_PARALLEL_TRANSPORT=1',
        '-I',source,'-O0','-target','spirv','-profile','spirv_1_5','-entry','main',
        '-stage','compute','-o',raw])
    run([tools['optimizer'],'--target-env=vulkan1.2','--ssa-rewrite',
        '--eliminate-dead-code-aggressive','--preserve-bindings','--preserve-interface',raw,'-o',opt])
    run([tools['disassembler'],opt,'-o',assembly])
    text=assembly.read_text()
    if label=='candidate': text=builder.portable_normal_sum_controls(text,upmultiply=True)
    else:
        assert 'OpTypeFloat' not in text and re.findall(r'OpCapability (\S+)',text)==['Shader']
    assembly.write_text(text)
    run([tools['assembler'],'--target-env','spv1.5',assembly,'-o',final])
    run([tools['validator'],'--target-env','vulkan1.2',final])
fixture=R/'bin/linux-gcc/tests/backend/resources/tests/backend/portable_binary32_reference.bin'
words=struct.unpack('<'+'I'*(fixture.stat().st_size//4),fixture.read_bytes())
assert words[:2]==(ref.MAGIC,ref.VERSION) and len(words)==3+6*words[2]
pairs={(words[i+1],words[i+2]) for i in range(3,len(words),6) if words[i]==2}
edges=[0,1,2,0x007ffffe,0x007fffff,0x00800000,0x00800001,0x3f7fffff,
    0x3f800000,0x3f800001,0x3f800003,0x3fc00000,0x5e7fffff,0x5e800000,
    0x5e800001,0x7f000000,0x7f7fffff,ref.INF,0x7f800001,ref.NAN]
edges=[v|sign for v in edges for sign in [0,ref.SIGN]]
pairs.update(itertools.product(edges,repeat=2))
mantissas=[0,1,2,3,0x3fffff,0x400000,0x7ffffe,0x7fffff]
for total in [127,128,129]:
    for ae in range(1,190):
        be=total-ae
        if not 1<=be<=189: continue
        for am,bm in itertools.product(mantissas,repeat=2):
            a=(ae<<23)|am; b=(be<<23)|bm
            if a<=0x5e800000 and b<=0x5e800000: pairs.add((a,b))
# Halfway products with either parity, exact powers and neighboring words at
# normal, ceiling and fallback exponents. Expectations use exact rationals.
for ae in [1,2,63,64,65,126,127,128,188,189]:
    for be in [1,2,63,64,65,126,127,128,188,189]:
        for am,bm in itertools.product(mantissas,repeat=2):
            pairs.add(((ae<<23)|am,(be<<23)|bm))
rng=random.Random(0x525055504d554c)
for _ in range(8192):
    pairs.add((rng.getrandbits(32),rng.getrandbits(32)))
    pairs.add((rng.randrange(0x00800000,0x5e800001),rng.randrange(0x00800000,0x5e800001)))
inputs=[len(pairs)]; expected=[]; classes=Counter(); sums=Counter(); ties=Counter()
for ordinal,(a,b) in enumerate(sorted(pairs)):
    inputs.extend([ordinal,a,b]); aa=a&ref.MASK; bb=b&ref.MASK
    domain=0x00800000<=a<=0x5e800000 and 0x00800000<=b<=0x5e800000 and (a>>23)+(b>>23)>=128
    if aa==0 or bb==0: result=0; category='zero_before_bounds'
    elif aa>0x5e800000 or bb>0x5e800000: result=0x7f000000; category='bounds_sentinel'
    else:
        exact=ref.exact(a)*ref.exact(b)
        result=(ref.round_exact(exact,(a^b)&ref.SIGN)+1)&0xffffffff
        category='positive_normal_guard' if domain else 'integer_fallback'
        if domain:
            unit=max(ref.exponent_of(exact)-23,-149)
            scaled=exact/(2**unit) if unit>=0 else exact*(1<<-unit)
            q,r=divmod(scaled.numerator,scaled.denominator)
            if 2*r==scaled.denominator: ties['odd' if q&1 else 'even']+=1
    expected.extend([result,int(domain),a,b,ordinal]); classes[category]+=1
    if a<ref.SIGN and b<ref.SIGN: sums[str((a>>23)+(b>>23))]+=1
for key in ['zero_before_bounds','bounds_sentinel','positive_normal_guard','integer_fallback']: assert classes[key]>0
assert ties['odd']>0 and ties['even']>0
assert all(sums[str(total)]>0 for total in [127,128,129])
for name,data in [('upward-input.bin',inputs),('upward-expected.bin',expected)]:
    (W/name).write_bytes(struct.pack('<'+'I'*len(data),*data))
assert (len(inputs)+len(expected))*4<8*1024*1024
cc=json.loads((R/'bin/linux-gcc/compile_commands.json').read_text())
template=next(item for item in cc if item['file'].endswith('/retained_compute.cpp'))
args=shlex.split(template['command']); args[args.index('-o')+1]=str(W/'gpu_probe.o'); args[-1]=str(W/'gpu_probe.cpp'); run(args)
line=subprocess.check_output(['ninja','-C',str(R/'bin/linux-gcc'),'-t','commands','sirius_backend_tests'],text=True).splitlines()[-1]
link=shlex.split(line.split('&&')[1]); link=[value for value in link if not value.endswith('.o')]
link.insert(1,str(W/'gpu_probe.o')); link[link.index('-o')+1]=str(W/'gpu_probe'); run(link)
record={'scope':'Finite independent exact-Fraction actual-function corpus; no universal equivalence, speed or qualification claim',
    'cases':len(pairs),'words_per_variant':len(expected),'classes':dict(classes),'admitted_halfway_parities':dict(ties),
    'exponent_sum_boundary_counts':{str(total):sums[str(total)] for total in [127,128,129]},
    'reference_sha256':hashlib.sha256((R/'tests/support/portable_binary32/reference.py').read_bytes()).hexdigest(),
    'fixture_sha256':hashlib.sha256(fixture.read_bytes()).hexdigest(),'commands':commands,
    'source_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()}
(W/'upward-oracle.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({key:record[key] for key in ['cases','classes','admitted_halfway_parities','exponent_sum_boundary_counts']}))
