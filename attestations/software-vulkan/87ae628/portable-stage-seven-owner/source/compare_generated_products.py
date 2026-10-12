"""Join the complete generated default products to their protected baseline."""
import hashlib, json, re, struct, sys
from pathlib import Path
sys.dont_write_bytecode = True
root = Path.cwd()
work = Path(__file__).resolve().parent
new = work / 'full-builder/products'
old = root / 'bin/linux-gcc/src/sirius/backend/retained'
seal = lambda p: {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
assert seal(old/'retained_kernels.h')['sha256'] == '3f7aca3372bb872934e978fb0270f2375ce8d2361555b35c8b891e42e07161cc'
assert (new/'retained_kernels.h').read_bytes() == (old/'retained_kernels.h').read_bytes()
def arrays(header):
    result = {}
    for count,name,body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',header.read_text(),re.S):
        values = tuple(int(v[:-1],0) for v in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',body))
        assert len(values) == int(count) and name not in result
        result[name] = values
    return result
assert arrays(new/'retained_kernels.h') == arrays(old/'retained_kernels.h')
assert len(arrays(new/'retained_kernels.h')) == 40
stock = {}; owners = {}
for binary in sorted(new.glob('*.spv')):
    if binary.name.endswith('_stage7_owner.spv'):
        owners[binary.name] = seal(binary)
    else:
        assert binary.read_bytes() == (old/binary.name).read_bytes(), binary.name
        stock[binary.name] = seal(binary)
assert len(stock) == 33 and len(owners) == 2
codes = arrays(new/'retained_owner_kernels.h')
assert set(codes) == {'kTransportShader','kEndpointShader'}
for kind in ('transport','endpoint'):
    data = (new/f'retained_{kind}_portable_stage7_owner.spv').read_bytes()
    assert codes[f'k{kind.title()}Shader'] == struct.unpack('<'+'I'*(len(data)//4), data)
for name in ('retained_owner_layout.h','retained_owner_layout.slang'):
    assert (new/name).read_bytes() == (work/'owner-layout'/name).read_bytes()
result = {'accepted':True,'stock_header':seal(new/'retained_kernels.h'), 'stock_arrays_exact':40,
          'stock_products_exact':stock,'new_owner_products':owners,
          'new_generated_headers':{name:seal(new/name) for name in ('retained_owner_kernels.h','retained_owner_layout.h','retained_owner_layout.slang')},
          'scope':'Complete source/build product comparison only. Canonical preset files remain untouched. No host test, real allocation, replay, gain, frame, science or release qualification.'}
(work/'generated-products-comparison.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('accepted','stock_header','stock_arrays_exact','new_owner_products','scope')},indent=2))
