from pathlib import Path
import hashlib,json,re,struct
root=Path.cwd(); work=Path(__file__).resolve().parent
new=root/'bin/linux-gcc/sealed-stage-seven-owner/products'
stock=root/'bin/linux-gcc/src/sirius/backend/retained'
old=root/'attestations/software-vulkan/87ae628/portable-stage-seven-owner/files/out/portable-stage-seven-owner-construction/full-builder/products'
seal=lambda p:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
assert seal(stock/'retained_kernels.h')['sha256']=='3f7aca3372bb872934e978fb0270f2375ce8d2361555b35c8b891e42e07161cc'
assert (new/'retained_kernels.h').read_bytes()==(stock/'retained_kernels.h').read_bytes()
def arrays(p):
    result={}
    for count,name,body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',p.read_text(),re.S):
        values=tuple(int(v[:-1],0) for v in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',body))
        assert len(values)==int(count) and name not in result
        result[name]=values
    return result
assert len(arrays(new/'retained_kernels.h'))==40
assert arrays(new/'retained_kernels.h')==arrays(stock/'retained_kernels.h')
products={}; owners={}; sealed={}
for p in sorted(new.glob('*.spv')):
    if p.name.endswith('_stage7_sealed.spv'):sealed[p.name]=seal(p)
    elif p.name.endswith('_stage7_owner.spv'):
        assert p.read_bytes()==(old/p.name).read_bytes(),p.name
        owners[p.name]=seal(p)
    else:
        assert p.read_bytes()==(stock/p.name).read_bytes(),p.name
        products[p.name]=seal(p)
assert (len(products),len(owners),len(sealed))==(33,2,2)
codes=arrays(new/'retained_owner_kernels.h')
assert set(codes)=={'kTransportShader','kEndpointShader','kTransportSealedShader','kEndpointSealedShader'}
for kind in ('transport','endpoint'):
    for mode,name in [('owner','Shader'),('sealed','SealedShader')]:
        data=(new/f'retained_{kind}_portable_stage7_{mode}.spv').read_bytes()
        assert codes[f'k{kind.title()}{name}']==struct.unpack('<'+'I'*(len(data)//4),data)
for n in ('retained_owner_layout.h','retained_owner_layout.slang'):
    assert (new/n).read_bytes()==(old/n).read_bytes(),n
result={'accepted':True,'stock_header':seal(new/'retained_kernels.h'),'stock_arrays_exact':40,'stock_products_exact':products,'checked_owner_products_exact':owners,'new_sealed_products':sealed,'scope':'Generated product identity/validator acceptance only; no numerical/runtime/performance/release claim.'}
(work/'products-verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('accepted','stock_header','stock_arrays_exact','checked_owner_products_exact','new_sealed_products','scope')},indent=2))
