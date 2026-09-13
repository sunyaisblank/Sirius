"""Independent cubic-basis/coordinate-variation dense-segment witnesses."""
import sys
sys.dont_write_bytecode = True
import importlib.util
import json
from pathlib import Path
import struct
import mpmath as mp

folder = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("endpoint_reference", folder / "endpoint_reference.py")
endpoint = importlib.util.module_from_spec(spec)
spec.loader.exec_module(endpoint)
ref, transport = endpoint.ref, endpoint.transport


def sample(row):
    start, end, increments = row[4:44], row[44:84], row[84:104]
    h, s = row[104:106]
    normal, moving = row[106:110], row[110]
    # Use the defining cubic Hermite basis; the device uses a factored secant
    # polynomial. Differentiation below acts on this independent expression.
    def cubic(first, v0, v1, delta, fraction):
        return first + (-2*fraction**3+3*fraction**2)*delta + h*(
            (fraction**3-2*fraction**2+fraction)*v0 + (fraction**3-fraction**2)*v1)
    x = mp.matrix([cubic(start[i],start[i+4],end[i+4],increments[i],s) for i in range(4)])
    k = mp.matrix([mp.diff(lambda z:cubic(start[i],start[i+4],end[i+4],increments[i],z),s)/h for i in range(4)])
    connections = [ref.geometry(list(position),False)[4] for position in [start[:4],end[:4],x]]
    def contract(connection,k,X):
        return mp.matrix([sum(connection[m][a][b]*k[a]*X[b]
                              for a in range(4) for b in range(4)) for m in range(4)])
    output = [*x,*k]
    for column in range(4):
        X0, V0 = mp.matrix(start[8+8*column:12+8*column]),mp.matrix(start[12+8*column:16+8*column])
        X1, V1 = mp.matrix(end[8+8*column:12+8*column]),mp.matrix(end[12+8*column:16+8*column])
        K0 = V0-contract(connections[0],start[4:8],X0)
        K1 = V1-contract(connections[1],end[4:8],X1)
        delta = increments[4+4*column:8+4*column]
        X = mp.matrix([cubic(X0[i],K0[i],K1[i],delta[i],s) for i in range(4)])
        K = mp.matrix([mp.diff(lambda z:cubic(X0[i],K0[i],K1[i],delta[i],z),s)/h for i in range(4)])
        if moving:
            shift = -sum(normal[i]*X[i] for i in range(4))/sum(normal[i]*k[i] for i in range(4))
            X += k*shift
            K -= contract(connections[2],k,k)*shift
        V = K+contract(connections[2],k,X)
        output.extend([*X,*V])
    return output


def main():
    mp.mp.dps=105
    originals=json.loads((folder/'reference_cases.json').read_text())['cases']
    cases=[]
    for original in originals:
        decoded=[]
        for i in range(46):
            hi,lo,tail,_,_=struct.unpack('<ffffI',struct.pack('<5I',*original['input'][5*i:5*i+5]))
            decoded.append(mp.mpf(hi)+mp.mpf(lo)+mp.mpf(tail))
        ref.params=decoded[:4]
        chart=decoded[44]
        ref.reflection=[chart,1,chart,1]
        # Exact encoded endpoint inputs define this independent interpolation
        # problem; a separate RK witness judges the ODE truncation error.
        _, start=endpoint.project(mp.matrix(decoded[4:44]))
        _, end=endpoint.project(mp.matrix(list(map(mp.mpf,original['reference'][:40]))))
        start,end=start[40:],end[40:]
        increments=[end[i]-start[i] for start_index in [0,8,16,24,32] for i in range(start_index,start_index+4)]
        for fraction,moving in [(mp.mpf('0.5'),0),(mp.mpf('0.75'),1),(mp.mpf(1),0)]:
            # A tilted plane avoids assuming radial or axis-aligned events.
            normal=list(map(mp.mpf,[0,1,mp.mpf('0.25'),mp.mpf('-0.125')])) if moving else [mp.mpf(0)]*4
            row=[*decoded[:4],*start,*end,*increments,decoded[45],fraction,*normal,mp.mpf(moving),chart]
            samples=[]
            for precision in [75,105]:
                mp.mp.dps=precision
                samples.append(sample(row))
            gaps=[abs(a-b) for a,b in zip(*samples)]
            assert max(gaps)<mp.mpf('1e-50')
            cases.append({'name':original['name']+'-'+('arrival' if moving else str(fraction)),
                          'input':[word for v in row for word in transport.pair(v)],
                          'reference':[str(v) for v in samples[1]], 'precision_gap':[str(v) for v in gaps]})
            print(cases[-1]['name'],'gap',max(gaps),flush=True)
    (folder/'dense_reference.json').write_text(json.dumps({'precision':[75,105],'cases':cases},indent=2)+'\n')
    lines=['// Independent cubic Hermite and coordinate arrival derivatives.', '// clang-format off',
           '#pragma once','#include <array>','#include <cstdint>','namespace sirius::test::retained_dense {',
           'struct Wide {double high,low;};',
           'struct Case {const char* name; std::array<std::uint32_t,560> input; std::array<Wide,40> reference;};',
           f'inline constexpr std::array<Case,{len(cases)}> cases{{{{']
    for case in cases:
        lines.append('{'+json.dumps(case['name'])+',{{')
        for i in range(0,560,16):lines.append(','.join(str(v)+'u' for v in case['input'][i:i+16])+',')
        lines.append('}},{{')
        for value in case['reference']:
            v=mp.mpf(value);hi=float(v);lo=float(v-mp.mpf(hi))
            lines.append('{'+hi.hex()+','+lo.hex()+'},')
        lines.append('}}},')
    lines.extend(['}};','} // namespace sirius::test::retained_dense','// clang-format on',''])
    (folder/'dense_reference.h').write_text('\n'.join(lines))


if __name__=='__main__':main()
