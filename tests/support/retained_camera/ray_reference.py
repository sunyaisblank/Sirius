"""Freeze physical lens seeds for the independent complete-camera witnesses."""
import sys
sys.dont_write_bytecode = True
import importlib.util
import json
from pathlib import Path
import struct
import mpmath as mp

folder=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('transport_reference',folder.parent/'retained_transport/reference.py')
transport=importlib.util.module_from_spec(spec)
spec.loader.exec_module(transport)


def seed(row):
    row=list(map(mp.mpf,row))
    thin=row[29]
    D=1+thin*(row[24]-1)
    def direction(x,y,right,up):
        q=[D,D*row[23]*(1-2*y/row[21])-up*thin,
           D*row[23]*(2*x/row[20]-1)*row[22]-right*thin]
        norm=mp.sqrt(sum(v*v for v in q))
        return [v/norm for v in q]
    sample=tuple(row[25:29])
    n=direction(*sample)
    derivatives=[]
    for axis in range(3):
        for column in range(4):
            order=[0]*4;order[column]=1
            derivatives.append(mp.diff(lambda *v:direction(*v)[axis],sample,tuple(order)))
    return [*row[:20],row[27]*thin,row[28]*thin,*n,*derivatives,0,0,thin,0,0,0,0,thin]


def main():
    originals=json.loads((folder/'reference_cases.json').read_text())['cases']
    continuous=json.loads((folder/'continuous_reference.json').read_text())['cases']
    cases=[]
    for original in originals+continuous:
        samples=[]
        for precision in [100,180]:
            mp.mp.dps=precision
            samples.append(seed(original['input']))
        assert max(abs(a-b) for a,b in zip(*samples))<mp.mpf('1e-90')
        packet=[word for value in samples[1] for word in transport.pair(mp.mpf(value))]
        assert len(packet)==225
        reference=original.get('scientific',original.get('reference'))
        gaps=original.get('precision_gap',original.get('reference_gap'))
        assert reference is not None and gaps is not None
        cases.append({'name':original['name'],'input':packet,'reference':reference,'reference_gap':gaps})
        print(original['name'],flush=True)
    (folder/'ray_reference.json').write_text(json.dumps({'precision':[100,180],'cases':cases},indent=2)+'\n')
    lines=['// Independent physical lens seeds and full camera witnesses. See ray_reference.py.',
           '// clang-format off','#pragma once','#include <array>','#include <cstdint>',
           'namespace sirius::test::retained_ray_camera {',
           'struct Case {const char* name; std::array<std::uint32_t,225> input; std::array<long double,104> reference,reference_gap;};',
           f'inline constexpr std::array<Case,{len(cases)}> cases{{{{']
    for case in cases:
        lines.append('{'+json.dumps(case['name'])+',{{')
        for i in range(0,225,16):lines.append(','.join(str(v)+'u' for v in case['input'][i:i+16])+',')
        lines.append('}},{{')
        lines.append(',\n'.join(v+'L' for v in case['reference']))
        lines.append('}},{{')
        lines.append(',\n'.join(v+'L' for v in case['reference_gap']))
        lines.append('}}},')
    lines.extend(['}};','} // namespace sirius::test::retained_ray_camera','// clang-format on',''])
    (folder/'ray_reference.h').write_text('\n'.join(lines))


if __name__=='__main__':main()
