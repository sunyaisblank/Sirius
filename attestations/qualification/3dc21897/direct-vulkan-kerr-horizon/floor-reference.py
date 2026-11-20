"""Re-evaluate the second Bardeen rejection from exact reflected state bytes."""

from decimal import Decimal as D, localcontext
import importlib.util
import json
from pathlib import Path
import struct
import time

ROOT=Path.cwd()
WORK=ROOT/'out/kerr-outgoing-trace'
PRIOR=ROOT/'attestations/scientific-core/02ee51f/kerr-cartesian-floor-diagnostic/cartesian_trial.py'
spec=importlib.util.spec_from_file_location('reference',PRIOR)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.MAX_STEPS=254
module.MAX_RHS=1016
module.MAX_METRICS=7130
raw=json.loads((WORK/'second-floor.json').read_text())
packet=struct.unpack('<70f',bytes.fromhex(raw['packet_hex']))
record=bytes.fromhex(raw['previous_record_hex'])
assert packet[3]==1 and packet[5]==packet[6]==0
assert raw['active']==1 and raw['physical_submissions']==2400
assert struct.unpack_from('<4I',record,240)==(1,0,936,725)
assert struct.unpack_from('<f',record,192)[0]==packet[23]
assert struct.unpack_from('<f',record,16)[0]>0  # Recorded reflected orientation.
measurements={}
with localcontext() as ctx:
    ctx.prec=80
    exact=module.exact
    start=tuple(map(exact,struct.unpack_from('<8f',record)))
    M,a,interval=map(exact,(packet[3],packet[4],packet[23]))
    horizon=M+(M*M-a*a).sqrt()
    null_scale=module.binary32(0x1e3ce508)
    null_tolerance=module.binary32(0x3a03126f)
    error_tolerance=module.binary32(0x37a7c5ac)
    for displacement in ('1e-18','1e-22'):
        stats={'started':time.monotonic(),'rhs':0,'step':0,'metric':0,
               'maximum_inverse_identity_residual':D(0)}
        metric=module.Cartesian(M,a,D(displacement),stats)
        init=metric.summary(start,D(0),horizon,null_scale,null_tolerance)
        last=None
        results={}
        for n in (2,4,8,16,32,64,128):
            end=start
            for _ in range(n):
                end=metric.rk4(end,interval/n)
            data=metric.summary(end,D(init['gkk']),horizon,null_scale,null_tolerance)
            results[n]={'relative_null':data['controller_style_unprojected_null_residual'],
                        'radius_minus_horizon':data['radius_minus_independent_buffered_capture'],
                        'endpoint_tangent':data['u'],
                        'endpoint_position':data['q'],
                        'error_ratio':None if last is None else module.trial_comparator(
                            start,last,end,error_tolerance)['diagnostic_ratio_to_binary32_2e_5']}
            last=end
        measurements[displacement]={'initial':init,'interval':str(interval),'endpoints':results}
artifact={'reflected_chart':True,'independent_backward_ray_reference':measurements,
          'archived_reference':{'path':str(PRIOR.relative_to(ROOT)),**module.pin(PRIOR.read_bytes())}}
(WORK/'second-floor-reference.json').write_text(json.dumps(artifact,indent=2)+'\n')
