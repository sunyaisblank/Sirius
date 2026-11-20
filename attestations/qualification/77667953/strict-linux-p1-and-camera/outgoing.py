"""Independent outgoing chart diagnostic for the captured exterior phase-space row."""

from decimal import Decimal, localcontext
import json
from pathlib import Path
import struct
import time

from reference import cartesian

WORK = Path(__file__).resolve().parent
raw = json.loads((WORK/'captured.json').read_text())
packet = struct.unpack('<70f', bytes.fromhex(raw['packet_hex']))
record = struct.unpack_from('<8f', bytes.fromhex(raw['previous_record_hex']))
cartesian.MAX_STEPS = 1024
cartesian.MAX_RHS = 4096
cartesian.MAX_METRICS = 29900

with localcontext() as context:
    context.prec = 80
    mass, spin, floor = map(cartesian.exact, (packet[3], packet[4], packet[23]))
    initial_ing = tuple(map(cartesian.exact, record))
    horizon = mass+(mass*mass-spin*spin).sqrt()
    relative_null_floor = cartesian.binary32(0x1e3ce508)
    null_tolerance = cartesian.binary32(0x3a03126f)
    rows = {}
    for displacement in ('1e-18','1e-22'):
        work = {'started':time.monotonic(), 'rhs':0, 'step':0, 'metric':0,
                'maximum_inverse_identity_residual':Decimal(0)}
        metric = cartesian.Cartesian(mass,spin,Decimal(displacement),work)
        q, u = initial_ing[:4], initial_ing[4:]
        _, r = metric.metric(q)
        _, x, y, z = q
        spin2, r2 = spin*spin, r*r
        sigma = r2+spin2*z*z/r2
        radial_speed = (r*x*u[1]+r*y*u[2]+(r2+spin2)*z*u[3]/r)/sigma
        delta = (r-horizon)*(r-(2*mass-horizon))
        d = (r2+spin2)/delta-1
        # The additive time and azimuth origins are selected at this point, so
        # the coordinate position is identical and only the Jacobian is needed.
        # The last reflection evaluates the outgoing metric as a pullback of
        # the existing defining ingoing metric, without a production RHS.
        time_rate = -2*d*radial_speed
        azimuth_rate = -2*spin*d*radial_speed/(r2+spin2)
        reflected = (-q[0], x, -y, z, -(u[0]+time_rate),
                     u[1]-y*azimuth_rate, -(u[2]+x*azimuth_rate), u[3])
        start = metric.summary(initial_ing,Decimal(0),horizon,relative_null_floor,null_tolerance)
        norm = Decimal(start['gkk'])
        mapped = metric.summary(reflected,norm,horizon,relative_null_floor,null_tolerance)
        assert abs(Decimal(mapped['gkk'])-norm) < Decimal('1e-65')
        assert radial_speed < 0
        resolutions = []
        for substeps in (4,8,16):
            h = floor/substeps
            state = reflected
            last = mapped
            for step in range(1,32*substeps+1):
                state = metric.rk4(state,h)
                current = metric.summary(state,norm,horizon,relative_null_floor,null_tolerance)
                if Decimal(current['radius']) < horizon:
                    assert Decimal(last['radius']) > horizon
                    resolutions.append({'substeps_per_original_interval':substeps,
                                        'bracket_affine':[str((step-1)*h), str(step*h)],
                                        'before':last, 'after':current})
                    break
                last = current
            else:
                raise ValueError(f'No crossing within 32 intervals: {substeps=}, {radial_speed=}, {current=} ')
        rows[displacement] = {'radial_speed_at_initial':str(radial_speed), 'ingoing':start,
                              'reflected_outgoing':mapped, 'resolutions':resolutions,
                              'work':{k:str(v) for k,v in work.items() if k!='started'}}

(WORK/'outgoing.json').write_text(json.dumps({
    'schema':'sirius-bardeen-captured-outgoing-chart-diagnostic-v1',
    'governing_reference':'src/sirius/core/metrics/outgoing_kerr_schild.h',
    'raw_failure':cartesian.pin((WORK/'captured.json').read_bytes()),
    'outer_horizon':str(horizon), 'original_interval':str(floor),
    'finite_difference_displacements':rows},indent=2)+'\n')
