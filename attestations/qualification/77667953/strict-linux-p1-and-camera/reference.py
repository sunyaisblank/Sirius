"""Replay the exact failed Bardeen interval with the archived off-null Cartesian reference."""

from decimal import Decimal, localcontext
import importlib.util
import json
from pathlib import Path
import struct
import sys
import time

ROOT = Path.cwd()
WORK = ROOT/'out/strict-p1-camera/floor'
PRIOR = ROOT/'attestations/scientific-core/02ee51f/kerr-cartesian-floor-diagnostic/cartesian_trial.py'
source = importlib.util.spec_from_file_location('archived_cartesian', PRIOR)
cartesian = importlib.util.module_from_spec(source)
source.loader.exec_module(cartesian)

log = WORK.parent/'device/runs/p1-legacy-bardeen/stdout.log'
lines = [line[line.index('{"schema"'):] for line in log.read_text().splitlines()
         if '{"schema":"sirius-trace-continuation-failure-v1"' in line]
assert lines and len(set(lines)) == 1
failure = json.loads(lines[0])
assert failure['scalar_bytes'] == failure['packet_scalar_bytes'] == 4
assert failure['record_bytes'] == 272 and failure['byte_order'] == 'little'
assert failure['action'] == 1 and failure['active'] == 1
old = cartesian.raw_hex(failure['previous_record_hex'], 272)
new = cartesian.raw_hex(failure['current_record_hex'], 272)
packet_raw = cartesian.raw_hex(failure['packet_hex'], 280)
packet = struct.unpack('<70f', packet_raw)
before = struct.unpack_from('<4I', old, 240)
after = struct.unpack_from('<4I', new, 240)
assert before == (1, 0, 287, 280) and after == (2, 5, 288, 280)
assert old[:240] == new[:240] and old[256:] == new[256:]
assert list(struct.unpack_from('<4I', old, 256)) == failure['expected_identity']
assert packet[2] == packet[5] == packet[6] == 0 and packet[3] == packet[26] == 1
assert packet[21] == 30000 and struct.unpack_from('<f', old, 192)[0] == packet[23]
state = struct.unpack_from('<8f', old)
sample = (*packet[2:7], *state[1:4], *state[4:], packet[23], 1.0, packet[23], state[0])
(WORK/'sample.bin').write_bytes(struct.pack('<16f', *sample))
(WORK/'captured.json').write_text(json.dumps(failure, indent=2)+'\n')

if len(sys.argv) > 1 and sys.argv[1] == '--extract':
    sys.exit(0)

# Cover both finite-difference displacements with at most 126 RK4 steps each.
# The finest discrete solution is used to measure convergence, without projection.
cartesian.MAX_STEPS = 126
cartesian.MAX_RHS = 504
cartesian.MAX_METRICS = 3600
results = {}
exact = cartesian.exact
with localcontext() as context:
    context.prec = 80
    M, a, interval, factor = map(exact, (packet[3], packet[4], packet[23], packet[26]))
    radius = factor*(M+(M*M-a*a).sqrt())
    tol, null_tol, null_floor = (cartesian.binary32(word)
                                for word in (0x37a7c5ac, 0x3a03126f, 0x1e3ce508))
    start = tuple(map(exact, state))
    for displacement in ('1e-18','1e-22'):
        work = {'started':time.monotonic(), 'rhs':0, 'step':0, 'metric':0,
                'maximum_inverse_identity_residual':Decimal(0)}
        model = cartesian.Cartesian(M,a,Decimal(displacement),work)
        initial = model.summary(start, Decimal(0), radius, null_floor,null_tol)
        norm = Decimal(initial['gkk'])
        endpoints = {}
        summaries = {}
        for count in (2,4,8,16,32,64):
            end = start
            for _ in range(count):
                end = model.rk4(end, interval/Decimal(count))
            endpoints[count] = end
            data = model.summary(end,norm,radius,null_floor,null_tol)
            if count//2 in endpoints:
                data['production_scaled_comparison'] = cartesian.trial_comparator(
                    start,endpoints[count//2],end,tol)
            summaries[count] = data
        results[displacement] = {'start':initial,'endpoints_by_substeps':summaries,
                                 'work':{key:str(value) for key,value in work.items() if key!='started'}}
output = {'schema':'sirius-bardeen-exact-floor-independent-cartesian-v1',
          'mass_shell_projection':False,'interval':str(interval),'buffered_capture_radius':str(radius),
          'archived_reference':{'path':str(PRIOR.relative_to(ROOT)),**cartesian.pin(PRIOR.read_bytes())},
          'raw_witness':cartesian.pin((WORK/'captured.json').read_bytes()),
          'sample_bytes': (WORK/'sample.bin').read_bytes().hex(),
          'finite_difference_displacements':results}
(WORK/'decimal.json').write_text(json.dumps(output, indent=2)+'\n')
