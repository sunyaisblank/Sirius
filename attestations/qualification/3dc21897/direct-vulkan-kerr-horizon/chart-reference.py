"""Independent high-precision oracle for the reflected outgoing chart and ellipsoid event."""

from decimal import Decimal as D, localcontext
import hashlib
import json
from pathlib import Path
import struct
import sys

WORK = Path.cwd() / 'out/kerr-outgoing-trace'
PI = D('3.141592653589793238462643383279502884197169399375105820974944592307816406286')
EXTENT = {'32': 'f', '64': 'd'}


def sha(path):
    return {'path':str(path.relative_to(Path.cwd())), 'bytes': path.stat().st_size,
            'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def trigonometry(angle):
    angle %= 2*PI
    sine = sine_term = angle
    cosine = cosine_term = D(1)
    for i in range(1,70):
        sine_term *= -(angle*angle)/(2*i*(2*i+1))
        cosine_term *= -(angle*angle)/(2*i*(2*i-1))
        sine += sine_term
        cosine += cosine_term
    return sine, cosine


def reference(raw):
    M, a, time, x, y, z, ut, ux, uy, uz = map(D.from_float,raw)
    spin_squared = a*a
    planar = x*x+y*y+z*z-spin_squared
    r_squared = (planar+(planar*planar+4*spin_squared*z*z).sqrt())/2
    r = r_squared.sqrt()
    pos, neg = M+(M*M-spin_squared).sqrt(), M-(M*M-spin_squared).sqrt()
    if r <= pos:
        return None
    # Partial fractions in dimensionful r, independent of the normalized
    # kernel decomposition. The time origin fixes the same constant.
    outer, inner = (r-pos)/M, (r-neg)/M
    if pos != neg:
        partial = (outer.ln()-inner.ln())/(pos-neg)
        time_shift = -4*M*(pos*outer.ln()-neg*inner.ln())/(pos-neg)
    else:
        partial = -1/(r-pos)
        time_shift = -4*M*(outer.ln()+M*partial)
    sine, cosine = trigonometry(-2*a*partial)
    # Rotation by -2 atan(a/r) from rational half-angle identities.
    base_sine,base_cosine=-2*a*r/(r*r+a*a),(r*r-a*a)/(r*r+a*a)
    sin_angle=sine*base_cosine+cosine*base_sine
    cos_angle=cosine*base_cosine-sine*base_sine
    X, Y = cos_angle*x-sin_angle*y, sin_angle*x+cos_angle*y

    # Differentiate the defining equation x^2/(r^2+a^2) +
    # y^2/(r^2+a^2)+z^2/r^2=1 implicitly.
    dr = ((x*ux+y*uy)/(r*r+a*a)+z*uz/(r*r))/(
        (x*x+y*y)*r/(r*r+a*a)**2+z*z/(r*r*r))
    time_rate=-4*M*r*dr/((r-pos)*(r-neg))
    rotation_rate=a*time_rate/(r*r+a*a)
    return (-time-time_shift, X,-Y,z, -ut-time_rate,
            cos_angle*ux-sin_angle*uy-Y*rotation_rate,
            -sin_angle*ux-cos_angle*uy-X*rotation_rate,uz)


def inputs(bits):
    samples = []
    for mass in (D('0.1'), D('1'), D('12.5')):
        for spin in (D(0), D('0.9'), D('-0.998'), D('0.999999'), D('1')):
            for offset in (D('0.0625'),D('0.5'),D(20)):
                # Make the spin exactly representable in this input precision
                # when defining the requested physical launch radius.
                form = '<'+EXTENT[bits]
                m = D.from_float(struct.unpack(form, struct.pack(form,float(mass)))[0])
                a = D.from_float(struct.unpack(form, struct.pack(form,float(m*spin)))[0])
                r = m+(m*m-a*a).sqrt()+m*offset
                z = r*D('0.3')
                x = ((r*r+a*a)*D('0.91')).sqrt()*D('0.6')
                y = ((r*r+a*a)*D('0.91')).sqrt()*D('0.8')
                samples.extend(map(float,(m,a,-m*2,x,y,z,D('-2.1'),D('0.19'),D('-0.43'),D('0.59'))))
    samples.extend((1,0,0,1,0,0,-1,0,0,0))
    (WORK/f'chart-{bits}.bin').write_bytes(struct.pack(f'<{len(samples)}{EXTENT[bits]}', *samples))


def validate(bits):
    source = WORK/f'chart-{bits}.bin'
    artifact = WORK/f'chart-{bits}.readback'
    samples = list(struct.iter_unpack('<10'+EXTENT[bits],source.read_bytes()))
    rows = list(struct.iter_unpack('<20'+EXTENT[bits],artifact.read_bytes()))
    assert len(rows)==len(samples)==46
    measure = {'forward':[D(0)]*8, 'round_trip':[D(0)]*8, 'capture':D(0)}
    root = (1-(1-4*(D('1.25')-D('0.75').sqrt())/3).sqrt())/2
    for incoming,result in zip(samples,rows):
        expected = reference(incoming)
        if expected is None:
            assert result==(0,)*20
            continue
        assert result[0]==result[9]==result[18]==1
        for key, measured, target in [('forward',result[1:9],expected),
                                      ('round_trip',result[10:18],map(D.from_float,incoming[2:]))]:
            measure[key] = [max(old,abs(D.from_float(got)-wanted)/max(1,abs(wanted)))
                            for old,got,wanted in zip(measure[key],measured,target)]
        measure['capture'] = max(measure['capture'], abs(D.from_float(result[19])-root))
    return {'rows':len(rows), 'input':sha(source), 'readback':sha(artifact),
            'maximum_normalized_errors':{key: ([str(x) for x in value] if isinstance(value,list)
                                               else str(value)) for key,value in measure.items()},
            'exact_capture_fraction':str(root)}


with localcontext() as context:
    context.prec = 75
    if sys.argv[1] == 'inputs':
        for rung in EXTENT:
            inputs(rung)
    elif sys.argv[1] == 'validate':
        findings = {rung:validate(rung) for rung in EXTENT}
        (WORK/'chart-reference.json').write_text(json.dumps(findings,sort_keys=True,indent=2)+'\n')
        print(json.dumps(findings,sort_keys=True,indent=2))
    else:
        raise ValueError('Expected inputs or validate')
