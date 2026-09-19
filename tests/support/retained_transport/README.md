# Retained device step witnesses

These fixtures compare the bounded device Hamiltonian RK stage with independent
75/105-digit geometry, generic metric inversion and numerical differentiation.
Ten inputs are the separately derived minimum-step CPU states preserved in
`../cpu_critical/projected_inputs.json`. An eleventh applies the exact chart
reflection. Three nonradial weak-field Kerr cases exercise half-unit and unit
steps in both regular charts, with independently projected angular and pupil
columns. An analytic flat-space trajectory supplies the dyadic control.
These fifteen witnesses are not replacements for the original launch packets
or a full trace test.

Every fixture retains the exact binary32 input high/low/tail/radius/status words and
all 160 scientific outputs: fifth-order phase state, fourth-order phase state,
the fifth-order increment, and embedded error. The phase order is position,
covariant momentum, then four position/covector variations. The step radius
describes finite arithmetic; the difference between reference precisions is an
observed stability witness. Neither certifies a global trajectory enclosure.

Transport uses three expansion terms because the earlier two-term candidate
lost small covariant components after projection at these critical events.
The camera stage retains its separately tested two-term representation. The
transport multiplication keeps the leading products explicitly and bounds
rounded or omitted lower-order products in its radius. The device test checks
every radius and rejects deletion of either retained correction term.
The frozen oracle also supplies twofold binary64 values so that these checks
preserve the third float term on hosts whose `long double` is only binary64.

Regenerate with a Python interpreter containing mpmath:

```
python tests/support/retained_transport/reference.py
```

Normal builds consume `reference_cases.h` without evaluating the reference.
The device code does not import these fixtures or the reference implementation.

`endpoint_reference.py` independently differentiates the coordinate null
quadratic, then converts its coordinate variations to physical covariant
variations and retained phase continuation. Its seventeen cases include the
exact linear horizon and an ergoregion with two spatial corrections separated
below one binary32 spacing. Comparing these corrections in scalar float picks
the wrong physical derivative. The device retains the comparison and declines
roots whose differentiated constraint is unrepresented.

`dense_reference.py` uses the defining cubic Hermite basis and independent
coordinate differentiation for forty-eight fixed-affine and moving-plane
samples, including the nonlinear flat controls. It compares the device's
factored secant polynomial, including small covariant variations and exact
endpoint ownership. For an arrival event,
geodesic flow changes X by k times the arrival shift; the corresponding
acceleration and connection terms cancel in V. The device preserves this
covariant identity without subtracting large rounded terms.

Both references use 75/105-digit stability checks and frozen binary64 pairs.
They are interpolation/projection witnesses, not full trajectory admission or
event-root localization tests. Regenerate them explicitly with the same mpmath
interpreter used above; normal builds consume the checked-in headers.
