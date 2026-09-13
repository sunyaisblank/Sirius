# Retained device step witnesses

These fixtures compare the bounded device Hamiltonian RK stage with independent
75/105-digit geometry, generic metric inversion and numerical differentiation.
Ten inputs are the separately derived minimum-step CPU states preserved in
`../cpu_critical/projected_inputs.json`. An eleventh applies the exact chart
reflection and a twelfth is an analytic flat-space trajectory with dyadic data.
They are not replacements for the original launch packets or a full trace test.

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
