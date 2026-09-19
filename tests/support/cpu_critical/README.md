# CPU critical transport witnesses

`original_launches.json` preserves all fourteen saved CPU camera launches. The
legacy image coordinates, directions and observer boosts are unchanged. The
original replay uses Kerr M=1, a=0.998, an observer at r=50 and theta=60 degrees,
a 1920x1080 pinhole camera with 100-degree FOV, no disk, escape radius 200,
30,000 attempts, initial step 0.02, maximum step 0.25, minimum step 1e-6,
absolute/relative tolerances 1e-9, and a 0.002 maximum step inside radius 5.
These are input witnesses, not a claim that every launch has passed.

`inputs.json` preserves ten failed interval readbacks from the original replay
at source revision 05c7261. These include the actual incoming four covariant
variation columns. They are not independently projected reference initial
states; their constraint and numerical defects must remain visible.

`reference.py` defines independent Hamiltonian phase evolution using the
separate high-precision defining Kerr metric, a general matrix inverse,
high-precision numerical differentiation and RK4 refinement. It also constructs
the outgoing chart by the explicit time/azimuth reflection. The frozen metric
fixtures cover the metric, inverse and all first derivatives at all ten events.
They compare 75- and 110-digit evaluations; the observed gap is a numerical
precision witness, not an interval certificate.

Run `bin/toolchains/python/bin/python tests/support/cpu_critical/reference.py --metric-fixtures` from Sirius to explicitly regenerate `metric_reference.json`
and its C++ header. This requires mpmath. Ordinary builds consume the frozen
header and never regenerate their oracle. The three Mandatory
`RetainedArithmetic` tests check exact small-term arithmetic, an independent
rational gradient and all 960 metric fields, including controls that delete the
retained low parts. These tests do not qualify complete ray transport or images.

`projected_inputs.json` is a separate derived set: the original readbacks are
independently projected onto the null constraint, including their derivatives,
at 75 digits and then materialized once as binary64 inputs. The original files
remain unchanged. `transport_reference.json` and its C++ header freeze an
independent RK4 step and two half steps from each derived input. Their observed
refinement gaps are below 1e-18. Regenerate these files explicitly with
`reference.py --transport-fixtures`; ordinary builds do not run this generator.
`CoupledTransport.MinimumStepCriticalColumnsMatchIndependentRefinedFlow` admits
all ten derived inputs at their original 1e-6 minimum step and checks the central
ray and all four columns against those independent results using the original
central and variation budgets. It is a local transport check, not full-image
qualification.

The original launch header is decoded from its saved JSON with `reference.py
--launch-header`; that operation does not recompute camera samples. The full
`OriginalCriticalLaunchesHavePhysicalFatesOrExplicitWorkExhaustion` Mandatory
test replays all fourteen packets, checks each generated direction exactly,
requires physical outcomes for twelve rays, and accepts only explicit original
work-budget exhaustion for the two longer rays when they do not finish. It
does not increase their work limit or turn a failure into a dark pixel.
