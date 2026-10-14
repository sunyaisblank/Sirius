// Bounded vacuum-Kerr exterior continuation to infinity. This module declares
// no descriptors and performs no hidden whole-tail loop: the owner schedules
// one adaptive attempt at a time and must govern all subsequent dispatches.
// Carter equations: Gralla & Lupsasca (2020), arXiv:1910.12881, equations3--9.
// The angular vectors use the public ingoing Cartesian KS orientation.
module gr_kerr_infinity;

import gr_types;
import gr_metrics;

public static const uint kKerrInfinityRunning = 0;
public static const uint kKerrInfinityComplete = 1;
public static const uint kKerrInfinityInvalidInput = 2;
public static const uint kKerrInfinityRadialUnresolved = 3;
public static const uint kKerrInfinityWorkLimit = 4;
public static const uint kKerrInfinityArithmetic = 5;

public struct KerrInfinityConfigGPU {
    public Real absoluteTolerance;
    public Real relativeTolerance;
    public Real inputConstraintTolerance;
    public uint maximumAttempts;
}

// Each Real3 is one first-order jet: x=nominal, y/z=the two variations
// per unit launch angle. All arithmetic, including derivatives, follows Real.
public struct KerrInfinityStateGPU {
    public Real3 values[6]; // angular n then canonical angular tangent v
    public Real3 compensation[6];
    public Real3 inverseRadius;
    public Real3 energy;
    public Real3 angularMomentum;
    public Real3 carterCombination; // Q+(L-aE)^2
    public Real3 frequency;
    public Real mass;
    public Real spin;
    public Real fraction;
    public Real step;
    public Real maximumLocalErrorRatio;
    public uint acceptedSteps;
    public uint attemptedSteps;
    public uint status;
    public KerrInfinityConfigGPU config;
}

public struct KerrInfinityResultGPU {
    public Real3 direction;
    public vector<Real, 2> jacobian[2]; // rows, columns are the launch angles
    public Real determinant;
    public Real frequency;
    public vector<Real, 2> frequencyDerivative;
    public Real maximumLocalErrorRatio;
    public uint acceptedSteps;
    public uint attemptedSteps;
}

Real3 KIConstant(Real value) { return Real3(value, Real(0), Real(0)); }
bool KIFinite(Real3 value) { return all(isfinite(value)); }
Real3 KIMultiply(Real3 a, Real3 b) {
    return Real3(a.x * b.x, a.y * b.x + a.x * b.y, a.z * b.x + a.x * b.z);
}
Real3 KIDivide(Real3 a, Real3 b) {
    Real value = a.x / b.x;
    return Real3(value, (a.y - value * b.y) / b.x, (a.z - value * b.z) / b.x);
}
Real3 KISqrt(Real3 a) {
    Real value = sqrt(a.x);
    return Real3(value, a.y / (Real(2) * value), a.z / (Real(2) * value));
}
Real3 KIDot(Real3 a[3], Real3 b[3]) {
    Real3 sum = KIConstant(Real(0));
    for (uint i = 0; i < 3; ++i) sum += KIMultiply(a[i], b[i]);
    return sum;
}
Real KIEpsilon() {
#ifdef SIRIUS_FP64
    return Real(2.2204460492503131e-16);
#else
    return Real(1.1920928955078125e-7);
#endif
}
Real4 KIComponents(Vec4Cart value) { return Real4(value.t, value.x, value.y, value.z); }

// Unlike the CPU critical-point test, the device certifies positivity through
// the convex hull of degree-four Bernstein coefficients. Unresolved is an
// explicit refusal, not an assertion that a physical turning point was found.
// At most256 nodes and depth12; the usual far exterior interval needs one.
bool KICertifyRadialInterval(KerrInfinityStateGPU state) {
    Real u = state.inverseRadius.x, e = state.energy.x;
    Real q = state.carterCombination.x, a = state.spin;
    Real z = a * a * e - a * state.angularMomentum.x;
    Real u2 = u * u;
    Real power[5] = {e * e, Real(0), (Real(2) * e * z - q) * u2,
                     Real(2) * state.mass * q * u2 * u,
                     (z * z - a * a * q) * u2 * u2};
    Real scale = Real(0);
    for (uint i = 0; i < 5; ++i) {
        if (!isfinite(power[i])) return false;
        scale += abs(power[i]);
    }
    if (!isfinite(scale) || !(scale > Real(0))) return false;
    Real margin = Real(128) * KIEpsilon() * scale;
    Real pending[13][5];
    uint depths[13];
    pending[0][0] = power[0];
    pending[0][1] = power[0];
    pending[0][2] = power[0] + power[2] / Real(6);
    pending[0][3] = power[0] + power[2] / Real(2) + power[3] / Real(4);
    pending[0][4] = power[0] + power[2] + power[3] + power[4];
    depths[0] = 0;
    uint count = 1;
    for (uint visited = 0; visited < 256; ++visited) {
        if (count == 0) return true;
        --count;
        uint depth = depths[count];
        Real coefficients[5];
        bool positive = true;
        for (uint i = 0; i < 5; ++i) {
            coefficients[i] = pending[count][i];
            if (!isfinite(coefficients[i])) return false;
            positive = positive && coefficients[i] > margin;
        }
        if (positive) continue;
        if (coefficients[0] <= margin || coefficients[4] <= margin || depth >= 12 || count + 2 > 13)
            return false;
        Real left[5], right[5];
        left[0] = coefficients[0];
        right[4] = coefficients[4];
        for (uint level = 1; level < 5; ++level) {
            for (uint i = 0; i < 5 - level; ++i)
                coefficients[i] = (coefficients[i] + coefficients[i + 1]) * Real(0.5);
            left[level] = coefficients[0];
            right[4 - level] = coefficients[4 - level];
        }
        for (uint i = 0; i < 5; ++i) {
            pending[count][i] = right[i];
            pending[count + 1][i] = left[i];
        }
        depths[count] = depth + 1;
        depths[count + 1] = depth + 1;
        count += 2;
    }
    return count == 0;
}

// X and V are geometric endpoint displacement and covariant tangent
// variation, already transformed into the public ingoing chart. They contain
// the launch angular seed. Metric compatibility provides an analytic JVP:
// K^mu=V^mu-Gamma^mu_ab k^a X^b,
// dp_mu=g_mu_nu V^nu+Gamma^sigma_mu_nu p_sigma X^nu.
// No production finite differences and no second metric-derivative authority.
public uint InitKerrInfinity(float mass, float spin, Vec4Cart position, Vec4Cart tangent,
                            Vec4Cart X[2], Vec4Cart V[2], Real angularSeed,
                            KerrInfinityConfigGPU config, out KerrInfinityStateGPU state) {
    state = (KerrInfinityStateGPU)0;
    state.status = kKerrInfinityInvalidInput;
    state.config = config;
    state.mass = Real(mass);
    state.spin = Real(spin);
    state.step = Real(0.125);
    Real4 x4 = KIComponents(position), k4 = KIComponents(tangent);
    if (!isfinite(mass) || mass < 0 || !isfinite(spin) ||
        (mass > 0 && abs(spin) > mass) || !all(isfinite(x4)) || !all(isfinite(k4)) ||
        !isfinite(angularSeed) || !(angularSeed > Real(0)) ||
        !isfinite(config.absoluteTolerance) || !(config.absoluteTolerance > Real(0)) ||
        !isfinite(config.relativeTolerance) || config.relativeTolerance < Real(0) ||
        !isfinite(config.inputConstraintTolerance) || config.inputConstraintTolerance < Real(0))
        return state.status;
    Real radialScale = max(max(abs(position.x), abs(position.y)), max(abs(position.z), abs(Real(spin))));
    if (!(radialScale > Real(0)) || !isfinite(radialScale)) return state.status;
    Real radiusValue = ComputeKerrRadiusCart(position.x, position.y, position.z, Real(spin));
    Real horizon = Real(mass) + sqrt(max(Real(0), Real(mass) * Real(mass) - Real(spin) * Real(spin)));
    if (!isfinite(radiusValue) || !(radiusValue > horizon)) return state.status;
    Real metric[4][4], inverse[4][4], connection[4][4][4];
    GetKerrSchildFamilyMetric(position, mass, spin, 0, 0, metric, inverse);
    GetKerrSchildFamilyChristoffel(position, mass, spin, 0, 0, connection);
    Real4 displacements[2], variations[2];
    for (uint column = 0; column < 2; ++column) {
        displacements[column] = KIComponents(X[column]) / angularSeed;
        variations[column] = KIComponents(V[column]) / angularSeed;
        if (!all(isfinite(displacements[column])) || !all(isfinite(variations[column]))) return state.status;
    }
    Real3 momentum[4], fourTangent[4];
    for (uint mu = 0; mu < 4; ++mu) {
        momentum[mu] = KIConstant(Real(0));
        fourTangent[mu] = KIConstant(k4[mu]);
        for (uint nu = 0; nu < 4; ++nu) {
            if (!isfinite(metric[mu][nu])) return state.status;
            momentum[mu].x += metric[mu][nu] * k4[nu];
        }
    }
    for (uint column = 0; column < 2; ++column) {
        for (uint mu = 0; mu < 4; ++mu) {
            Real coordinateVariation = variations[column][mu];
            Real covectorVariation = Real(0);
            for (uint nu = 0; nu < 4; ++nu) {
                covectorVariation += metric[mu][nu] * variations[column][nu];
                for (uint sigma = 0; sigma < 4; ++sigma) {
                    coordinateVariation -= connection[mu][sigma][nu] * k4[sigma] * displacements[column][nu];
                    covectorVariation += connection[sigma][mu][nu] * momentum[sigma].x * displacements[column][nu];
                }
            }
            fourTangent[mu][column + 1] = coordinateVariation;
            momentum[mu][column + 1] = covectorVariation;
        }
    }
    state.frequency = momentum[0];
    if (!KIFinite(state.frequency) || !(state.frequency.x > Real(0))) return state.status;
    Real3 nullResidual = KIConstant(Real(0)), nullScale = KIConstant(Real(0));
    for (uint mu = 0; mu < 4; ++mu) {
        if (!KIFinite(momentum[mu]) || !KIFinite(fourTangent[mu])) return state.status;
        Real3 term = KIMultiply(momentum[mu], fourTangent[mu]);
        nullResidual += term;
        nullScale.x += abs(term.x);
        for (uint column = 1; column < 3; ++column)
            nullScale[column] += abs(momentum[mu][column] * fourTangent[mu].x) +
                                 abs(momentum[mu].x * fourTangent[mu][column]);
    }
    if (!KIFinite(nullResidual) || !KIFinite(nullScale) || !(nullScale.x > Real(0))) return state.status;
    for (uint component = 0; component < 3; ++component)
        if (abs(nullResidual[component]) > config.inputConstraintTolerance * max(nullScale.x, nullScale[component]))
            return state.status;
    for (uint mu = 0; mu < 4; ++mu) {
        momentum[mu] = KIDivide(momentum[mu], state.frequency);
        fourTangent[mu] = KIDivide(fourTangent[mu], state.frequency);
    }
    Real3 x[3], scaled[3];
    for (uint i = 0; i < 3; ++i) {
        x[i] = Real3(x4[i + 1], displacements[0][i + 1], displacements[1][i + 1]);
        scaled[i] = x[i] / radialScale;
    }
    Real a = Real(spin), scaledA = a / radialScale;
    Real3 difference = KIDot(scaled, scaled) - KIConstant(scaledA * scaledA);
    Real3 discriminant = KISqrt(KIMultiply(difference, difference) + KIMultiply(scaled[2], scaled[2]) * (Real(4) * scaledA * scaledA));
    // Rationalize the oblate root on the negative-difference branch, as the
    // shared nominal radius authority does. This avoids cancellation near the
    // disk while differentiating the same positive root, without a radius floor.
    Real3 radiusSquared = difference.x >= Real(0) ? (difference + discriminant) * Real(0.5) :
        KIDivide(KIMultiply(scaled[2], scaled[2]) * (Real(2) * scaledA * scaledA), discriminant - difference);
    Real3 radius = KISqrt(radiusSquared) * radialScale;
    if (!KIFinite(radius) || !(radius.x > horizon)) return state.status;
    Real3 r2 = KIMultiply(radius, radius), denominator = r2 + KIConstant(a * a);
    Real3 n[3] = {KIDivide(KIMultiply(radius, x[0]) + x[1] * a, denominator),
                  KIDivide(KIMultiply(radius, x[1]) - x[0] * a, denominator), KIDivide(x[2], radius)};
    Real3 w[3] = {KIMultiply(radius, momentum[1]) + momentum[2] * a,
                  KIMultiply(radius, momentum[2]) - momentum[1] * a, KIMultiply(radius, momentum[3])};
    Real3 projection = KIDot(n, w), v[3];
    for (uint i = 0; i < 3; ++i) {
        v[i] = w[i] - KIMultiply(n[i], projection);
        state.values[i] = n[i];
        state.values[i + 3] = v[i];
        if (!KIFinite(n[i]) || !KIFinite(v[i])) return state.status;
    }
    Real3 radialRate = KIDivide(
        KIMultiply(radius, KIMultiply(x[0], fourTangent[1]) + KIMultiply(x[1], fourTangent[2])) +
        KIDivide(KIMultiply(denominator, KIMultiply(x[2], fourTangent[3])), radius),
        discriminant * (radialScale * radialScale));
    if (!KIFinite(radialRate) || !(radialRate.x > Real(0))) return state.status;
    state.inverseRadius = KIDivide(KIConstant(Real(1)), radius);
    state.energy = -momentum[0];
    state.angularMomentum = KIMultiply(n[0], v[1]) - KIMultiply(n[1], v[0]);
    state.carterCombination = KIDot(v, v) - KIMultiply(state.energy, state.angularMomentum) * (Real(2) * a) +
        KIMultiply(KIMultiply(state.energy, state.energy), KIConstant(Real(1)) - KIMultiply(n[2], n[2])) * (a * a);
    if (!KIFinite(state.inverseRadius) || !KIFinite(state.energy) || !KIFinite(state.angularMomentum) ||
        !KIFinite(state.carterCombination)) return state.status;
    state.status = KICertifyRadialInterval(state) ? kKerrInfinityRunning : kKerrInfinityRadialUnresolved;
    return state.status;
}

bool KIDerivative(KerrInfinityStateGPU state, Real fraction, Real3 values[6], out Real3 derivative[6]) {
    Real3 u = state.inverseRadius * (Real(1) - fraction), u2 = KIMultiply(u, u);
    Real a = state.spin, M = state.mass;
    Real3 D = KIConstant(Real(1)) - u * (Real(2) * M) + u2 * (a * a);
    Real3 P = state.energy + KIMultiply(state.energy * (a * a) - state.angularMomentum * a, u2);
    Real3 H = KIMultiply(P, P) - KIMultiply(KIMultiply(D, u2), state.carterCombination);
    if (!KIFinite(H) || !KIFinite(D) || !(H.x > Real(0)) || !(D.x > Real(0))) return false;
    Real3 root = KISqrt(H);
    Real3 drag = KIDivide(KIMultiply(state.energy, u) * (Real(2) * M * a) -
                          KIMultiply(state.angularMomentum, u2) * (a * a), D);
    Real3 twist = KIDivide(drag, root) + KIDivide(KIConstant(a), D);
    Real3 n[3] = {values[0], values[1], values[2]}, v[3] = {values[3], values[4], values[5]};
    Real3 crossN[3] = {-n[1], n[0], KIConstant(Real(0))};
    Real3 crossV[3] = {-v[1], v[0], KIConstant(Real(0))};
    Real3 v2 = KIDot(v, v);
    for (uint i = 0; i < 3; ++i) {
        Real3 force = KIMultiply(KIMultiply(KIMultiply(state.energy, state.energy), n[2]),
                         KIConstant(i == 2 ? Real(1) : Real(0)) - KIMultiply(n[2], n[i])) * (a * a);
        derivative[i] = KIMultiply(state.inverseRadius, KIDivide(v[i], root) + KIMultiply(twist, crossN[i]));
        derivative[i + 3] = KIMultiply(state.inverseRadius,
            KIDivide(force - KIMultiply(v2, n[i]), root) + KIMultiply(twist, crossV[i]));
        if (!KIFinite(derivative[i]) || !KIFinite(derivative[i + 3])) return false;
    }
    return true;
}

// Explicit coefficients preserve the DP5(4) tableau while avoiding the
// dynamic local-array lookup that lost trial increments on the tested Dozen
// path. Actual axis and independent Kerr continuation controls cover this
// representation; unroll annotations alone did not repair that execution.
Real KICoefficient(uint stage, uint previous) {
    if (stage == 1 && previous == 0) return Real(1) / Real(5);
    if (stage == 2 && previous == 0) return Real(3) / Real(40);
    if (stage == 2 && previous == 1) return Real(9) / Real(40);
    if (stage == 3 && previous == 0) return Real(44) / Real(45);
    if (stage == 3 && previous == 1) return Real(-56) / Real(15);
    if (stage == 3 && previous == 2) return Real(32) / Real(9);
    if (stage == 4 && previous == 0) return Real(19372) / Real(6561);
    if (stage == 4 && previous == 1) return Real(-25360) / Real(2187);
    if (stage == 4 && previous == 2) return Real(64448) / Real(6561);
    if (stage == 4 && previous == 3) return Real(-212) / Real(729);
    if (stage == 5 && previous == 0) return Real(9017) / Real(3168);
    if (stage == 5 && previous == 1) return Real(-355) / Real(33);
    if (stage == 5 && previous == 2) return Real(46732) / Real(5247);
    if (stage == 5 && previous == 3) return Real(49) / Real(176);
    if (stage == 5 && previous == 4) return Real(-5103) / Real(18656);
    if (stage == 6 && previous == 0) return Real(35) / Real(384);
    if (stage == 6 && previous == 2) return Real(500) / Real(1113);
    if (stage == 6 && previous == 3) return Real(125) / Real(192);
    if (stage == 6 && previous == 4) return Real(-2187) / Real(6784);
    if (stage == 6 && previous == 5) return Real(11) / Real(84);
    return Real(0);
}
Real3 KIAccumulate(Real3 value, inout Real3 compensation, Real3 increment) {
#ifdef SIRIUS_FP32_COMP
    Real3 corrected = increment - compensation;
    Real3 next = value + corrected;
    compensation = (next - value) - corrected;
    return next;
#else
    return value + increment;
#endif
}

// Solve factor^5*error=.9^5 on the clamped interval using Real arithmetic.
// This avoids a float-narrowed transcendental on the fp64 rung. Controller
// approximation only selects the next attempt; acceptance uses the Real error.
Real KINextStepFactor(Real error) {
    if (error == Real(0)) return Real(5);
    Real target = (Real(9)/Real(10)) * (Real(9)/Real(10)) * (Real(9)/Real(10)) *
                  (Real(9)/Real(10)) * (Real(9)/Real(10));
    Real lower = Real(1)/Real(10), upper = Real(5);
    for (uint iteration = 0; iteration < 48; ++iteration) {
        Real middle = (lower + upper) * Real(0.5);
        Real squared = middle * middle;
        if (squared * squared * middle * error > target) upper = middle;
        else lower = middle;
    }
    return (lower + upper) * Real(0.5);
}

// Exactly one attempted DP5(4) step, shared by nominal and both JVP columns.
// Rejection retains values, compensation and fraction. Every attempt counts
// toward the caller's explicit work budget; terminal status remains sticky.
public uint AttemptKerrInfinityStep(inout KerrInfinityStateGPU state) {
    if (state.status != kKerrInfinityRunning) return state.status;
    if (state.attemptedSteps >= state.config.maximumAttempts) {
        state.status = kKerrInfinityWorkLimit;
        return state.status;
    }
    ++state.attemptedSteps;
    Real step = min(state.step, Real(1) - state.fraction);
    if (!isfinite(step) || !(step > Real(0)) || !(state.fraction + step > state.fraction)) {
        state.status = kKerrInfinityArithmetic;
        return state.status;
    }
    const Real nodes[7] = {Real(0),Real(1)/Real(5),Real(3)/Real(10),Real(4)/Real(5),Real(8)/Real(9),Real(1),Real(1)};
    const Real fourth[7] = {Real(5179)/Real(57600),Real(0),Real(7571)/Real(16695),Real(393)/Real(640),
                            -Real(92097)/Real(339200),Real(187)/Real(2100),Real(1)/Real(40)};
    Real3 stages[7][6], fifth[6], candidateCompensation[6];
    for (uint stage = 0; stage < 7; ++stage) {
        Real3 trial[6];
        for (uint i = 0; i < 6; ++i) {
            Real3 increment = KIConstant(Real(0)), sumCompensation = KIConstant(Real(0));
            for (uint previous = 0; previous < stage; ++previous) {
                Real3 term = stages[previous][i] * (step * KICoefficient(stage, previous));
                increment = KIAccumulate(increment, sumCompensation, term);
            }
            Real3 compensation = state.compensation[i];
            trial[i] = KIAccumulate(state.values[i], compensation, increment);
            if (stage == 6) {
                fifth[i] = trial[i];
                candidateCompensation[i] = compensation;
            }
        }
        Real3 derivative[6];
        if (!KIDerivative(state, state.fraction + nodes[stage] * step, trial, derivative)) {
            state.status = kKerrInfinityArithmetic;
            return state.status;
        }
        for (uint i = 0; i < 6; ++i) stages[stage][i] = derivative[i];
    }
    Real error = Real(0);
    for (uint i = 0; i < 6; ++i) {
        Real3 increment = KIConstant(Real(0)), sumCompensation = KIConstant(Real(0));
        for (uint stage = 0; stage < 7; ++stage)
            increment = KIAccumulate(increment, sumCompensation, stages[stage][i] * (step * fourth[stage]));
        Real3 compensation = state.compensation[i];
        Real3 lower = KIAccumulate(state.values[i], compensation, increment);
        if (!KIFinite(lower) || !KIFinite(fifth[i])) {
            state.status = kKerrInfinityArithmetic;
            return state.status;
        }
        for (uint component = 0; component < 3; ++component) {
            Real denominator = state.config.absoluteTolerance + state.config.relativeTolerance *
                max(abs(state.values[i][component]), abs(fifth[i][component]));
            Real difference = abs(fifth[i][component] - lower[component]);
            Real ratio = difference / denominator;
            if (!isfinite(denominator) || !(denominator > Real(0)) || !isfinite(difference) || !isfinite(ratio)) {
                state.status = kKerrInfinityArithmetic;
                return state.status;
            }
            error = max(error, ratio);
        }
    }
    if (error <= Real(1)) {
        for (uint i = 0; i < 6; ++i) {
            state.values[i] = fifth[i];
            state.compensation[i] = candidateCompensation[i];
        }
        state.fraction += step;
        ++state.acceptedSteps;
        state.maximumLocalErrorRatio = max(state.maximumLocalErrorRatio, error);
        if (state.fraction >= Real(1)) state.status = kKerrInfinityComplete;
    }
    state.step = step * KINextStepFactor(error);
    if (!isfinite(state.step)) state.status = kKerrInfinityArithmetic;
    return state.status;
}

// No result is available from a partial, failed, or budget-exhausted state.
public bool FinishKerrInfinity(KerrInfinityStateGPU state, out KerrInfinityResultGPU result) {
    result = (KerrInfinityResultGPU)0;
    if (state.status != kKerrInfinityComplete || state.fraction != Real(1)) return false;
    KerrInfinityResultGPU candidate = (KerrInfinityResultGPU)0;
    Real3 n[3] = {state.values[0], state.values[1], state.values[2]};
    Real3 length = KISqrt(KIDot(n, n));
    if (!KIFinite(length) || !(length.x > Real(0))) return false;
    for (uint i = 0; i < 3; ++i) {
        n[i] = KIDivide(n[i], length);
        if (!KIFinite(n[i])) return false;
        candidate.direction[i] = n[i].x;
    }
    uint least = 0;
    for (uint i = 1; i < 3; ++i)
        if (abs(candidate.direction[i]) < abs(candidate.direction[least])) least = i;
    Real3 first = -candidate.direction[least] * candidate.direction;
    first[least] += Real(1);
    Real norm = sqrt(dot(first, first));
    if (!isfinite(norm) || !(norm > Real(0))) return false;
    first /= norm;
    Real3 second = cross(candidate.direction, first);
    for (uint column = 0; column < 2; ++column) {
        for (uint i = 0; i < 3; ++i) {
            candidate.jacobian[0][column] += first[i] * n[i][column + 1];
            candidate.jacobian[1][column] += second[i] * n[i][column + 1];
        }
        candidate.frequencyDerivative[column] = state.frequency[column + 1];
    }
    candidate.frequency = state.frequency.x;
    candidate.determinant = candidate.jacobian[0][0] * candidate.jacobian[1][1] -
                            candidate.jacobian[0][1] * candidate.jacobian[1][0];
    candidate.maximumLocalErrorRatio = state.maximumLocalErrorRatio;
    candidate.acceptedSteps = state.acceptedSteps;
    candidate.attemptedSteps = state.attemptedSteps;
    if (!all(isfinite(candidate.direction)) || !all(isfinite(candidate.jacobian[0])) ||
        !all(isfinite(candidate.jacobian[1])) || !isfinite(candidate.determinant) ||
        !isfinite(candidate.frequency) || !all(isfinite(candidate.frequencyDerivative))) return false;
    result = candidate;
    return true;
}
