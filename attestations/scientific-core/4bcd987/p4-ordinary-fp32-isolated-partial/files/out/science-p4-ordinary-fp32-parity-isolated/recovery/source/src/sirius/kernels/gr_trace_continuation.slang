// Persistent per-active-pixel trace state. Physical integration fields retain
// the selected Real rung; source radiance retains the existing float arithmetic.
module gr_trace_continuation;
import gr_types;

public static const uint kTraceRunning = 1;
public static const uint kTraceTerminal = 2;
public static const uint kTraceFinalized = 3;
public static const uint kTraceCapture = 1;
public static const uint kTraceDisk = 2;
public static const uint kTraceEscape = 3;
public static const uint kTraceWorkLimit = 4;
public static const uint kTraceMinimumStep = 5;
public static const uint kTraceEventFailure = 6;
public static const uint kTraceInvalidInput = 7;
public static const uint kTraceIdentityMismatch = 8;
public static const uint kTraceArithmetic = 9;
public static const uint kTraceOppositeEscape = 10;

// Storage strides: 272 bytes fp32/fp32-comp; 480 bytes fp64.
public struct TraceContinuation {
    public Real4 x;
    public Real4 u;
    // Fixed canonical order: direction A/B, then position A/B.
    public Real4 X[4];
    public Real4 V[4];
    public Real4 compX;
    public Real4 compU;
    public Real4 integration; // adaptive step, affine length, last dt, seed
    public float4 source;    // colour RGB, transmittance
    public float4 volume;    // emission RGB, optical depth
    public uint4 control;    // phase, reason, attempts, accepts
    public uint4 identity;   // absolute x, absolute y, sample index, epoch
};
public Real4 TracePack(Vec4Cart v) { return Real4(v.t, v.x, v.y, v.z); }
public Vec4Cart TraceUnpack(Real4 v) { return MakeVec4Cart(v.x, v.y, v.z, v.w); }
public bool TraceStateFinite(TraceContinuation record) {
    for (int i = 0; i < 4; ++i)
        if (!all(isfinite(record.X[i])) || !all(isfinite(record.V[i]))) return false;
    return all(isfinite(record.x)) && all(isfinite(record.u)) &&
        all(isfinite(record.compX)) && all(isfinite(record.compU)) &&
        all(isfinite(record.integration)) && all(isfinite(record.source)) &&
        all(isfinite(record.volume));
}
public bool TraceReasonFinalizable(uint reason) {
    return reason == kTraceCapture || reason == kTraceDisk ||
        reason == kTraceEscape || reason == kTraceOppositeEscape;
}
