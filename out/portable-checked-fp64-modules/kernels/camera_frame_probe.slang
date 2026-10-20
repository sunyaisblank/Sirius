// Observe the camera authority used at launch and escape, without integration.
// Inputs: 4 family parameters, event[4], forward/down/right[3 each], beta[3],
// rest direction[3], padding. Camera/family parameters retain the float ABI.
// Outputs: event[4], observer/forward/up/right[4 each], past ray[4],
// labelled Sachs screens[4 each], completion marker, float catalogue basis[6].
// Frame exports retain Real; catalogue values widen only for observation.
import gr_types;
import gr_metrics;
import gr_camera;
import gr_deviation;

[[vk::binding(0, 0)]] StructuredBuffer<Real> inputs;
[[vk::binding(1, 0)]] RWStructuredBuffer<Real> outputs;

float3 CameraInput3(uint offset) {
    return float3(float(inputs[offset]), float(inputs[offset + 1]),
                  float(inputs[offset + 2]));
}

void ExportCameraVector(uint offset, Vec4Cart value) {
    outputs[offset] = value.t;
    outputs[offset + 1] = value.x;
    outputs[offset + 2] = value.y;
    outputs[offset + 3] = value.z;
}

[shader("compute")]
[numthreads(1, 1, 1)]
void ComputeMain(uint3 id : SV_DispatchThreadID) {
    if (id.x != 0 || id.y != 0 || id.z != 0) return;
    outputs[32] = Real(0);
    float M = float(inputs[0]), a = float(inputs[1]);
    float Q = float(inputs[2]), Lambda = float(inputs[3]);
    Vec4Cart position = MakeVec4Cart(inputs[4], inputs[5], inputs[6], inputs[7]);
    if (!IsKerrSchildFamilyEventRepresented(position, M, a, Q, Lambda)) return;
    Real metric[4][4], inverse[4][4];
    GetKerrSchildFamilyMetric(position, M, a, Q, Lambda, metric, inverse);
    CameraObserverFrame frame = BuildCameraObserverFrame(
        metric, inverse, CameraInput3(8), CameraInput3(11), CameraInput3(14),
        CameraInput3(17));
    float3 direction = CameraInput3(20);
    Vec4Cart first, second;
    SachsScreenBasis(frame, direction, first, second);
    ExportCameraVector(0, position);
    ExportCameraVector(4, frame.time);
    ExportCameraVector(8, frame.forward);
    ExportCameraVector(12, frame.up);
    ExportCameraVector(16, frame.right);
    ExportCameraVector(20, PastDirectedCameraRay(frame, direction));
    ExportCameraVector(24, first);
    ExportCameraVector(28, second);
    float3 catalogueFirst, catalogueSecond;
    MakeCelestialTangentBasis(direction, catalogueFirst, catalogueSecond);
    for (int i = 0; i < 3; ++i) {
        outputs[33 + i] = Real(catalogueFirst[i]);
        outputs[36 + i] = Real(catalogueSecond[i]);
    }
    outputs[32] = Real(1);
}
