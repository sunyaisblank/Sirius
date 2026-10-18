#include <metal_stdlib>
#include <metal_math>
#include <metal_texture>
using namespace metal;

#line 334 "/home/astra/.project/Sirius/out/trace-jacobi-shared-segment-rhs-screen/inputs/candidate/gr_deviation.slang"
void MakeCelestialTangentBasis_0(float3 inputDirection_0, float3 thread* first_0, float3 thread* second_0)
{
    float3 direction_0 = normalize(inputDirection_0);



    float _S1 = abs(inputDirection_0.x);

#line 340
    float _S2 = abs(inputDirection_0.y);

#line 340
    bool _S3;

#line 340
    if(_S1 <= _S2)
    {

#line 340
        _S3 = _S1 <= (abs(inputDirection_0.z));

#line 340
    }
    else
    {

#line 340
        _S3 = false;

#line 340
    }

#line 340
    float3 reference_0;

#line 340
    if(_S3)
    {

#line 340
        reference_0 = float3(1.0, 0.0, 0.0);

#line 340
    }
    else
    {

#line 342
        if(_S2 <= (abs(inputDirection_0.z)))
        {

#line 342
            reference_0 = float3(0.0, 1.0, 0.0);

#line 342
        }
        else
        {

#line 342
            reference_0 = float3(0.0, 0.0, 1.0);

#line 342
        }

#line 340
    }

#line 340
    float3 _S4 = reference_0 - float3(dot(reference_0, direction_0))  * direction_0;

#line 340
    float3 _S5 = float3((1.0 / sqrt(dot(_S4, _S4)))) ;

#line 349
    *first_0 = _S4 * _S5;

    *second_0 = cross(direction_0, reference_0) * _S5;
    return;
}


#line 12253 "hlsl.meta.slang"
float ldexp_0(float x_0, float exp_0)
{

#line 12261
    return x_0 * exp2(exp_0);
}


#line 150 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_disk.slang"
constant array<float, int(16)> PT_GL_WEIGHTS_0 = { 0.02715245820581913, 0.06225352361798286, 0.09515850991010666, 0.12462896853685379, 0.14959599077701569, 0.1691565215587616, 0.18260350823402405, 0.18945060670375824, 0.18945060670375824, 0.18260350823402405, 0.1691565215587616, 0.14959599077701569, 0.12462896853685379, 0.09515850991010666, 0.06225352361798286, 0.02715245820581913 };

#line 144
constant array<float, int(16)> PT_GL_NODES_0 = { -0.98940092325210571, -0.94457501173019409, -0.86563122272491455, -0.75540441274642944, -0.61787623167037964, -0.4580167829990387, -0.28160354495048523, -0.09501250833272934, 0.09501250833272934, 0.28160354495048523, 0.4580167829990387, 0.61787623167037964, 0.75540441274642944, 0.86563122272491455, 0.94457501173019409, 0.98940092325210571 };

#line 232 "/home/astra/.project/Sirius/out/trace-jacobi-shared-segment-rhs-screen/inputs/candidate/trace.slang"
constant array<float3, int(12)> gradients_0 = { float3(1.0, 1.0, 0.0), float3(-1.0, 1.0, 0.0), float3(1.0, -1.0, 0.0), float3(-1.0, -1.0, 0.0), float3(1.0, 0.0, 1.0), float3(-1.0, 0.0, 1.0), float3(1.0, 0.0, -1.0), float3(-1.0, 0.0, -1.0), float3(0.0, 1.0, 1.0), float3(0.0, -1.0, 1.0), float3(0.0, 1.0, -1.0), float3(0.0, -1.0, -1.0) };

#line 421
struct _Array_natural_vectorx3Cfloatx2C4x3E4_0
{
    array<packed_float4, int(4)> data_0;
};


#line 421
struct TraceContinuation_natural_0
{
    packed_float4 x_1;
    packed_float4 u_0;
    _Array_natural_vectorx3Cfloatx2C4x3E4_0 X_0;
    _Array_natural_vectorx3Cfloatx2C4x3E4_0 V_0;
    packed_float4 compX_0;
    packed_float4 compU_0;
    packed_float4 integration_0;
    packed_float4 source_0;
    packed_float4 volume_0;
    packed_uint4 control_0;
    packed_uint4 identity_0;
};


#line 140
struct KernelContext_0
{
    float device* params_0;
    TraceContinuation_natural_0 device* continuation_0;
    float device* radiance_0;
    uint device* pointOffsets_0;
    uint device* pointIndices_0;
    float device* pointStars_0;
    uint device* starfield_0;
};


#line 65
float3 P3_0(uint base_0, KernelContext_0 thread* kernelContext_0)
{

#line 65
    return float3(kernelContext_0->params_0[base_0], kernelContext_0->params_0[base_0 + 1U], kernelContext_0->params_0[base_0 + 2U]);
}


#line 45 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_trace_continuation.slang"
bool TraceReasonFinalizable_0(uint reason_0)
{

#line 45
    bool _S6;
    if(reason_0 == 1U)
    {

#line 46
        _S6 = true;

#line 46
    }
    else
    {

#line 46
        _S6 = reason_0 == 2U;

#line 46
    }

#line 46
    if(_S6)
    {

#line 46
        _S6 = true;

#line 46
    }
    else
    {

#line 46
        _S6 = reason_0 == 3U;

#line 46
    }
    if(_S6)
    {

#line 47
        _S6 = true;

#line 47
    }
    else
    {

#line 47
        _S6 = reason_0 == 10U;

#line 47
    }

#line 46
    return _S6;
}


#line 48 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_types.slang"
struct Vec4Cart_0
{
    float t_0;
    float x_2;
    float y_0;
    float z_0;
};


#line 61
Vec4Cart_0 MakeVec4Cart_0(float t_1, float x_3, float y_1, float z_1)
{

#line 62
    thread Vec4Cart_0 v_0;
    (&v_0)->t_0 = t_1;

#line 63
    (&v_0)->x_2 = x_3;

#line 63
    (&v_0)->y_0 = y_1;

#line 63
    (&v_0)->z_0 = z_1;
    return v_0;
}


#line 10 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_camera.slang"
struct ThinLensProjectionSample_0
{
    float pupilRight_0;
    float pupilUp_0;
    float3 direction_1;
};


#line 16
ThinLensProjectionSample_0 ProjectThinLensSample_0(float imageRight_0, float imageUp_0, float tanHalfFov_0, float focalLength_0, float aperture_0, float focusDistance_0, float sampleU_0, float sampleV_0)
{

    float focusRight_0 = imageRight_0 * tanHalfFov_0 * focusDistance_0;
    float focusUp_0 = imageUp_0 * tanHalfFov_0 * focusDistance_0;



    float pupilAngle_0 = 6.28318548202514648 * sampleU_0;
    float pupilRadius_0 = focalLength_0 / 50.0 / (2.0 * aperture_0) * sqrt(sampleV_0);

    thread ThinLensProjectionSample_0 sample_0;
    float _S7 = pupilRadius_0 * cos(pupilAngle_0);

#line 28
    (&sample_0)->pupilRight_0 = _S7;
    float _S8 = pupilRadius_0 * sin(pupilAngle_0);

#line 29
    (&sample_0)->pupilUp_0 = _S8;
    (&sample_0)->direction_1 = normalize(float3(focusDistance_0, focusUp_0 - _S8, focusRight_0 - _S7));



    return sample_0;
}


#line 440 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_metrics.slang"
bool IsWarpDriveEventRepresented_0(const Vec4Cart_0 thread* pos_0, float vs_0, float sigma_0, float R_0, float xs_0, float ys_0, float zs_0)
{

#line 441
    bool _S9;

    if(!isfinite(pos_0->t_0))
    {

#line 443
        _S9 = true;

#line 443
    }
    else
    {

#line 443
        _S9 = !isfinite(pos_0->x_2);

#line 443
    }

#line 443
    if(_S9)
    {

#line 443
        _S9 = true;

#line 443
    }
    else
    {

#line 443
        _S9 = !isfinite(pos_0->y_0);

#line 443
    }

#line 443
    if(_S9)
    {

#line 443
        _S9 = true;

#line 443
    }
    else
    {

#line 443
        _S9 = !isfinite(pos_0->z_0);

#line 443
    }

#line 443
    if(_S9)
    {

#line 443
        _S9 = true;

#line 443
    }
    else
    {

#line 443
        _S9 = !isfinite(vs_0);

#line 443
    }
    if(_S9)
    {

#line 444
        _S9 = true;

#line 444
    }
    else
    {

#line 444
        _S9 = (abs(vs_0)) > 10.0;

#line 444
    }

#line 444
    if(_S9)
    {

#line 444
        _S9 = true;

#line 444
    }
    else
    {

#line 444
        _S9 = !isfinite(sigma_0);

#line 444
    }

#line 444
    if(_S9)
    {

#line 444
        _S9 = true;

#line 444
    }
    else
    {

#line 444
        _S9 = !(sigma_0 > 0.0);

#line 444
    }

#line 444
    if(_S9)
    {

#line 444
        _S9 = true;

#line 444
    }
    else
    {

#line 444
        _S9 = sigma_0 > 1000.0;

#line 444
    }
    if(_S9)
    {

#line 445
        _S9 = true;

#line 445
    }
    else
    {

#line 445
        _S9 = !isfinite(R_0);

#line 445
    }

#line 445
    if(_S9)
    {

#line 445
        _S9 = true;

#line 445
    }
    else
    {

#line 445
        _S9 = !(R_0 > 0.0);

#line 445
    }

#line 445
    if(_S9)
    {

#line 445
        _S9 = true;

#line 445
    }
    else
    {

#line 445
        _S9 = R_0 > 1000.0;

#line 445
    }

#line 445
    if(_S9)
    {

#line 445
        _S9 = true;

#line 445
    }
    else
    {

#line 445
        _S9 = !isfinite(xs_0);

#line 445
    }
    if(_S9)
    {

#line 446
        _S9 = true;

#line 446
    }
    else
    {

#line 446
        _S9 = !isfinite(ys_0);

#line 446
    }

#line 446
    if(_S9)
    {

#line 446
        _S9 = true;

#line 446
    }
    else
    {

#line 446
        _S9 = !isfinite(zs_0);

#line 446
    }

#line 443
    if(_S9)
    {


        return false;
    }
    float sigmaR_0 = sigma_0 * R_0;
    if(sigmaR_0 >= 0.10000000149011612)
    {

#line 450
        _S9 = sigmaR_0 <= 100.0;

#line 450
    }
    else
    {

#line 450
        _S9 = false;

#line 450
    }

#line 450
    return _S9;
}

float AlcubierreShapeFunction_0(float rs_0, float sigma_1, float R_1)
{


    return (tanh(sigma_1 * (rs_0 + R_1)) - tanh(sigma_1 * (rs_0 - R_1))) / max(2.0 * tanh(sigma_1 * R_1), 1.00000001335143196e-10);
}


#line 469
void GetWarpDriveMetric_0(const Vec4Cart_0 thread* pos_1, float vs_1, float sigma_2, float R_2, float xs_1, float ys_1, float zs_1, array<array<float, int(4)>, int(4)> thread* g_0, array<array<float, int(4)>, int(4)> thread* g_inv_0)
{

    int i_0 = int(0);

    for(;;)
    {

#line 474
        if(i_0 < int(4))
        {
        }
        else
        {

#line 474
            break;
        }

#line 474
        int j_0 = int(0);
        for(;;)
        {

#line 475
            if(j_0 < int(4))
            {
            }
            else
            {

#line 475
                break;
            }

#line 475
            (*g_0)[i_0][j_0] = 0.0;

#line 475
            (*g_inv_0)[i_0][j_0] = 0.0;

#line 475
            j_0 = j_0 + int(1);

#line 475
        }

#line 474
        i_0 = i_0 + int(1);

#line 474
    }

#line 474
    bool _S10 = IsWarpDriveEventRepresented_0(pos_1, vs_1, sigma_2, R_2, xs_1, ys_1, zs_1);

    if(!_S10)
    {

#line 476
        return;
    }

#line 485
    float dx_0 = pos_1->x_2 - (xs_1 + vs_1 * pos_1->t_0);
    float dy_0 = pos_1->y_0 - ys_1;
    float dz_0 = pos_1->z_0 - zs_1;

#line 495
    float vsf_0 = vs_1 * AlcubierreShapeFunction_0(max(sqrt(dx_0 * dx_0 + dy_0 * dy_0 + dz_0 * dz_0), 1.00000001335143196e-10), sigma_2, R_2);


    float _S11 = 1.0 - vsf_0 * vsf_0;

#line 498
    (*g_0)[int(0)][int(0)] = - _S11;
    float _S12 = - vsf_0;

#line 499
    (*g_0)[int(0)][int(1)] = _S12;
    (*g_0)[int(1)][int(0)] = _S12;
    (*g_0)[int(1)][int(1)] = 1.0;
    (*g_0)[int(2)][int(2)] = 1.0;
    (*g_0)[int(3)][int(3)] = 1.0;


    (*g_inv_0)[int(0)][int(0)] = -1.0;
    (*g_inv_0)[int(0)][int(1)] = _S12;
    (*g_inv_0)[int(1)][int(0)] = _S12;
    (*g_inv_0)[int(1)][int(1)] = _S11;
    (*g_inv_0)[int(2)][int(2)] = 1.0;
    (*g_inv_0)[int(3)][int(3)] = 1.0;
    return;
}


#line 363
bool IsMorrisThorneCartesianEventRepresented_0(const Vec4Cart_0 thread* pos_2, float b0_0, float Phi0_0, int shapeType_0)
{

#line 364
    bool _S13;

#line 370
    if(shapeType_0 != int(0))
    {

#line 370
        _S13 = true;

#line 370
    }
    else
    {

#line 370
        _S13 = !isfinite(b0_0);

#line 370
    }

#line 370
    if(_S13)
    {

#line 370
        _S13 = true;

#line 370
    }
    else
    {

#line 370
        _S13 = !(b0_0 > 0.0);

#line 370
    }

#line 370
    if(_S13)
    {

#line 370
        _S13 = true;

#line 370
    }
    else
    {

#line 370
        _S13 = Phi0_0 != 0.0;

#line 370
    }

#line 370
    if(_S13)
    {

#line 370
        _S13 = true;

#line 370
    }
    else
    {

#line 370
        _S13 = !isfinite(pos_2->t_0);

#line 370
    }
    if(_S13)
    {

#line 371
        _S13 = true;

#line 371
    }
    else
    {

#line 371
        _S13 = !isfinite(pos_2->x_2);

#line 371
    }

#line 371
    if(_S13)
    {

#line 371
        _S13 = true;

#line 371
    }
    else
    {

#line 371
        _S13 = !isfinite(pos_2->y_0);

#line 371
    }

#line 371
    if(_S13)
    {

#line 371
        _S13 = true;

#line 371
    }
    else
    {

#line 371
        _S13 = !isfinite(pos_2->z_0);

#line 371
    }

#line 370
    if(_S13)
    {
        return false;
    }
    float rho2_0 = pos_2->x_2 * pos_2->x_2 + pos_2->y_0 * pos_2->y_0 + pos_2->z_0 * pos_2->z_0;
    if(isfinite(rho2_0))
    {

#line 375
        _S13 = rho2_0 > 0.0;

#line 375
    }
    else
    {

#line 375
        _S13 = false;

#line 375
    }

#line 375
    return _S13;
}

void GetMorrisThorneCartesianMetric_0(const Vec4Cart_0 thread* pos_3, float b0_1, float Phi0_1, int shapeType_1, array<array<float, int(4)>, int(4)> thread* g_1, array<array<float, int(4)>, int(4)> thread* g_inv_1)
{
    int i_1 = int(0);

    for(;;)
    {

#line 382
        if(i_1 < int(4))
        {
        }
        else
        {

#line 382
            break;
        }

#line 382
        int j_1 = int(0);
        for(;;)
        {

#line 383
            if(j_1 < int(4))
            {
            }
            else
            {

#line 383
                break;
            }

#line 383
            (*g_1)[i_1][j_1] = 0.0;

#line 383
            (*g_inv_1)[i_1][j_1] = 0.0;

#line 383
            j_1 = j_1 + int(1);

#line 383
        }

#line 382
        i_1 = i_1 + int(1);

#line 382
    }

#line 382
    bool _S14 = IsMorrisThorneCartesianEventRepresented_0(pos_3, b0_1, Phi0_1, shapeType_1);

    if(!_S14)
    {

#line 384
        return;
    }

    float halfB0OverRho_0 = 0.5 * b0_1 / sqrt(pos_3->x_2 * pos_3->x_2 + pos_3->y_0 * pos_3->y_0 + pos_3->z_0 * pos_3->z_0);

    float conformalBase_0 = 1.0 + halfB0OverRho_0 * halfB0OverRho_0;
    float conformal_0 = conformalBase_0 * conformalBase_0;
    float _S15 = 1.0 / conformal_0;
    (*g_1)[int(0)][int(0)] = -1.0;
    (*g_inv_1)[int(0)][int(0)] = -1.0;

#line 393
    i_1 = int(0);
    for(;;)
    {

#line 394
        if(i_1 < int(3))
        {
        }
        else
        {

#line 394
            break;
        }

#line 395
        int _S16 = i_1 + int(1);

#line 395
        (*g_1)[_S16][_S16] = conformal_0;
        (*g_inv_1)[_S16][_S16] = _S15;

#line 394
        i_1 = _S16;

#line 394
    }



    return;
}


#line 55
float ComputeKerrRadiusCart_0(float x_4, float y_2, float z_2, float a_0)
{

#line 56
    float _S17 = max(max(abs(x_4), abs(y_2)), max(abs(z_2), abs(a_0)));
    if(_S17 == 0.0)
    {

#line 57
        return 0.0;
    }
    float scaledX_0 = x_4 / _S17;
    float scaledY_0 = y_2 / _S17;
    float scaledZ_0 = z_2 / _S17;
    float scaledA_0 = a_0 / _S17;

    float reduced_0 = scaledX_0 * scaledX_0 + scaledY_0 * scaledY_0 + scaledZ_0 * scaledZ_0 - scaledA_0 * scaledA_0;
    float scaledCross_0 = 2.0 * scaledA_0 * scaledZ_0;
    float discriminantRoot_0 = sqrt(reduced_0 * reduced_0 + scaledCross_0 * scaledCross_0);

#line 66
    float scaledRadius_0;

    if(reduced_0 >= 0.0)
    {

#line 68
        scaledRadius_0 = sqrt(0.5 * (reduced_0 + discriminantRoot_0));

#line 68
    }
    else
    {

#line 68
        bool _S18;

        if(scaledA_0 != 0.0)
        {

#line 70
            _S18 = scaledZ_0 != 0.0;

#line 70
        }
        else
        {

#line 70
            _S18 = false;

#line 70
        }

#line 70
        if(_S18)
        {

#line 70
            scaledRadius_0 = sqrt(2.0) * abs(scaledA_0 * scaledZ_0) / sqrt(discriminantRoot_0 - reduced_0);

#line 70
        }
        else
        {

#line 70
            scaledRadius_0 = 0.0;

#line 70
        }

#line 68
    }

#line 74
    return _S17 * scaledRadius_0;
}


#line 98
void GetKerrSchildFamilyMetric_0(const Vec4Cart_0 thread* pos_4, float M_0, float a_1, float Q_0, float Lambda_0, array<array<float, int(4)>, int(4)> thread* g_2, array<array<float, int(4)>, int(4)> thread* g_inv_2)
{
    int j_2;

#line 100
    int nu_0;

#line 100
    float _S19 = pos_4->x_2;

#line 100
    float _S20 = pos_4->y_0;

#line 100
    float _S21 = pos_4->z_0;

#line 108
    float a2_0 = a_1 * a_1;
    float Q2_0 = Q_0 * Q_0;

#line 109
    int i_2 = int(0);

    for(;;)
    {

#line 111
        if(i_2 < int(4))
        {
        }
        else
        {

#line 111
            break;
        }

#line 111
        j_2 = int(0);
        for(;;)
        {

#line 112
            if(j_2 < int(4))
            {
            }
            else
            {

#line 112
                break;
            }

#line 112
            (*g_2)[i_2][j_2] = 0.0;

#line 112
            (*g_inv_2)[i_2][j_2] = 0.0;

#line 112
            j_2 = j_2 + int(1);

#line 112
        }

#line 111
        i_2 = i_2 + int(1);

#line 111
    }


    (*g_2)[int(0)][int(0)] = -1.0;

#line 114
    (*g_inv_2)[int(0)][int(0)] = -1.0;
    (*g_2)[int(1)][int(1)] = 1.0;

#line 115
    (*g_inv_2)[int(1)][int(1)] = 1.0;
    (*g_2)[int(2)][int(2)] = 1.0;

#line 116
    (*g_inv_2)[int(2)][int(2)] = 1.0;
    (*g_2)[int(3)][int(3)] = 1.0;

#line 117
    (*g_inv_2)[int(3)][int(3)] = 1.0;

    bool _S22 = M_0 == 0.0;

#line 119
    bool _S23;

#line 119
    if(_S22)
    {

#line 119
        _S23 = a_1 == 0.0;

#line 119
    }
    else
    {

#line 119
        _S23 = false;

#line 119
    }

#line 119
    if(_S23)
    {

#line 119
        _S23 = Q_0 == 0.0;

#line 119
    }
    else
    {

#line 119
        _S23 = false;

#line 119
    }

#line 119
    if(_S23)
    {

#line 119
        _S23 = Lambda_0 == 0.0;

#line 119
    }
    else
    {

#line 119
        _S23 = false;

#line 119
    }

#line 119
    if(_S23)
    {

#line 120
        return;
    }

    float r_0 = ComputeKerrRadiusCart_0(_S19, _S20, _S21, a_1);
    float r2_0 = r_0 * r_0;

    if(_S22)
    {

#line 126
        _S23 = a_1 == 0.0;

#line 126
    }
    else
    {

#line 126
        _S23 = false;

#line 126
    }

#line 126
    if(_S23)
    {

#line 126
        _S23 = Q_0 == 0.0;

#line 126
    }
    else
    {

#line 126
        _S23 = false;

#line 126
    }

#line 126
    if(_S23)
    {

#line 126
        _S23 = r_0 == 0.0;

#line 126
    }
    else
    {

#line 126
        _S23 = false;

#line 126
    }

#line 126
    if(_S23)
    {

#line 127
        return;
    }


    if(M_0 != 0.0)
    {

#line 131
        _S23 = true;

#line 131
    }
    else
    {

#line 131
        _S23 = Q_0 != 0.0;

#line 131
    }

#line 131
    float H_0;

#line 131
    if(_S23)
    {

#line 132
        float cosine_0 = _S21 / r_0;

#line 132
        H_0 = (2.0 * M_0 * r_0 - Q2_0) / (r2_0 + a2_0 * cosine_0 * cosine_0);

#line 131
    }
    else
    {

#line 131
        H_0 = 0.0;

#line 131
    }

#line 136
    if(a_1 == 0.0)
    {

#line 136
        _S23 = Lambda_0 != 0.0;

#line 136
    }
    else
    {

#line 136
        _S23 = false;

#line 136
    }

#line 136
    if(_S23)
    {

#line 136
        H_0 = H_0 + Lambda_0 * r2_0 / 3.0;

#line 136
    }



    float denom_0 = r2_0 + a2_0;


    float _S24 = (r_0 * _S19 + a_1 * _S20) / denom_0;
    float _S25 = (r_0 * _S20 - a_1 * _S19) / denom_0;
    float _S26 = _S21 / r_0;

#line 141
    array<float, int(4)> l_0 = { 1.0, _S24, _S25, _S26 };

#line 141
    j_2 = int(0);

#line 148
    for(;;)
    {

#line 148
        if(j_2 < int(4))
        {
        }
        else
        {

#line 148
            break;
        }

#line 148
        nu_0 = int(0);
        for(;;)
        {

#line 149
            if(nu_0 < int(4))
            {
            }
            else
            {

#line 149
                break;
            }

#line 150
            (*g_2)[j_2][nu_0] = (*g_2)[j_2][nu_0] + H_0 * l_0[j_2] * l_0[nu_0];

#line 149
            nu_0 = nu_0 + int(1);

#line 149
        }

#line 148
        j_2 = j_2 + int(1);

#line 148
    }



    array<float, int(4)> _S27 = { 1.0, - _S24, - _S25, - _S26 };

#line 152
    j_2 = int(0);
    for(;;)
    {

#line 153
        if(j_2 < int(4))
        {
        }
        else
        {

#line 153
            break;
        }

#line 153
        nu_0 = int(0);
        for(;;)
        {

#line 154
            if(nu_0 < int(4))
            {
            }
            else
            {

#line 154
                break;
            }

#line 155
            (*g_inv_2)[j_2][nu_0] = (*g_inv_2)[j_2][nu_0] - H_0 * _S27[j_2] * _S27[nu_0];

#line 154
            nu_0 = nu_0 + int(1);

#line 154
        }

#line 153
        j_2 = j_2 + int(1);

#line 153
    }



    return;
}


#line 44 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_camera.slang"
float MetricDotCart_0(const array<array<float, int(4)>, int(4)> thread* metric_0, const Vec4Cart_0 thread* lhs_0, const Vec4Cart_0 thread* rhs_0)
{

#line 45
    array<float, int(4)> _S28 = { lhs_0->t_0, lhs_0->x_2, lhs_0->y_0, lhs_0->z_0 };
    array<float, int(4)> _S29 = { rhs_0->t_0, rhs_0->x_2, rhs_0->y_0, rhs_0->z_0 };

#line 46
    int mu_0 = int(0);

#line 46
    float result_0 = 0.0;

    for(;;)
    {

#line 48
        if(mu_0 < int(4))
        {
        }
        else
        {

#line 48
            break;
        }

#line 48
        int nu_1 = int(0);
        for(;;)
        {

#line 49
            if(nu_1 < int(4))
            {
            }
            else
            {

#line 49
                break;
            }

#line 50
            float result_1 = result_0 + (*metric_0)[mu_0][nu_1] * _S28[mu_0] * _S29[nu_1];

#line 49
            nu_1 = nu_1 + int(1);

#line 49
            result_0 = result_1;

#line 49
        }

#line 48
        mu_0 = mu_0 + int(1);

#line 48
    }

#line 53
    return result_0;
}


#line 44
float MetricDotCart_1(const array<array<float, int(4)>, int(4)> thread* metric_1, const Vec4Cart_0 thread* lhs_1, const Vec4Cart_0 thread* rhs_1)
{

#line 45
    array<float, int(4)> _S30 = { lhs_1->t_0, lhs_1->x_2, lhs_1->y_0, lhs_1->z_0 };
    array<float, int(4)> _S31 = { rhs_1->t_0, rhs_1->x_2, rhs_1->y_0, rhs_1->z_0 };

#line 46
    int mu_1 = int(0);

#line 46
    float result_2 = 0.0;

    for(;;)
    {

#line 48
        if(mu_1 < int(4))
        {
        }
        else
        {

#line 48
            break;
        }

#line 48
        int nu_2 = int(0);
        for(;;)
        {

#line 49
            if(nu_2 < int(4))
            {
            }
            else
            {

#line 49
                break;
            }

#line 50
            float result_3 = result_2 + (*metric_1)[mu_1][nu_2] * _S30[mu_1] * _S31[nu_2];

#line 49
            nu_2 = nu_2 + int(1);

#line 49
            result_2 = result_3;

#line 49
        }

#line 48
        mu_1 = mu_1 + int(1);

#line 48
    }

#line 53
    return result_2;
}


#line 79 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_types.slang"
Vec4Cart_0 ScaleVec4Cart_0(const Vec4Cart_0 thread* a_2, float s_0)
{

#line 80
    return MakeVec4Cart_0(a_2->t_0 * s_0, a_2->x_2 * s_0, a_2->y_0 * s_0, a_2->z_0 * s_0);
}


#line 79
Vec4Cart_0 ScaleVec4Cart_1(const Vec4Cart_0 thread* a_3, float s_1)
{

#line 80
    return MakeVec4Cart_0(a_3->t_0 * s_1, a_3->x_2 * s_1, a_3->y_0 * s_1, a_3->z_0 * s_1);
}


#line 75
Vec4Cart_0 AddVec4Cart_0(const Vec4Cart_0 thread* a_4, const Vec4Cart_0 thread* b_0)
{

#line 76
    return MakeVec4Cart_0(a_4->t_0 + b_0->t_0, a_4->x_2 + b_0->x_2, a_4->y_0 + b_0->y_0, a_4->z_0 + b_0->z_0);
}


#line 56 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_camera.slang"
Vec4Cart_0 NormaliseSpacelike_0(const array<array<float, int(4)>, int(4)> thread* metric_2, const Vec4Cart_0 thread* vector_0)
{

#line 56
    float _S32 = MetricDotCart_1(metric_2, vector_0, vector_0);

#line 56
    Vec4Cart_0 _S33 = ScaleVec4Cart_1(vector_0, 1.0 / sqrt(_S32));
    return _S33;
}


#line 37
struct CameraObserverFrame_0
{
    Vec4Cart_0 time_0;
    Vec4Cart_0 forward_0;
    Vec4Cart_0 up_0;
    Vec4Cart_0 right_0;
};


#line 63
CameraObserverFrame_0 BuildCameraObserverFrame_0(const array<array<float, int(4)>, int(4)> thread* metric_3, const array<array<float, int(4)>, int(4)> thread* inverseMetric_0, float3 camDir_0, float3 camUp_0, float3 camRight_0, float3 cameraBeta_0)
{


    thread CameraObserverFrame_0 reference_1;

    float _S34 = - (1.0 / sqrt(- (*inverseMetric_0)[int(0)][int(0)]));

#line 69
    Vec4Cart_0 _S35 = MakeVec4Cart_0(_S34 * (*inverseMetric_0)[int(0)][int(0)], _S34 * (*inverseMetric_0)[int(1)][int(0)], _S34 * (*inverseMetric_0)[int(2)][int(0)], _S34 * (*inverseMetric_0)[int(3)][int(0)]);

#line 69
    (&reference_1)->time_0 = _S35;

#line 74
    Vec4Cart_0 _S36 = MakeVec4Cart_0(0.0, camDir_0.x, camDir_0.y, camDir_0.z);

#line 74
    thread Vec4Cart_0 _S37 = _S36;

#line 74
    thread Vec4Cart_0 _S38 = _S35;

#line 74
    float _S39 = MetricDotCart_0(metric_3, &_S37, &_S38);

#line 74
    thread Vec4Cart_0 _S40 = _S35;

#line 74
    Vec4Cart_0 _S41 = ScaleVec4Cart_0(&_S40, _S39);

#line 74
    thread Vec4Cart_0 _S42 = _S36;

#line 74
    thread Vec4Cart_0 _S43 = _S41;

#line 74
    Vec4Cart_0 _S44 = AddVec4Cart_0(&_S42, &_S43);

#line 74
    thread Vec4Cart_0 _S45 = _S44;

#line 74
    Vec4Cart_0 _S46 = NormaliseSpacelike_0(metric_3, &_S45);

#line 79
    (&reference_1)->forward_0 = _S46;



    Vec4Cart_0 _S47 = MakeVec4Cart_0(0.0, - camUp_0.x, - camUp_0.y, - camUp_0.z);

#line 83
    thread Vec4Cart_0 _S48 = _S47;

#line 83
    thread Vec4Cart_0 _S49 = _S35;

#line 83
    float _S50 = MetricDotCart_0(metric_3, &_S48, &_S49);

#line 83
    thread Vec4Cart_0 _S51 = _S35;

#line 83
    Vec4Cart_0 _S52 = ScaleVec4Cart_0(&_S51, _S50);

#line 83
    thread Vec4Cart_0 _S53 = _S47;

#line 83
    thread Vec4Cart_0 _S54 = _S52;

#line 83
    Vec4Cart_0 _S55 = AddVec4Cart_0(&_S53, &_S54);

#line 83
    thread Vec4Cart_0 _S56 = _S55;

#line 83
    thread Vec4Cart_0 _S57 = _S46;

#line 83
    float _S58 = MetricDotCart_0(metric_3, &_S56, &_S57);

#line 90
    float _S59 = - _S58;

#line 90
    thread Vec4Cart_0 _S60 = _S46;

#line 90
    Vec4Cart_0 _S61 = ScaleVec4Cart_0(&_S60, _S59);

#line 90
    thread Vec4Cart_0 _S62 = _S55;

#line 90
    thread Vec4Cart_0 _S63 = _S61;

#line 90
    Vec4Cart_0 _S64 = AddVec4Cart_0(&_S62, &_S63);

#line 90
    thread Vec4Cart_0 _S65 = _S64;

#line 90
    Vec4Cart_0 _S66 = NormaliseSpacelike_0(metric_3, &_S65);
    (&reference_1)->up_0 = _S66;


    Vec4Cart_0 _S67 = MakeVec4Cart_0(0.0, camRight_0.x, camRight_0.y, camRight_0.z);

#line 94
    thread Vec4Cart_0 _S68 = _S67;

#line 94
    thread Vec4Cart_0 _S69 = _S35;

#line 94
    float _S70 = MetricDotCart_0(metric_3, &_S68, &_S69);

#line 94
    thread Vec4Cart_0 _S71 = _S35;

#line 94
    Vec4Cart_0 _S72 = ScaleVec4Cart_0(&_S71, _S70);

#line 94
    thread Vec4Cart_0 _S73 = _S67;

#line 94
    thread Vec4Cart_0 _S74 = _S72;

#line 94
    Vec4Cart_0 _S75 = AddVec4Cart_0(&_S73, &_S74);

#line 94
    thread Vec4Cart_0 _S76 = _S75;

#line 94
    thread Vec4Cart_0 _S77 = _S46;

#line 94
    float _S78 = MetricDotCart_0(metric_3, &_S76, &_S77);

#line 101
    float _S79 = - _S78;

#line 101
    thread Vec4Cart_0 _S80 = _S46;

#line 101
    Vec4Cart_0 _S81 = ScaleVec4Cart_0(&_S80, _S79);

#line 101
    thread Vec4Cart_0 _S82 = _S75;

#line 101
    thread Vec4Cart_0 _S83 = _S81;

#line 101
    Vec4Cart_0 _S84 = AddVec4Cart_0(&_S82, &_S83);

#line 101
    thread Vec4Cart_0 _S85 = _S84;

#line 101
    thread Vec4Cart_0 _S86 = _S66;

#line 101
    float _S87 = MetricDotCart_0(metric_3, &_S85, &_S86);


    float _S88 = - _S87;

#line 104
    thread Vec4Cart_0 _S89 = _S66;

#line 104
    Vec4Cart_0 _S90 = ScaleVec4Cart_0(&_S89, _S88);

#line 104
    thread Vec4Cart_0 _S91 = _S84;

#line 104
    thread Vec4Cart_0 _S92 = _S90;

#line 104
    Vec4Cart_0 _S93 = AddVec4Cart_0(&_S91, &_S92);

#line 104
    thread Vec4Cart_0 _S94 = _S93;

#line 104
    Vec4Cart_0 _S95 = NormaliseSpacelike_0(metric_3, &_S94);
    (&reference_1)->right_0 = _S95;


    float _S96 = cameraBeta_0.x;

#line 108
    float _S97 = cameraBeta_0.y;

#line 108
    float _S98 = cameraBeta_0.z;

#line 108
    float3 representedBeta_0 = float3(_S96, _S97, _S98);
    float betaSquared_0 = dot(representedBeta_0, representedBeta_0);
    if(betaSquared_0 == 0.0)
    {

#line 110
        return reference_1;
    }

#line 111
    float gamma_0 = 1.0 / sqrt(1.0 - betaSquared_0);

#line 111
    thread Vec4Cart_0 _S99 = (&reference_1)->forward_0;

#line 111
    Vec4Cart_0 _S100 = ScaleVec4Cart_0(&_S99, _S96);

#line 111
    thread Vec4Cart_0 _S101 = (&reference_1)->up_0;

#line 111
    Vec4Cart_0 _S102 = ScaleVec4Cart_0(&_S101, _S97);

#line 111
    thread Vec4Cart_0 _S103 = _S100;

#line 111
    thread Vec4Cart_0 _S104 = _S102;

#line 111
    Vec4Cart_0 _S105 = AddVec4Cart_0(&_S103, &_S104);

#line 111
    thread Vec4Cart_0 _S106 = (&reference_1)->right_0;

#line 111
    Vec4Cart_0 _S107 = ScaleVec4Cart_0(&_S106, _S98);

#line 111
    thread Vec4Cart_0 _S108 = _S105;

#line 111
    thread Vec4Cart_0 _S109 = _S107;

#line 111
    Vec4Cart_0 _S110 = AddVec4Cart_0(&_S108, &_S109);

#line 116
    thread CameraObserverFrame_0 boosted_0;

#line 116
    thread Vec4Cart_0 _S111 = (&reference_1)->time_0;

#line 116
    thread Vec4Cart_0 _S112 = _S110;

#line 116
    Vec4Cart_0 _S113 = AddVec4Cart_0(&_S111, &_S112);

#line 116
    thread Vec4Cart_0 _S114 = _S113;

#line 116
    Vec4Cart_0 _S115 = ScaleVec4Cart_0(&_S114, gamma_0);
    (&boosted_0)->time_0 = _S115;
    float coefficient_0 = (gamma_0 - 1.0) / betaSquared_0;


    float _S116 = coefficient_0 * _S96;

#line 121
    thread Vec4Cart_0 _S117 = _S110;

#line 121
    Vec4Cart_0 _S118 = ScaleVec4Cart_0(&_S117, _S116);

#line 121
    thread Vec4Cart_0 _S119 = (&reference_1)->forward_0;

#line 121
    thread Vec4Cart_0 _S120 = _S118;

#line 121
    Vec4Cart_0 _S121 = AddVec4Cart_0(&_S119, &_S120);
    float _S122 = gamma_0 * _S96;

#line 122
    thread Vec4Cart_0 _S123 = (&reference_1)->time_0;

#line 122
    Vec4Cart_0 _S124 = ScaleVec4Cart_0(&_S123, _S122);

#line 122
    thread Vec4Cart_0 _S125 = _S121;

#line 122
    thread Vec4Cart_0 _S126 = _S124;

#line 122
    Vec4Cart_0 _S127 = AddVec4Cart_0(&_S125, &_S126);

#line 119
    (&boosted_0)->forward_0 = _S127;

#line 125
    float _S128 = coefficient_0 * _S97;

#line 125
    thread Vec4Cart_0 _S129 = _S110;

#line 125
    Vec4Cart_0 _S130 = ScaleVec4Cart_0(&_S129, _S128);

#line 125
    thread Vec4Cart_0 _S131 = (&reference_1)->up_0;

#line 125
    thread Vec4Cart_0 _S132 = _S130;

#line 125
    Vec4Cart_0 _S133 = AddVec4Cart_0(&_S131, &_S132);
    float _S134 = gamma_0 * _S97;

#line 126
    thread Vec4Cart_0 _S135 = (&reference_1)->time_0;

#line 126
    Vec4Cart_0 _S136 = ScaleVec4Cart_0(&_S135, _S134);

#line 126
    thread Vec4Cart_0 _S137 = _S133;

#line 126
    thread Vec4Cart_0 _S138 = _S136;

#line 126
    Vec4Cart_0 _S139 = AddVec4Cart_0(&_S137, &_S138);

#line 123
    (&boosted_0)->up_0 = _S139;

#line 129
    float _S140 = coefficient_0 * _S98;

#line 129
    thread Vec4Cart_0 _S141 = _S110;

#line 129
    Vec4Cart_0 _S142 = ScaleVec4Cart_0(&_S141, _S140);

#line 129
    thread Vec4Cart_0 _S143 = (&reference_1)->right_0;

#line 129
    thread Vec4Cart_0 _S144 = _S142;

#line 129
    Vec4Cart_0 _S145 = AddVec4Cart_0(&_S143, &_S144);
    float _S146 = gamma_0 * _S98;

#line 130
    thread Vec4Cart_0 _S147 = (&reference_1)->time_0;

#line 130
    Vec4Cart_0 _S148 = ScaleVec4Cart_0(&_S147, _S146);

#line 130
    thread Vec4Cart_0 _S149 = _S145;

#line 130
    thread Vec4Cart_0 _S150 = _S148;

#line 130
    Vec4Cart_0 _S151 = AddVec4Cart_0(&_S149, &_S150);

#line 127
    (&boosted_0)->right_0 = _S151;



    return boosted_0;
}

Vec4Cart_0 PastDirectedCameraRay_0(const CameraObserverFrame_0 thread* frame_0, float3 screenDirection_0)
{

#line 135
    float3 direction_2 = float3(screenDirection_0.x, screenDirection_0.y, screenDirection_0.z);

    float3 direction_3 = direction_2 / float3(sqrt(dot(direction_2, direction_2))) ;

#line 137
    Vec4Cart_0 _S152 = ScaleVec4Cart_0(&frame_0->time_0, -1.0);

#line 137
    Vec4Cart_0 _S153 = ScaleVec4Cart_0(&frame_0->forward_0, direction_3.x);

#line 137
    thread Vec4Cart_0 _S154 = _S152;

#line 137
    thread Vec4Cart_0 _S155 = _S153;

#line 137
    Vec4Cart_0 _S156 = AddVec4Cart_0(&_S154, &_S155);

#line 137
    Vec4Cart_0 _S157 = ScaleVec4Cart_0(&frame_0->up_0, direction_3.y);

#line 137
    thread Vec4Cart_0 _S158 = _S156;

#line 137
    thread Vec4Cart_0 _S159 = _S157;

#line 137
    Vec4Cart_0 _S160 = AddVec4Cart_0(&_S158, &_S159);

#line 137
    Vec4Cart_0 _S161 = ScaleVec4Cart_0(&frame_0->right_0, direction_3.z);

#line 137
    thread Vec4Cart_0 _S162 = _S160;

#line 137
    thread Vec4Cart_0 _S163 = _S161;

#line 137
    Vec4Cart_0 _S164 = AddVec4Cart_0(&_S162, &_S163);
    return _S164;
}


#line 146
void SachsScreenBasis_0(const CameraObserverFrame_0 thread* frame_1, float3 restDirection_0, Vec4Cart_0 thread* first_1, Vec4Cart_0 thread* second_1)
{
    thread float3 firstComponents_0;
    thread float3 secondComponents_0;
    MakeCelestialTangentBasis_0(float3(restDirection_0.x, restDirection_0.y, restDirection_0.z), &firstComponents_0, &secondComponents_0);

#line 150
    Vec4Cart_0 _S165 = ScaleVec4Cart_0(&frame_1->forward_0, firstComponents_0.x);

#line 150
    Vec4Cart_0 _S166 = ScaleVec4Cart_0(&frame_1->up_0, firstComponents_0.y);

#line 150
    thread Vec4Cart_0 _S167 = _S165;

#line 150
    thread Vec4Cart_0 _S168 = _S166;

#line 150
    Vec4Cart_0 _S169 = AddVec4Cart_0(&_S167, &_S168);

#line 150
    Vec4Cart_0 _S170 = ScaleVec4Cart_0(&frame_1->right_0, firstComponents_0.z);

#line 150
    thread Vec4Cart_0 _S171 = _S169;

#line 150
    thread Vec4Cart_0 _S172 = _S170;

#line 150
    Vec4Cart_0 _S173 = AddVec4Cart_0(&_S171, &_S172);

    *first_1 = _S173;

#line 152
    Vec4Cart_0 _S174 = ScaleVec4Cart_0(&frame_1->forward_0, secondComponents_0.x);

#line 152
    Vec4Cart_0 _S175 = ScaleVec4Cart_0(&frame_1->up_0, secondComponents_0.y);

#line 152
    thread Vec4Cart_0 _S176 = _S174;

#line 152
    thread Vec4Cart_0 _S177 = _S175;

#line 152
    Vec4Cart_0 _S178 = AddVec4Cart_0(&_S176, &_S177);

#line 152
    Vec4Cart_0 _S179 = ScaleVec4Cart_0(&frame_1->right_0, secondComponents_0.z);

#line 152
    thread Vec4Cart_0 _S180 = _S178;

#line 152
    thread Vec4Cart_0 _S181 = _S179;

#line 152
    Vec4Cart_0 _S182 = AddVec4Cart_0(&_S180, &_S181);



    *second_1 = _S182;



    return;
}


#line 35 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_trace_continuation.slang"
float4 TracePack_0(const Vec4Cart_0 thread* v_1)
{

#line 35
    return float4(v_1->t_0, v_1->x_2, v_1->y_0, v_1->z_0);
}


#line 24 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_ellis_topology.slang"
float EllisOppositeEscapeRadius_0(float b0_2, float escapeRadius_0)
{

#line 25
    return b0_2 * b0_2 / (4.0 * escapeRadius_0);
}


#line 105 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_types.slang"
struct GeodesicStateCart_0
{
    Vec4Cart_0 x_5;
    Vec4Cart_0 u_1;
};


#line 14 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_ellis_topology.slang"
struct EllisTopologyTraceResult_0
{
    GeodesicStateCart_0 state_0;
    int termination_0;
    int attempts_0;
    float affineLength_0;
};


#line 55
struct EllisTopologyContinuation_0
{
    EllisTopologyTraceResult_0 trace_0;
    float adaptiveStep_0;
    Vec4Cart_0 compensationX_0;
    Vec4Cart_0 compensationU_0;
    uint accepted_0;
    float lastStep_0;
    uint failure_0;
};


#line 65
EllisTopologyContinuation_0 InitEllisTwoSheet_0(const GeodesicStateCart_0 thread* launch_0, int maxAttempts_0, float stepScale_0, float minStep_0, float maxStep_0, float escapeRadius_1, float b0_3, int wormShape_0)
{

    Vec4Cart_0 _S183 = { 0.0, 0.0, 0.0, 0.0 };

#line 68
    GeodesicStateCart_0 _S184 = { _S183, _S183 };

#line 68
    EllisTopologyTraceResult_0 _S185 = { _S184, int(0), int(0), 0.0 };

#line 68
    thread EllisTopologyContinuation_0 continuation_1;

#line 68
    (&continuation_1)->trace_0 = _S185;

#line 68
    (&continuation_1)->adaptiveStep_0 = 0.0;

#line 68
    (&continuation_1)->compensationX_0 = _S183;

#line 68
    (&continuation_1)->compensationU_0 = _S183;

#line 68
    (&continuation_1)->accepted_0 = 0U;

#line 68
    (&continuation_1)->lastStep_0 = 0.0;

#line 68
    (&continuation_1)->failure_0 = 0U;
    (&(&continuation_1)->trace_0)->state_0 = *launch_0;
    (&(&continuation_1)->trace_0)->termination_0 = int(0);
    (&continuation_1)->adaptiveStep_0 = maxStep_0;
    float throatRadius_0 = 0.5 * b0_3;
    float oppositeRadius_0 = EllisOppositeEscapeRadius_0(b0_3, escapeRadius_1);

#line 73
    bool _S186;
    if(maxAttempts_0 > int(0))
    {

#line 74
        _S186 = isfinite(stepScale_0);

#line 74
    }
    else
    {

#line 74
        _S186 = false;

#line 74
    }

#line 74
    if(_S186)
    {

#line 74
        _S186 = stepScale_0 > 0.0;

#line 74
    }
    else
    {

#line 74
        _S186 = false;

#line 74
    }

#line 74
    if(_S186)
    {

#line 74
        _S186 = isfinite(minStep_0);

#line 74
    }
    else
    {

#line 74
        _S186 = false;

#line 74
    }
    if(_S186)
    {

#line 75
        _S186 = minStep_0 > 0.0;

#line 75
    }
    else
    {

#line 75
        _S186 = false;

#line 75
    }

#line 75
    if(_S186)
    {

#line 75
        _S186 = isfinite(maxStep_0);

#line 75
    }
    else
    {

#line 75
        _S186 = false;

#line 75
    }

#line 75
    if(_S186)
    {

#line 75
        _S186 = maxStep_0 >= minStep_0;

#line 75
    }
    else
    {

#line 75
        _S186 = false;

#line 75
    }
    if(_S186)
    {

#line 76
        _S186 = isfinite(escapeRadius_1);

#line 76
    }
    else
    {

#line 76
        _S186 = false;

#line 76
    }

#line 76
    if(_S186)
    {

#line 76
        _S186 = escapeRadius_1 > throatRadius_0;

#line 76
    }
    else
    {

#line 76
        _S186 = false;

#line 76
    }
    if(_S186)
    {

#line 77
        _S186 = isfinite(b0_3);

#line 77
    }
    else
    {

#line 77
        _S186 = false;

#line 77
    }

#line 77
    if(_S186)
    {

#line 77
        _S186 = b0_3 > 0.0;

#line 77
    }
    else
    {

#line 77
        _S186 = false;

#line 77
    }

#line 77
    if(_S186)
    {

#line 77
        _S186 = isfinite(oppositeRadius_0);

#line 77
    }
    else
    {

#line 77
        _S186 = false;

#line 77
    }
    if(_S186)
    {

#line 78
        _S186 = oppositeRadius_0 > 0.0;

#line 78
    }
    else
    {

#line 78
        _S186 = false;

#line 78
    }

#line 78
    if(_S186)
    {

#line 78
        _S186 = oppositeRadius_0 < throatRadius_0;

#line 78
    }
    else
    {

#line 78
        _S186 = false;

#line 78
    }
    if(_S186)
    {

#line 79
        _S186 = wormShape_0 == int(0);

#line 79
    }
    else
    {

#line 79
        _S186 = false;

#line 79
    }

#line 74
    if(!_S186)
    {

#line 80
        (&continuation_1)->failure_0 = 1U;
        return continuation_1;
    }

    return continuation_1;
}


#line 21 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_trace_continuation.slang"
struct TraceContinuation_0
{
    float4 x_1;
    float4 u_0;
    array<float4, int(4)> X_0;
    array<float4, int(4)> V_0;
    float4 compX_0;
    float4 compU_0;
    float4 integration_0;
    float4 source_0;
    float4 volume_0;
    uint4 control_0;
    uint4 identity_0;
};


bool TraceStateFinite_0(const TraceContinuation_0 thread* record_0)
{

#line 37
    bool _S187;

#line 37
    int i_3 = int(0);
    for(;;)
    {

#line 38
        if(i_3 < int(4))
        {
        }
        else
        {

#line 38
            break;
        }

#line 39
        if(!all(isfinite(record_0->X_0[i_3])))
        {

#line 39
            _S187 = true;

#line 39
        }
        else
        {

#line 39
            _S187 = !all(isfinite(record_0->V_0[i_3]));

#line 39
        }

#line 39
        if(_S187)
        {

#line 39
            return false;
        }

#line 38
        i_3 = i_3 + int(1);

#line 38
    }

    if(all(isfinite(record_0->x_1)))
    {

#line 40
        _S187 = all(isfinite(record_0->u_0));

#line 40
    }
    else
    {

#line 40
        _S187 = false;

#line 40
    }

#line 40
    if(_S187)
    {

#line 40
        _S187 = all(isfinite(record_0->compX_0));

#line 40
    }
    else
    {

#line 40
        _S187 = false;

#line 40
    }
    if(_S187)
    {

#line 41
        _S187 = all(isfinite(record_0->compU_0));

#line 41
    }
    else
    {

#line 41
        _S187 = false;

#line 41
    }

#line 41
    if(_S187)
    {

#line 41
        _S187 = all(isfinite(record_0->integration_0));

#line 41
    }
    else
    {

#line 41
        _S187 = false;

#line 41
    }
    if(_S187)
    {

#line 42
        _S187 = all(isfinite(record_0->source_0));

#line 42
    }
    else
    {

#line 42
        _S187 = false;

#line 42
    }

#line 42
    if(_S187)
    {

#line 42
        _S187 = all(isfinite(record_0->volume_0));

#line 42
    }
    else
    {

#line 42
        _S187 = false;

#line 42
    }

#line 40
    return _S187;
}


#line 36
Vec4Cart_0 TraceUnpack_0(float4 v_2)
{

#line 36
    return MakeVec4Cart_0(v_2.x, v_2.y, v_2.z, v_2.w);
}


#line 417 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_integrator.slang"
float AdaptiveIntegrationTolerance_0()
{


    return 0.00001999999949476;
}


float AdaptiveNullTolerance_0()
{

#line 433
    return 0.00050000002374873;
}


#line 460 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_metrics.slang"
float AlcubierreShapeDerivative_0(float rs_1, float sigma_3, float R_3)
{

#line 461
    float cosh_plus_0 = cosh(sigma_3 * (rs_1 + R_3));
    float cosh_minus_0 = cosh(sigma_3 * (rs_1 - R_3));



    return sigma_3 * (1.0 / (cosh_plus_0 * cosh_plus_0) - 1.0 / (cosh_minus_0 * cosh_minus_0)) / max(2.0 * tanh(sigma_3 * R_3), 1.00000001335143196e-10);
}


#line 514
void GetWarpDriveChristoffel_0(const Vec4Cart_0 thread* pos_5, float vs_2, float sigma_4, float R_4, float xs_2, float ys_2, float zs_2, array<array<array<float, int(4)>, int(4)>, int(4)> thread* Gamma_0)
{

    int j_3;

#line 517
    int k_0;

#line 517
    int mu_2;

#line 517
    int nu_3;

#line 517
    int i_4 = int(0);

    for(;;)
    {

#line 519
        if(i_4 < int(4))
        {
        }
        else
        {

#line 519
            break;
        }

#line 519
        j_3 = int(0);
        for(;;)
        {

#line 520
            if(j_3 < int(4))
            {
            }
            else
            {

#line 520
                break;
            }

#line 520
            k_0 = int(0);
            for(;;)
            {

#line 521
                if(k_0 < int(4))
                {
                }
                else
                {

#line 521
                    break;
                }

#line 522
                (*Gamma_0)[i_4][j_3][k_0] = 0.0;

#line 521
                k_0 = k_0 + int(1);

#line 521
            }

#line 520
            j_3 = j_3 + int(1);

#line 520
        }

#line 519
        i_4 = i_4 + int(1);

#line 519
    }

#line 519
    bool _S188 = IsWarpDriveEventRepresented_0(pos_5, vs_2, sigma_4, R_4, xs_2, ys_2, zs_2);



    if(!_S188)
    {

#line 523
        return;
    }

#line 529
    float dx_1 = pos_5->x_2 - (xs_2 + vs_2 * pos_5->t_0);
    float dy_1 = pos_5->y_0 - ys_2;
    float dz_1 = pos_5->z_0 - zs_2;

    float _S189 = max(sqrt(dx_1 * dx_1 + dy_1 * dy_1 + dz_1 * dz_1), 1.00000001335143196e-10);


    float _S190 = AlcubierreShapeFunction_0(_S189, sigma_4, R_4);
    float df_drs_0 = AlcubierreShapeDerivative_0(_S189, sigma_4, R_4);

    float _S191 = - vs_2;
    array<float, int(4)> _S192 = { df_drs_0 * (_S191 * dx_1 / _S189), df_drs_0 * (dx_1 / _S189), df_drs_0 * (dy_1 / _S189), df_drs_0 * (dz_1 / _S189) };



    thread array<array<float, int(4)>, int(4)> metric_4;
    thread array<array<float, int(4)>, int(4)> inverseMetric_1;

#line 545
    GetWarpDriveMetric_0(pos_5, vs_2, sigma_4, R_4, xs_2, ys_2, zs_2, &metric_4, &inverseMetric_1);

    thread array<array<array<float, int(4)>, int(4)>, int(4)> dg_0;

#line 547
    j_3 = int(0);
    for(;;)
    {

#line 548
        if(j_3 < int(4))
        {
        }
        else
        {

#line 548
            break;
        }

#line 548
        mu_2 = int(0);
        for(;;)
        {

#line 549
            if(mu_2 < int(4))
            {
            }
            else
            {

#line 549
                break;
            }

#line 549
            nu_3 = int(0);
            for(;;)
            {

#line 550
                if(nu_3 < int(4))
                {
                }
                else
                {

#line 550
                    break;
                }

#line 551
                dg_0[j_3][mu_2][nu_3] = 0.0;

#line 550
                nu_3 = nu_3 + int(1);

#line 550
            }

#line 549
            mu_2 = mu_2 + int(1);

#line 549
        }

#line 554
        dg_0[j_3][int(0)][int(0)] = 2.0 * vs_2 * vs_2 * _S190 * _S192[j_3];
        float _S193 = _S191 * _S192[j_3];

#line 555
        dg_0[j_3][int(0)][int(1)] = _S193;
        dg_0[j_3][int(1)][int(0)] = _S193;

#line 548
        j_3 = j_3 + int(1);

#line 548
    }

#line 548
    k_0 = int(0);

#line 559
    for(;;)
    {

#line 559
        if(k_0 < int(4))
        {
        }
        else
        {

#line 559
            break;
        }

#line 559
        mu_2 = int(0);
        for(;;)
        {

#line 560
            if(mu_2 < int(4))
            {
            }
            else
            {

#line 560
                break;
            }

#line 560
            nu_3 = int(0);
            for(;;)
            {

#line 561
                if(nu_3 < int(4))
                {
                }
                else
                {

#line 561
                    break;
                }

#line 561
                int sigmaIndex_0 = int(0);

#line 561
                float sum_0 = 0.0;

                for(;;)
                {

#line 563
                    if(sigmaIndex_0 < int(4))
                    {
                    }
                    else
                    {

#line 563
                        break;
                    }

#line 564
                    float sum_1 = sum_0 + inverseMetric_1[k_0][sigmaIndex_0] * (dg_0[mu_2][sigmaIndex_0][nu_3] + dg_0[nu_3][sigmaIndex_0][mu_2] - dg_0[sigmaIndex_0][mu_2][nu_3]);

#line 563
                    sigmaIndex_0 = sigmaIndex_0 + int(1);

#line 563
                    sum_0 = sum_1;

#line 563
                }

#line 568
                (*Gamma_0)[k_0][mu_2][nu_3] = 0.5 * sum_0;

#line 561
                nu_3 = nu_3 + int(1);

#line 561
            }

#line 560
            mu_2 = mu_2 + int(1);

#line 560
        }

#line 559
        k_0 = k_0 + int(1);

#line 559
    }

#line 572
    return;
}


#line 170 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_integrator.slang"
bool TryGeodesicAccelerationWarpCart_0(const GeodesicStateCart_0 thread* state_1, float vs_3, float sigma_5, float R_5, float xs_3, float ys_3, float zs_3, Vec4Cart_0 thread* acceleration_0)
{

    *acceleration_0 = MakeVec4Cart_0(0.0, 0.0, 0.0, 0.0);

#line 173
    bool _S194 = IsWarpDriveEventRepresented_0(&state_1->x_5, vs_3, sigma_5, R_5, xs_3, ys_3, zs_3);
    if(!_S194)
    {

#line 174
        return false;
    }

#line 175
    thread array<array<array<float, int(4)>, int(4)>, int(4)> Gamma_1;

#line 175
    GetWarpDriveChristoffel_0(&state_1->x_5, vs_3, sigma_5, R_5, xs_3, ys_3, zs_3, &Gamma_1);

    array<float, int(4)> _S195 = { (&state_1->u_1)->t_0, (&state_1->u_1)->x_2, (&state_1->u_1)->y_0, (&state_1->u_1)->z_0 };
    thread array<float, int(4)> acc_0;

#line 178
    acc_0[int(0)] = 0.0;

#line 178
    acc_0[int(1)] = 0.0;

#line 178
    acc_0[int(2)] = 0.0;

#line 178
    acc_0[int(3)] = 0.0;

#line 178
    int mu_3 = int(0);
    for(;;)
    {

#line 179
        if(mu_3 < int(4))
        {
        }
        else
        {

#line 179
            break;
        }

#line 179
        int nu_4 = int(0);
        for(;;)
        {

#line 180
            if(nu_4 < int(4))
            {
            }
            else
            {

#line 180
                break;
            }

#line 180
            int rho_0 = int(0);
            for(;;)
            {

#line 181
                if(rho_0 < int(4))
                {
                }
                else
                {

#line 181
                    break;
                }

#line 182
                acc_0[mu_3] = acc_0[mu_3] - Gamma_1[mu_3][nu_4][rho_0] * _S195[nu_4] * _S195[rho_0];

#line 181
                rho_0 = rho_0 + int(1);

#line 181
            }

#line 180
            nu_4 = nu_4 + int(1);

#line 180
        }

#line 179
        mu_3 = mu_3 + int(1);

#line 179
    }



    *acceleration_0 = MakeVec4Cart_0(acc_0[int(0)], acc_0[int(1)], acc_0[int(2)], acc_0[int(3)]);

#line 183
    bool _S196;
    if(isfinite(acceleration_0->t_0))
    {

#line 184
        _S196 = isfinite(acceleration_0->x_2);

#line 184
    }
    else
    {

#line 184
        _S196 = false;

#line 184
    }

#line 184
    if(_S196)
    {

#line 184
        _S196 = isfinite(acceleration_0->y_0);

#line 184
    }
    else
    {

#line 184
        _S196 = false;

#line 184
    }
    if(_S196)
    {

#line 185
        _S196 = isfinite(acceleration_0->z_0);

#line 185
    }
    else
    {

#line 185
        _S196 = false;

#line 185
    }

#line 184
    return _S196;
}


bool IntegrateGeodesicRK4WarpCartDelta_0(const GeodesicStateCart_0 thread* state_2, float dt_0, float vs_4, float sigma_6, float R_6, float xs_4, float ys_4, float zs_4, Vec4Cart_0 thread* dx_2, Vec4Cart_0 thread* du_0)
{

    Vec4Cart_0 _S197 = MakeVec4Cart_0(0.0, 0.0, 0.0, 0.0);

#line 191
    *dx_2 = _S197;
    *du_0 = _S197;
    thread Vec4Cart_0 acceleration_1;

#line 193
    bool _S198 = TryGeodesicAccelerationWarpCart_0(state_2, vs_4, sigma_6, R_6, xs_4, ys_4, zs_4, &acceleration_1);
    if(!_S198)
    {

#line 195
        return false;
    }

#line 195
    Vec4Cart_0 _S199 = ScaleVec4Cart_0(&state_2->u_1, dt_0);

#line 195
    thread Vec4Cart_0 _S200 = acceleration_1;

#line 195
    Vec4Cart_0 _S201 = ScaleVec4Cart_0(&_S200, dt_0);



    thread GeodesicStateCart_0 s2_0;

#line 199
    thread Vec4Cart_0 _S202 = _S199;

#line 199
    Vec4Cart_0 _S203 = ScaleVec4Cart_0(&_S202, 0.5);

#line 199
    thread Vec4Cart_0 _S204 = _S203;

#line 199
    Vec4Cart_0 _S205 = AddVec4Cart_0(&state_2->x_5, &_S204);
    (&s2_0)->x_5 = _S205;

#line 200
    thread Vec4Cart_0 _S206 = _S201;

#line 200
    Vec4Cart_0 _S207 = ScaleVec4Cart_0(&_S206, 0.5);

#line 200
    thread Vec4Cart_0 _S208 = _S207;

#line 200
    Vec4Cart_0 _S209 = AddVec4Cart_0(&state_2->u_1, &_S208);
    (&s2_0)->u_1 = _S209;

#line 201
    thread GeodesicStateCart_0 _S210 = s2_0;

#line 201
    bool _S211 = TryGeodesicAccelerationWarpCart_0(&_S210, vs_4, sigma_6, R_6, xs_4, ys_4, zs_4, &acceleration_1);
    if(!_S211)
    {

#line 203
        return false;
    }

#line 203
    thread Vec4Cart_0 _S212 = (&s2_0)->u_1;

#line 203
    Vec4Cart_0 _S213 = ScaleVec4Cart_0(&_S212, dt_0);

#line 203
    thread Vec4Cart_0 _S214 = acceleration_1;

#line 203
    Vec4Cart_0 _S215 = ScaleVec4Cart_0(&_S214, dt_0);



    thread GeodesicStateCart_0 s3_0;

#line 207
    thread Vec4Cart_0 _S216 = _S213;

#line 207
    Vec4Cart_0 _S217 = ScaleVec4Cart_0(&_S216, 0.5);

#line 207
    thread Vec4Cart_0 _S218 = _S217;

#line 207
    Vec4Cart_0 _S219 = AddVec4Cart_0(&state_2->x_5, &_S218);
    (&s3_0)->x_5 = _S219;

#line 208
    thread Vec4Cart_0 _S220 = _S215;

#line 208
    Vec4Cart_0 _S221 = ScaleVec4Cart_0(&_S220, 0.5);

#line 208
    thread Vec4Cart_0 _S222 = _S221;

#line 208
    Vec4Cart_0 _S223 = AddVec4Cart_0(&state_2->u_1, &_S222);
    (&s3_0)->u_1 = _S223;

#line 209
    thread GeodesicStateCart_0 _S224 = s3_0;

#line 209
    bool _S225 = TryGeodesicAccelerationWarpCart_0(&_S224, vs_4, sigma_6, R_6, xs_4, ys_4, zs_4, &acceleration_1);
    if(!_S225)
    {

#line 211
        return false;
    }

#line 211
    thread Vec4Cart_0 _S226 = (&s3_0)->u_1;

#line 211
    Vec4Cart_0 _S227 = ScaleVec4Cart_0(&_S226, dt_0);

#line 211
    thread Vec4Cart_0 _S228 = acceleration_1;

#line 211
    Vec4Cart_0 _S229 = ScaleVec4Cart_0(&_S228, dt_0);



    thread GeodesicStateCart_0 s4_0;

#line 215
    thread Vec4Cart_0 _S230 = _S227;

#line 215
    Vec4Cart_0 _S231 = AddVec4Cart_0(&state_2->x_5, &_S230);
    (&s4_0)->x_5 = _S231;

#line 216
    thread Vec4Cart_0 _S232 = _S229;

#line 216
    Vec4Cart_0 _S233 = AddVec4Cart_0(&state_2->u_1, &_S232);
    (&s4_0)->u_1 = _S233;

#line 217
    thread GeodesicStateCart_0 _S234 = s4_0;

#line 217
    bool _S235 = TryGeodesicAccelerationWarpCart_0(&_S234, vs_4, sigma_6, R_6, xs_4, ys_4, zs_4, &acceleration_1);
    if(!_S235)
    {

#line 219
        return false;
    }

#line 219
    thread Vec4Cart_0 _S236 = (&s4_0)->u_1;

#line 219
    Vec4Cart_0 _S237 = ScaleVec4Cart_0(&_S236, dt_0);

#line 219
    thread Vec4Cart_0 _S238 = acceleration_1;

#line 219
    Vec4Cart_0 _S239 = ScaleVec4Cart_0(&_S238, dt_0);

#line 219
    thread Vec4Cart_0 _S240 = _S213;

#line 219
    Vec4Cart_0 _S241 = ScaleVec4Cart_0(&_S240, 2.0);

#line 219
    thread Vec4Cart_0 _S242 = _S199;

#line 219
    thread Vec4Cart_0 _S243 = _S241;

#line 219
    Vec4Cart_0 _S244 = AddVec4Cart_0(&_S242, &_S243);

#line 219
    thread Vec4Cart_0 _S245 = _S227;

#line 219
    Vec4Cart_0 _S246 = ScaleVec4Cart_0(&_S245, 2.0);

#line 219
    thread Vec4Cart_0 _S247 = _S246;

#line 219
    thread Vec4Cart_0 _S248 = _S237;

#line 219
    Vec4Cart_0 _S249 = AddVec4Cart_0(&_S247, &_S248);

#line 219
    thread Vec4Cart_0 _S250 = _S244;

#line 219
    thread Vec4Cart_0 _S251 = _S249;

#line 219
    Vec4Cart_0 _S252 = AddVec4Cart_0(&_S250, &_S251);

#line 219
    thread Vec4Cart_0 _S253 = _S252;

#line 219
    Vec4Cart_0 _S254 = ScaleVec4Cart_0(&_S253, 0.1666666716337204);

#line 224
    *dx_2 = _S254;

#line 224
    thread Vec4Cart_0 _S255 = _S215;

#line 224
    Vec4Cart_0 _S256 = ScaleVec4Cart_0(&_S255, 2.0);

#line 224
    thread Vec4Cart_0 _S257 = _S201;

#line 224
    thread Vec4Cart_0 _S258 = _S256;

#line 224
    Vec4Cart_0 _S259 = AddVec4Cart_0(&_S257, &_S258);

#line 224
    thread Vec4Cart_0 _S260 = _S229;

#line 224
    Vec4Cart_0 _S261 = ScaleVec4Cart_0(&_S260, 2.0);

#line 224
    thread Vec4Cart_0 _S262 = _S261;

#line 224
    thread Vec4Cart_0 _S263 = _S239;

#line 224
    Vec4Cart_0 _S264 = AddVec4Cart_0(&_S262, &_S263);

#line 224
    thread Vec4Cart_0 _S265 = _S259;

#line 224
    thread Vec4Cart_0 _S266 = _S264;

#line 224
    Vec4Cart_0 _S267 = AddVec4Cart_0(&_S265, &_S266);

#line 224
    thread Vec4Cart_0 _S268 = _S267;

#line 224
    Vec4Cart_0 _S269 = ScaleVec4Cart_0(&_S268, 0.1666666716337204);

    *du_0 = _S269;

#line 226
    thread Vec4Cart_0 _S270 = *dx_2;

#line 226
    Vec4Cart_0 _S271 = AddVec4Cart_0(&state_2->x_5, &_S270);

#line 226
    thread Vec4Cart_0 _S272 = _S269;

#line 226
    Vec4Cart_0 _S273 = AddVec4Cart_0(&state_2->u_1, &_S272);

#line 226
    thread Vec4Cart_0 _S274 = _S271;

#line 226
    bool _S275 = IsWarpDriveEventRepresented_0(&_S274, vs_4, sigma_6, R_6, xs_4, ys_4, zs_4);

#line 226
    bool _S276;



    if(_S275)
    {

#line 230
        _S276 = isfinite(_S273.t_0);

#line 230
    }
    else
    {

#line 230
        _S276 = false;

#line 230
    }
    if(_S276)
    {

#line 231
        _S276 = isfinite(_S273.x_2);

#line 231
    }
    else
    {

#line 231
        _S276 = false;

#line 231
    }

#line 231
    if(_S276)
    {

#line 231
        _S276 = isfinite(_S273.y_0);

#line 231
    }
    else
    {

#line 231
        _S276 = false;

#line 231
    }
    if(_S276)
    {

#line 232
        _S276 = isfinite(_S273.z_0);

#line 232
    }
    else
    {

#line 232
        _S276 = false;

#line 232
    }

#line 230
    return _S276;
}



bool IntegrateGeodesicRK4WarpCart_0(const GeodesicStateCart_0 thread* state_3, float dt_1, float vs_5, float sigma_7, float R_7, float xs_5, float ys_5, float zs_5, GeodesicStateCart_0 thread* result_4)
{

    *result_4 = *state_3;
    thread Vec4Cart_0 dx_3;

#line 239
    thread Vec4Cart_0 du_1;

#line 239
    bool _S277 = IntegrateGeodesicRK4WarpCartDelta_0(state_3, dt_1, vs_5, sigma_7, R_7, xs_5, ys_5, zs_5, &dx_3, &du_1);
    if(!_S277)
    {

#line 241
        return false;
    }

#line 241
    thread Vec4Cart_0 _S278 = dx_3;

#line 241
    Vec4Cart_0 _S279 = AddVec4Cart_0(&state_3->x_5, &_S278);
    result_4->x_5 = _S279;

#line 242
    thread Vec4Cart_0 _S280 = du_1;

#line 242
    Vec4Cart_0 _S281 = AddVec4Cart_0(&state_3->u_1, &_S280);
    result_4->u_1 = _S281;
    return true;
}


#line 403 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_metrics.slang"
void GetMorrisThorneCartesianChristoffel_0(const Vec4Cart_0 thread* pos_6, float b0_4, float Phi0_2, int shapeType_2, array<array<array<float, int(4)>, int(4)>, int(4)> thread* Gamma_2)
{
    int k_1;

#line 405
    int j_4;

#line 405
    int i_5 = int(0);

    for(;;)
    {

#line 407
        if(i_5 < int(4))
        {
        }
        else
        {

#line 407
            break;
        }

#line 407
        j_4 = int(0);
        for(;;)
        {

#line 408
            if(j_4 < int(4))
            {
            }
            else
            {

#line 408
                break;
            }

#line 408
            k_1 = int(0);
            for(;;)
            {

#line 409
                if(k_1 < int(4))
                {
                }
                else
                {

#line 409
                    break;
                }

#line 410
                (*Gamma_2)[i_5][j_4][k_1] = 0.0;

#line 409
                k_1 = k_1 + int(1);

#line 409
            }

#line 408
            j_4 = j_4 + int(1);

#line 408
        }

#line 407
        i_5 = i_5 + int(1);

#line 407
    }

#line 407
    bool _S282 = IsMorrisThorneCartesianEventRepresented_0(pos_6, b0_4, Phi0_2, shapeType_2);



    if(!_S282)
    {

#line 411
        return;
    }
    float rho2_1 = pos_6->x_2 * pos_6->x_2 + pos_6->y_0 * pos_6->y_0 + pos_6->z_0 * pos_6->z_0;
    float halfB0OverRho_1 = 0.5 * b0_4 / sqrt(rho2_1);
    float q_0 = halfB0OverRho_1 * halfB0OverRho_1;
    float _S283 = 1.0 + q_0;
    array<float, int(3)> _S284 = { pos_6->x_2, pos_6->y_0, pos_6->z_0 };
    thread array<float, int(3)> w_0;

#line 418
    k_1 = int(0);
    for(;;)
    {

#line 419
        if(k_1 < int(3))
        {
        }
        else
        {

#line 419
            break;
        }

#line 420
        w_0[k_1] = -2.0 * q_0 * _S284[k_1] / (_S283 * rho2_1);

#line 419
        k_1 = k_1 + int(1);

#line 419
    }

#line 419
    int l_1 = int(0);



    for(;;)
    {

#line 423
        if(l_1 < int(3))
        {
        }
        else
        {

#line 423
            break;
        }

#line 423
        i_5 = int(0);
        for(;;)
        {

#line 424
            if(i_5 < int(3))
            {
            }
            else
            {

#line 424
                break;
            }

#line 424
            j_4 = int(0);
            for(;;)
            {

#line 425
                if(j_4 < int(3))
                {
                }
                else
                {

#line 425
                    break;
                }

#line 425
                float _S285;
                if(l_1 == j_4)
                {

#line 426
                    _S285 = w_0[i_5];

#line 426
                }
                else
                {

#line 426
                    _S285 = 0.0;

#line 426
                }

#line 426
                float _S286;
                if(l_1 == i_5)
                {

#line 427
                    _S286 = w_0[j_4];

#line 427
                }
                else
                {

#line 427
                    _S286 = 0.0;

#line 427
                }

#line 427
                float _S287 = _S285 + _S286;

#line 427
                float _S288;
                if(i_5 == j_4)
                {

#line 428
                    _S288 = w_0[l_1];

#line 428
                }
                else
                {

#line 428
                    _S288 = 0.0;

#line 428
                }
                int _S289 = j_4 + int(1);

#line 429
                (*Gamma_2)[l_1 + int(1)][i_5 + int(1)][_S289] = _S287 - _S288;

#line 425
                j_4 = _S289;

#line 425
            }

#line 424
            i_5 = i_5 + int(1);

#line 424
        }

#line 423
        l_1 = l_1 + int(1);

#line 423
    }

#line 433
    return;
}


#line 93 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_integrator.slang"
bool TryGeodesicAccelerationMorrisThorneCart_0(const GeodesicStateCart_0 thread* state_4, float b0_5, float Phi0_3, int shapeType_3, Vec4Cart_0 thread* acceleration_2)
{

    *acceleration_2 = MakeVec4Cart_0(0.0, 0.0, 0.0, 0.0);

#line 96
    bool _S290 = IsMorrisThorneCartesianEventRepresented_0(&state_4->x_5, b0_5, Phi0_3, shapeType_3);
    if(!_S290)
    {

#line 97
        return false;
    }

#line 98
    thread array<array<array<float, int(4)>, int(4)>, int(4)> Gamma_3;

#line 98
    GetMorrisThorneCartesianChristoffel_0(&state_4->x_5, b0_5, Phi0_3, shapeType_3, &Gamma_3);

    array<float, int(4)> _S291 = { (&state_4->u_1)->t_0, (&state_4->u_1)->x_2, (&state_4->u_1)->y_0, (&state_4->u_1)->z_0 };
    thread array<float, int(4)> acc_1;

#line 101
    acc_1[int(0)] = 0.0;

#line 101
    acc_1[int(1)] = 0.0;

#line 101
    acc_1[int(2)] = 0.0;

#line 101
    acc_1[int(3)] = 0.0;

#line 101
    int mu_4 = int(0);
    for(;;)
    {

#line 102
        if(mu_4 < int(4))
        {
        }
        else
        {

#line 102
            break;
        }

#line 102
        int nu_5 = int(0);
        for(;;)
        {

#line 103
            if(nu_5 < int(4))
            {
            }
            else
            {

#line 103
                break;
            }

#line 103
            int rho_1 = int(0);
            for(;;)
            {

#line 104
                if(rho_1 < int(4))
                {
                }
                else
                {

#line 104
                    break;
                }

#line 105
                acc_1[mu_4] = acc_1[mu_4] - Gamma_3[mu_4][nu_5][rho_1] * _S291[nu_5] * _S291[rho_1];

#line 104
                rho_1 = rho_1 + int(1);

#line 104
            }

#line 103
            nu_5 = nu_5 + int(1);

#line 103
        }

#line 102
        mu_4 = mu_4 + int(1);

#line 102
    }



    *acceleration_2 = MakeVec4Cart_0(acc_1[int(0)], acc_1[int(1)], acc_1[int(2)], acc_1[int(3)]);
    return true;
}

bool IntegrateGeodesicRK4MorrisThorneCartDelta_0(const GeodesicStateCart_0 thread* state_5, float dt_2, float b0_6, float Phi0_4, int shapeType_4, Vec4Cart_0 thread* dx_4, Vec4Cart_0 thread* du_2)
{

    Vec4Cart_0 _S292 = MakeVec4Cart_0(0.0, 0.0, 0.0, 0.0);

#line 113
    *dx_4 = _S292;
    *du_2 = _S292;
    thread Vec4Cart_0 acceleration_3;

#line 115
    bool _S293 = TryGeodesicAccelerationMorrisThorneCart_0(state_5, b0_6, Phi0_4, shapeType_4, &acceleration_3);
    if(!_S293)
    {

#line 117
        return false;
    }

#line 117
    Vec4Cart_0 _S294 = ScaleVec4Cart_0(&state_5->u_1, dt_2);

#line 117
    thread Vec4Cart_0 _S295 = acceleration_3;

#line 117
    Vec4Cart_0 _S296 = ScaleVec4Cart_0(&_S295, dt_2);



    thread GeodesicStateCart_0 s2_1;

#line 121
    thread Vec4Cart_0 _S297 = _S294;

#line 121
    Vec4Cart_0 _S298 = ScaleVec4Cart_0(&_S297, 0.5);

#line 121
    thread Vec4Cart_0 _S299 = _S298;

#line 121
    Vec4Cart_0 _S300 = AddVec4Cart_0(&state_5->x_5, &_S299);
    (&s2_1)->x_5 = _S300;

#line 122
    thread Vec4Cart_0 _S301 = _S296;

#line 122
    Vec4Cart_0 _S302 = ScaleVec4Cart_0(&_S301, 0.5);

#line 122
    thread Vec4Cart_0 _S303 = _S302;

#line 122
    Vec4Cart_0 _S304 = AddVec4Cart_0(&state_5->u_1, &_S303);
    (&s2_1)->u_1 = _S304;

#line 123
    thread GeodesicStateCart_0 _S305 = s2_1;

#line 123
    bool _S306 = TryGeodesicAccelerationMorrisThorneCart_0(&_S305, b0_6, Phi0_4, shapeType_4, &acceleration_3);
    if(!_S306)
    {

#line 125
        return false;
    }

#line 125
    thread Vec4Cart_0 _S307 = (&s2_1)->u_1;

#line 125
    Vec4Cart_0 _S308 = ScaleVec4Cart_0(&_S307, dt_2);

#line 125
    thread Vec4Cart_0 _S309 = acceleration_3;

#line 125
    Vec4Cart_0 _S310 = ScaleVec4Cart_0(&_S309, dt_2);



    thread GeodesicStateCart_0 s3_1;

#line 129
    thread Vec4Cart_0 _S311 = _S308;

#line 129
    Vec4Cart_0 _S312 = ScaleVec4Cart_0(&_S311, 0.5);

#line 129
    thread Vec4Cart_0 _S313 = _S312;

#line 129
    Vec4Cart_0 _S314 = AddVec4Cart_0(&state_5->x_5, &_S313);
    (&s3_1)->x_5 = _S314;

#line 130
    thread Vec4Cart_0 _S315 = _S310;

#line 130
    Vec4Cart_0 _S316 = ScaleVec4Cart_0(&_S315, 0.5);

#line 130
    thread Vec4Cart_0 _S317 = _S316;

#line 130
    Vec4Cart_0 _S318 = AddVec4Cart_0(&state_5->u_1, &_S317);
    (&s3_1)->u_1 = _S318;

#line 131
    thread GeodesicStateCart_0 _S319 = s3_1;

#line 131
    bool _S320 = TryGeodesicAccelerationMorrisThorneCart_0(&_S319, b0_6, Phi0_4, shapeType_4, &acceleration_3);
    if(!_S320)
    {

#line 133
        return false;
    }

#line 133
    thread Vec4Cart_0 _S321 = (&s3_1)->u_1;

#line 133
    Vec4Cart_0 _S322 = ScaleVec4Cart_0(&_S321, dt_2);

#line 133
    thread Vec4Cart_0 _S323 = acceleration_3;

#line 133
    Vec4Cart_0 _S324 = ScaleVec4Cart_0(&_S323, dt_2);



    thread GeodesicStateCart_0 s4_1;

#line 137
    thread Vec4Cart_0 _S325 = _S322;

#line 137
    Vec4Cart_0 _S326 = AddVec4Cart_0(&state_5->x_5, &_S325);
    (&s4_1)->x_5 = _S326;

#line 138
    thread Vec4Cart_0 _S327 = _S324;

#line 138
    Vec4Cart_0 _S328 = AddVec4Cart_0(&state_5->u_1, &_S327);
    (&s4_1)->u_1 = _S328;

#line 139
    thread GeodesicStateCart_0 _S329 = s4_1;

#line 139
    bool _S330 = TryGeodesicAccelerationMorrisThorneCart_0(&_S329, b0_6, Phi0_4, shapeType_4, &acceleration_3);
    if(!_S330)
    {

#line 141
        return false;
    }

#line 141
    thread Vec4Cart_0 _S331 = (&s4_1)->u_1;

#line 141
    Vec4Cart_0 _S332 = ScaleVec4Cart_0(&_S331, dt_2);

#line 141
    thread Vec4Cart_0 _S333 = acceleration_3;

#line 141
    Vec4Cart_0 _S334 = ScaleVec4Cart_0(&_S333, dt_2);

#line 141
    thread Vec4Cart_0 _S335 = _S308;

#line 141
    Vec4Cart_0 _S336 = ScaleVec4Cart_0(&_S335, 2.0);

#line 141
    thread Vec4Cart_0 _S337 = _S294;

#line 141
    thread Vec4Cart_0 _S338 = _S336;

#line 141
    Vec4Cart_0 _S339 = AddVec4Cart_0(&_S337, &_S338);

#line 141
    thread Vec4Cart_0 _S340 = _S322;

#line 141
    Vec4Cart_0 _S341 = ScaleVec4Cart_0(&_S340, 2.0);

#line 141
    thread Vec4Cart_0 _S342 = _S341;

#line 141
    thread Vec4Cart_0 _S343 = _S332;

#line 141
    Vec4Cart_0 _S344 = AddVec4Cart_0(&_S342, &_S343);

#line 141
    thread Vec4Cart_0 _S345 = _S339;

#line 141
    thread Vec4Cart_0 _S346 = _S344;

#line 141
    Vec4Cart_0 _S347 = AddVec4Cart_0(&_S345, &_S346);

#line 141
    thread Vec4Cart_0 _S348 = _S347;

#line 141
    Vec4Cart_0 _S349 = ScaleVec4Cart_0(&_S348, 0.1666666716337204);

#line 146
    *dx_4 = _S349;

#line 146
    thread Vec4Cart_0 _S350 = _S310;

#line 146
    Vec4Cart_0 _S351 = ScaleVec4Cart_0(&_S350, 2.0);

#line 146
    thread Vec4Cart_0 _S352 = _S296;

#line 146
    thread Vec4Cart_0 _S353 = _S351;

#line 146
    Vec4Cart_0 _S354 = AddVec4Cart_0(&_S352, &_S353);

#line 146
    thread Vec4Cart_0 _S355 = _S324;

#line 146
    Vec4Cart_0 _S356 = ScaleVec4Cart_0(&_S355, 2.0);

#line 146
    thread Vec4Cart_0 _S357 = _S356;

#line 146
    thread Vec4Cart_0 _S358 = _S334;

#line 146
    Vec4Cart_0 _S359 = AddVec4Cart_0(&_S357, &_S358);

#line 146
    thread Vec4Cart_0 _S360 = _S354;

#line 146
    thread Vec4Cart_0 _S361 = _S359;

#line 146
    Vec4Cart_0 _S362 = AddVec4Cart_0(&_S360, &_S361);

#line 146
    thread Vec4Cart_0 _S363 = _S362;

#line 146
    Vec4Cart_0 _S364 = ScaleVec4Cart_0(&_S363, 0.1666666716337204);

    *du_2 = _S364;

#line 148
    thread Vec4Cart_0 _S365 = *dx_4;

#line 148
    Vec4Cart_0 _S366 = AddVec4Cart_0(&state_5->x_5, &_S365);

#line 148
    thread Vec4Cart_0 _S367 = _S366;

#line 148
    bool _S368 = IsMorrisThorneCartesianEventRepresented_0(&_S367, b0_6, Phi0_4, shapeType_4);


    return _S368;
}

bool IntegrateGeodesicRK4MorrisThorneCart_0(const GeodesicStateCart_0 thread* state_6, float dt_3, float b0_7, float Phi0_5, int shapeType_5, GeodesicStateCart_0 thread* result_5)
{

    *result_5 = *state_6;
    thread Vec4Cart_0 dx_5;

#line 158
    thread Vec4Cart_0 du_3;

#line 158
    bool _S369 = IntegrateGeodesicRK4MorrisThorneCartDelta_0(state_6, dt_3, b0_7, Phi0_5, shapeType_5, &dx_5, &du_3);
    if(!_S369)
    {

#line 160
        return false;
    }

#line 160
    thread Vec4Cart_0 _S370 = dx_5;

#line 160
    Vec4Cart_0 _S371 = AddVec4Cart_0(&state_6->x_5, &_S370);
    result_5->x_5 = _S371;

#line 161
    thread Vec4Cart_0 _S372 = du_3;

#line 161
    Vec4Cart_0 _S373 = AddVec4Cart_0(&state_6->u_1, &_S372);
    result_5->u_1 = _S373;
    return true;
}


#line 82 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_metrics.slang"
bool IsKerrSchildFamilyEventRepresented_0(const Vec4Cart_0 thread* event_0, float M_1, float a_5, float Q_1, float Lambda_1)
{

#line 83
    bool _S374;
    if(!isfinite(event_0->t_0))
    {

#line 84
        _S374 = true;

#line 84
    }
    else
    {

#line 84
        _S374 = !isfinite(event_0->x_2);

#line 84
    }

#line 84
    if(_S374)
    {

#line 84
        _S374 = true;

#line 84
    }
    else
    {

#line 84
        _S374 = !isfinite(event_0->y_0);

#line 84
    }

#line 84
    if(_S374)
    {

#line 84
        _S374 = true;

#line 84
    }
    else
    {

#line 84
        _S374 = !isfinite(event_0->z_0);

#line 84
    }
    if(_S374)
    {

#line 85
        _S374 = true;

#line 85
    }
    else
    {

#line 85
        _S374 = !isfinite(M_1);

#line 85
    }

#line 85
    if(_S374)
    {

#line 85
        _S374 = true;

#line 85
    }
    else
    {

#line 85
        _S374 = !isfinite(a_5);

#line 85
    }

#line 85
    if(_S374)
    {

#line 85
        _S374 = true;

#line 85
    }
    else
    {

#line 85
        _S374 = !isfinite(Q_1);

#line 85
    }
    if(_S374)
    {

#line 86
        _S374 = true;

#line 86
    }
    else
    {

#line 86
        _S374 = !isfinite(Lambda_1);

#line 86
    }

#line 84
    if(_S374)
    {

        return false;
    }
    bool _S375 = M_1 == 0.0;

#line 89
    if(_S375)
    {

#line 89
        _S374 = a_5 == 0.0;

#line 89
    }
    else
    {

#line 89
        _S374 = false;

#line 89
    }

#line 89
    if(_S374)
    {

#line 89
        _S374 = Q_1 == 0.0;

#line 89
    }
    else
    {

#line 89
        _S374 = false;

#line 89
    }

#line 89
    if(_S374)
    {

#line 89
        _S374 = Lambda_1 == 0.0;

#line 89
    }
    else
    {

#line 89
        _S374 = false;

#line 89
    }

#line 89
    if(_S374)
    {

#line 89
        return true;
    }

#line 90
    if(_S375)
    {

#line 90
        _S374 = a_5 == 0.0;

#line 90
    }
    else
    {

#line 90
        _S374 = false;

#line 90
    }

#line 90
    if(_S374)
    {

#line 90
        _S374 = Q_1 == 0.0;

#line 90
    }
    else
    {

#line 90
        _S374 = false;

#line 90
    }

#line 90
    if(_S374)
    {

#line 90
        _S374 = (event_0->x_2) == 0.0;

#line 90
    }
    else
    {

#line 90
        _S374 = false;

#line 90
    }
    if(_S374)
    {

#line 91
        _S374 = (event_0->y_0) == 0.0;

#line 91
    }
    else
    {

#line 91
        _S374 = false;

#line 91
    }

#line 91
    if(_S374)
    {

#line 91
        _S374 = (event_0->z_0) == 0.0;

#line 91
    }
    else
    {

#line 91
        _S374 = false;

#line 91
    }

#line 90
    if(_S374)
    {
        return true;
    }
    float r_1 = ComputeKerrRadiusCart_0(event_0->x_2, event_0->y_0, event_0->z_0, a_5);
    if(isfinite(r_1))
    {

#line 95
        _S374 = r_1 > 0.0;

#line 95
    }
    else
    {

#line 95
        _S374 = false;

#line 95
    }

#line 95
    return _S374;
}


#line 162
void GetKerrSchildFamilyChristoffel_0(const Vec4Cart_0 thread* pos_7, float M_2, float a_6, float Q_2, float Lambda_2, array<array<array<float, int(4)>, int(4)>, int(4)> thread* Gamma_4)
{
    int j_5;

#line 164
    int k_2;

#line 164
    float sum_2;

#line 164
    int nu_6;

#line 164
    int i_6 = int(0);

    for(;;)
    {

#line 166
        if(i_6 < int(4))
        {
        }
        else
        {

#line 166
            break;
        }

#line 166
        j_5 = int(0);
        for(;;)
        {

#line 167
            if(j_5 < int(4))
            {
            }
            else
            {

#line 167
                break;
            }

#line 167
            k_2 = int(0);
            for(;;)
            {

#line 168
                if(k_2 < int(4))
                {
                }
                else
                {

#line 168
                    break;
                }

#line 169
                (*Gamma_4)[i_6][j_5][k_2] = 0.0;

#line 168
                k_2 = k_2 + int(1);

#line 168
            }

#line 167
            j_5 = j_5 + int(1);

#line 167
        }

#line 166
        i_6 = i_6 + int(1);

#line 166
    }

#line 171
    bool _S376 = M_2 == 0.0;

#line 171
    bool _S377;

#line 171
    if(_S376)
    {

#line 171
        _S377 = a_6 == 0.0;

#line 171
    }
    else
    {

#line 171
        _S377 = false;

#line 171
    }

#line 171
    if(_S377)
    {

#line 171
        _S377 = Q_2 == 0.0;

#line 171
    }
    else
    {

#line 171
        _S377 = false;

#line 171
    }

#line 171
    if(_S377)
    {

#line 171
        _S377 = Lambda_2 == 0.0;

#line 171
    }
    else
    {

#line 171
        _S377 = false;

#line 171
    }

#line 171
    if(_S377)
    {

#line 172
        return;
    }

#line 172
    float _S378 = pos_7->x_2;

#line 172
    float _S379 = pos_7->y_0;

#line 172
    float _S380 = pos_7->z_0;

#line 179
    float a2_1 = a_6 * a_6;
    float Q2_1 = Q_2 * Q_2;
    float R2_0 = pos_7->x_2 * pos_7->x_2 + pos_7->y_0 * pos_7->y_0 + pos_7->z_0 * pos_7->z_0;

    float r_2 = ComputeKerrRadiusCart_0(pos_7->x_2, pos_7->y_0, pos_7->z_0, a_6);
    if(_S376)
    {

#line 184
        _S377 = a_6 == 0.0;

#line 184
    }
    else
    {

#line 184
        _S377 = false;

#line 184
    }

#line 184
    if(_S377)
    {

#line 184
        _S377 = Q_2 == 0.0;

#line 184
    }
    else
    {

#line 184
        _S377 = false;

#line 184
    }

#line 184
    if(_S377)
    {

#line 184
        _S377 = r_2 == 0.0;

#line 184
    }
    else
    {

#line 184
        _S377 = false;

#line 184
    }

#line 184
    if(_S377)
    {

#line 185
        return;
    }
    float r2_1 = r_2 * r_2;

    float Rm2_0 = R2_0 - a2_1;

#line 197
    float _S381 = 4.0 * sqrt(Rm2_0 * Rm2_0 + 4.0 * a2_1 * _S380 * _S380);

#line 197
    float d_r2_dx_0 = _S378 + 4.0 * _S378 * Rm2_0 / _S381;
    float d_r2_dy_0 = _S379 + 4.0 * _S379 * Rm2_0 / _S381;
    float d_r2_dz_0 = _S380 + 4.0 * _S380 * (R2_0 + a2_1) / _S381;

    thread array<float, int(4)> dr_0;
    dr_0[int(0)] = 0.0;
    float _S382 = 2.0 * r_2;

#line 203
    dr_0[int(1)] = d_r2_dx_0 / _S382;
    dr_0[int(2)] = d_r2_dy_0 / _S382;
    dr_0[int(3)] = d_r2_dz_0 / _S382;

    float denom_1 = r2_1 + a2_1;
    float denom2_0 = denom_1 * denom_1;


    float _S383 = r_2 * _S378 + a_6 * _S379;

#line 211
    float _S384 = _S383 / denom_1;
    float _S385 = r_2 * _S379 - a_6 * _S378;

#line 212
    float _S386 = _S385 / denom_1;
    float _S387 = _S380 / r_2;

#line 209
    array<float, int(4)> l_2 = { 1.0, _S384, _S386, _S387 };

#line 215
    array<float, int(4)> _S388 = { 1.0, - _S384, - _S386, - _S387 };

    thread array<array<float, int(4)>, int(4)> dl_0;

#line 217
    i_6 = int(0);
    for(;;)
    {

#line 218
        if(i_6 < int(4))
        {
        }
        else
        {

#line 218
            break;
        }

#line 218
        j_5 = int(0);
        for(;;)
        {

#line 219
            if(j_5 < int(4))
            {
            }
            else
            {

#line 219
                break;
            }

#line 219
            dl_0[i_6][j_5] = 0.0;

#line 219
            j_5 = j_5 + int(1);

#line 219
        }

#line 218
        i_6 = i_6 + int(1);

#line 218
    }

#line 224
    dl_0[int(1)][int(1)] = (r_2 + _S378 * dr_0[int(1)]) / denom_1 - _S383 * (_S382 * dr_0[int(1)]) / denom2_0;
    dl_0[int(2)][int(1)] = (_S378 * dr_0[int(2)] + a_6) / denom_1 - _S383 * (_S382 * dr_0[int(2)]) / denom2_0;
    dl_0[int(3)][int(1)] = _S378 * dr_0[int(3)] / denom_1 - _S383 * (_S382 * dr_0[int(3)]) / denom2_0;



    dl_0[int(1)][int(2)] = (_S379 * dr_0[int(1)] - a_6) / denom_1 - _S385 * (_S382 * dr_0[int(1)]) / denom2_0;
    dl_0[int(2)][int(2)] = (r_2 + _S379 * dr_0[int(2)]) / denom_1 - _S385 * (_S382 * dr_0[int(2)]) / denom2_0;
    dl_0[int(3)][int(2)] = _S379 * dr_0[int(3)] / denom_1 - _S385 * (_S382 * dr_0[int(3)]) / denom2_0;

    float _S389 = - _S380;

#line 234
    dl_0[int(1)][int(3)] = _S389 * dr_0[int(1)] / r2_1;
    dl_0[int(2)][int(3)] = _S389 * dr_0[int(2)] / r2_1;
    dl_0[int(3)][int(3)] = 1.0 / r_2 - _S380 * dr_0[int(3)] / r2_1;

    float cosine_1 = _S380 / r_2;
    float sigma_8 = r2_1 + a2_1 * cosine_1 * cosine_1;


    bool _S390 = M_2 != 0.0;

#line 242
    if(_S390)
    {

#line 242
        _S377 = true;

#line 242
    }
    else
    {

#line 242
        _S377 = Q_2 != 0.0;

#line 242
    }

#line 242
    float asymptoticH_0;

#line 242
    if(_S377)
    {

#line 242
        asymptoticH_0 = (2.0 * M_2 * r_2 - Q2_1) / sigma_8;

#line 242
    }
    else
    {

#line 242
        asymptoticH_0 = 0.0;

#line 242
    }



    bool _S391 = a_6 == 0.0;

#line 246
    if(_S391)
    {

#line 246
        _S377 = Lambda_2 != 0.0;

#line 246
    }
    else
    {

#line 246
        _S377 = false;

#line 246
    }

#line 246
    float H_1;

#line 246
    if(_S377)
    {

#line 246
        H_1 = asymptoticH_0 + Lambda_2 * r2_1 / 3.0;

#line 246
    }
    else
    {

#line 246
        H_1 = asymptoticH_0;

#line 246
    }



    thread array<float, int(4)> dH_0;

#line 250
    dH_0[int(0)] = 0.0;

#line 250
    dH_0[int(1)] = 0.0;

#line 250
    dH_0[int(2)] = 0.0;

#line 250
    dH_0[int(3)] = 0.0;

#line 250
    j_5 = int(1);
    for(;;)
    {

#line 251
        if(j_5 <= int(3))
        {
        }
        else
        {

#line 251
            break;
        }

#line 252
        if(_S390)
        {

#line 252
            _S377 = true;

#line 252
        }
        else
        {

#line 252
            _S377 = Q_2 != 0.0;

#line 252
        }

#line 252
        if(_S377)
        {

#line 253
            if(j_5 == int(3))
            {

#line 253
                sum_2 = 1.0 / r_2;

#line 253
            }
            else
            {

#line 253
                sum_2 = 0.0;

#line 253
            }



            dH_0[j_5] = (2.0 * M_2 * dr_0[j_5] - asymptoticH_0 * (_S382 * dr_0[j_5] + 2.0 * a2_1 * cosine_1 * (sum_2 - _S380 * dr_0[j_5] / r2_1))) / sigma_8;

#line 252
        }

#line 252
        bool _S392;

#line 259
        if(_S391)
        {

#line 259
            _S392 = Lambda_2 != 0.0;

#line 259
        }
        else
        {

#line 259
            _S392 = false;

#line 259
        }

#line 259
        if(_S392)
        {

#line 260
            dH_0[j_5] = dH_0[j_5] + 2.0 * Lambda_2 / 3.0 * r_2 * dr_0[j_5];

#line 259
        }

#line 251
        j_5 = j_5 + int(1);

#line 251
    }

#line 264
    thread array<array<array<float, int(4)>, int(4)>, int(4)> dg_1;

#line 264
    j_5 = int(0);
    for(;;)
    {

#line 265
        if(j_5 < int(4))
        {
        }
        else
        {

#line 265
            break;
        }

#line 265
        k_2 = int(0);
        for(;;)
        {

#line 266
            if(k_2 < int(4))
            {
            }
            else
            {

#line 266
                break;
            }

#line 266
            nu_6 = int(0);
            for(;;)
            {

#line 267
                if(nu_6 < int(4))
                {
                }
                else
                {

#line 267
                    break;
                }

#line 268
                dg_1[j_5][k_2][nu_6] = dH_0[j_5] * l_2[k_2] * l_2[nu_6] + H_1 * dl_0[j_5][k_2] * l_2[nu_6] + H_1 * l_2[k_2] * dl_0[j_5][nu_6];

#line 267
                nu_6 = nu_6 + int(1);

#line 267
            }

#line 266
            k_2 = k_2 + int(1);

#line 266
        }

#line 265
        j_5 = j_5 + int(1);

#line 265
    }

#line 272
    array<float, int(4)> _S393 = { -1.0, 1.0, 1.0, 1.0 };

#line 272
    j_5 = int(0);

    for(;;)
    {

#line 274
        if(j_5 < int(4))
        {
        }
        else
        {

#line 274
            break;
        }

#line 274
        k_2 = int(0);
        for(;;)
        {

#line 275
            if(k_2 < int(4))
            {
            }
            else
            {

#line 275
                break;
            }

#line 275
            nu_6 = k_2;
            for(;;)
            {

#line 276
                if(nu_6 < int(4))
                {
                }
                else
                {

#line 276
                    break;
                }

#line 276
                int s_2 = int(0);

#line 276
                sum_2 = 0.0;

                for(;;)
                {

#line 278
                    if(s_2 < int(4))
                    {
                    }
                    else
                    {

#line 278
                        break;
                    }

#line 278
                    float g_inv_lam_s_0;
                    if(j_5 == s_2)
                    {

#line 279
                        g_inv_lam_s_0 = _S393[j_5];

#line 279
                    }
                    else
                    {

#line 279
                        g_inv_lam_s_0 = 0.0;

#line 279
                    }


                    float sum_3 = sum_2 + (g_inv_lam_s_0 - H_1 * _S388[j_5] * _S388[s_2]) * (dg_1[k_2][s_2][nu_6] + dg_1[nu_6][s_2][k_2] - dg_1[s_2][k_2][nu_6]);

#line 278
                    s_2 = s_2 + int(1);

#line 278
                    sum_2 = sum_3;

#line 278
                }

#line 284
                float _S394 = 0.5 * sum_2;

#line 284
                (*Gamma_4)[j_5][k_2][nu_6] = _S394;
                (*Gamma_4)[j_5][nu_6][k_2] = _S394;

#line 276
                nu_6 = nu_6 + int(1);

#line 276
            }

#line 275
            k_2 = k_2 + int(1);

#line 275
        }

#line 274
        j_5 = j_5 + int(1);

#line 274
    }

#line 290
    return;
}


#line 13 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_integrator.slang"
bool TryGeodesicAccelerationCart_0(const GeodesicStateCart_0 thread* state_7, float M_3, float a_7, float Q_3, float Lambda_3, Vec4Cart_0 thread* acceleration_4)
{
    *acceleration_4 = MakeVec4Cart_0(0.0, 0.0, 0.0, 0.0);

#line 15
    bool _S395 = IsKerrSchildFamilyEventRepresented_0(&state_7->x_5, M_3, a_7, Q_3, Lambda_3);
    if(!_S395)
    {

#line 16
        return false;
    }

#line 17
    thread array<array<array<float, int(4)>, int(4)>, int(4)> Gamma_5;

#line 17
    GetKerrSchildFamilyChristoffel_0(&state_7->x_5, M_3, a_7, Q_3, Lambda_3, &Gamma_5);

    array<float, int(4)> _S396 = { (&state_7->u_1)->t_0, (&state_7->u_1)->x_2, (&state_7->u_1)->y_0, (&state_7->u_1)->z_0 };
    thread array<float, int(4)> acc_2;

#line 20
    acc_2[int(0)] = 0.0;

#line 20
    acc_2[int(1)] = 0.0;

#line 20
    acc_2[int(2)] = 0.0;

#line 20
    acc_2[int(3)] = 0.0;

#line 20
    int mu_5 = int(0);
    for(;;)
    {

#line 21
        if(mu_5 < int(4))
        {
        }
        else
        {

#line 21
            break;
        }

#line 21
        int nu_7 = int(0);
        for(;;)
        {

#line 22
            if(nu_7 < int(4))
            {
            }
            else
            {

#line 22
                break;
            }

#line 22
            int rho_2 = int(0);
            for(;;)
            {

#line 23
                if(rho_2 < int(4))
                {
                }
                else
                {

#line 23
                    break;
                }

#line 24
                acc_2[mu_5] = acc_2[mu_5] - Gamma_5[mu_5][nu_7][rho_2] * _S396[nu_7] * _S396[rho_2];

#line 23
                rho_2 = rho_2 + int(1);

#line 23
            }

#line 22
            nu_7 = nu_7 + int(1);

#line 22
        }

#line 21
        mu_5 = mu_5 + int(1);

#line 21
    }



    *acceleration_4 = MakeVec4Cart_0(acc_2[int(0)], acc_2[int(1)], acc_2[int(2)], acc_2[int(3)]);
    return true;
}


bool IntegrateGeodesicRK4CartDelta_0(const GeodesicStateCart_0 thread* state_8, float dt_4, float M_4, float a_8, float Q_4, float Lambda_4, Vec4Cart_0 thread* dx_6, Vec4Cart_0 thread* du_4)
{

    Vec4Cart_0 _S397 = MakeVec4Cart_0(0.0, 0.0, 0.0, 0.0);

#line 33
    *dx_6 = _S397;
    *du_4 = _S397;

#line 40
    thread Vec4Cart_0 acceleration_5;

#line 40
    Vec4Cart_0 previousDx_0 = *dx_6;

#line 40
    Vec4Cart_0 previousDu_0 = _S397;

#line 40
    Vec4Cart_0 dx1_0 = *dx_6;

#line 40
    Vec4Cart_0 dx2_0 = *dx_6;

#line 40
    Vec4Cart_0 dx3_0 = *dx_6;

#line 40
    Vec4Cart_0 dx4_0 = *dx_6;

#line 40
    Vec4Cart_0 du1_0 = _S397;

#line 40
    Vec4Cart_0 du2_0 = _S397;

#line 40
    Vec4Cart_0 du3_0 = _S397;

#line 40
    Vec4Cart_0 du4_0 = _S397;

#line 40
    int stage_0 = int(0);
    for(;;)
    {

#line 41
        if(stage_0 < int(4))
        {
        }
        else
        {

#line 41
            break;
        }

#line 42
        thread GeodesicStateCart_0 stageState_0 = *state_8;
        bool _S398 = stage_0 == int(1);

#line 43
        bool _S399;

#line 43
        if(_S398)
        {

#line 43
            _S399 = true;

#line 43
        }
        else
        {

#line 43
            _S399 = stage_0 == int(2);

#line 43
        }

#line 43
        if(_S399)
        {

#line 43
            thread Vec4Cart_0 _S400 = previousDx_0;

#line 43
            Vec4Cart_0 _S401 = ScaleVec4Cart_0(&_S400, 0.5);

#line 43
            thread Vec4Cart_0 _S402 = _S401;

#line 43
            Vec4Cart_0 _S403 = AddVec4Cart_0(&state_8->x_5, &_S402);
            (&stageState_0)->x_5 = _S403;

#line 44
            thread Vec4Cart_0 _S404 = previousDu_0;

#line 44
            Vec4Cart_0 _S405 = ScaleVec4Cart_0(&_S404, 0.5);

#line 44
            thread Vec4Cart_0 _S406 = _S405;

#line 44
            Vec4Cart_0 _S407 = AddVec4Cart_0(&state_8->u_1, &_S406);
            (&stageState_0)->u_1 = _S407;

#line 43
        }
        else
        {
            if(stage_0 == int(3))
            {

#line 46
                thread Vec4Cart_0 _S408 = previousDx_0;

#line 46
                Vec4Cart_0 _S409 = AddVec4Cart_0(&state_8->x_5, &_S408);
                (&stageState_0)->x_5 = _S409;

#line 47
                thread Vec4Cart_0 _S410 = previousDu_0;

#line 47
                Vec4Cart_0 _S411 = AddVec4Cart_0(&state_8->u_1, &_S410);
                (&stageState_0)->u_1 = _S411;

#line 46
            }

#line 43
        }

#line 43
        thread GeodesicStateCart_0 _S412 = stageState_0;

#line 43
        bool _S413 = TryGeodesicAccelerationCart_0(&_S412, M_4, a_8, Q_4, Lambda_4, &acceleration_5);

#line 50
        if(!_S413)
        {

#line 51
            return false;
        }

#line 51
        thread Vec4Cart_0 _S414 = (&stageState_0)->u_1;

#line 51
        Vec4Cart_0 _S415 = ScaleVec4Cart_0(&_S414, dt_4);

#line 51
        thread Vec4Cart_0 _S416 = acceleration_5;

#line 51
        Vec4Cart_0 _S417 = ScaleVec4Cart_0(&_S416, dt_4);


        if(stage_0 == int(0))
        {

#line 54
            dx1_0 = _S415;

#line 54
            du1_0 = _S417;

#line 54
        }
        else
        {

#line 54
            Vec4Cart_0 dx2_1;

#line 54
            Vec4Cart_0 dx3_1;

#line 54
            Vec4Cart_0 dx4_1;

#line 54
            Vec4Cart_0 du2_1;

#line 54
            Vec4Cart_0 du3_1;

#line 54
            Vec4Cart_0 du4_1;


            if(_S398)
            {

#line 57
                dx2_1 = _S415;

#line 57
                dx3_1 = dx3_0;

#line 57
                dx4_1 = dx4_0;

#line 57
                du2_1 = _S417;

#line 57
                du3_1 = du3_0;

#line 57
                du4_1 = du4_0;

#line 57
            }
            else
            {
                if(stage_0 == int(2))
                {

#line 60
                    dx2_1 = _S415;

#line 60
                    dx3_1 = dx4_0;

#line 60
                    dx4_1 = _S417;

#line 60
                    du2_1 = du4_0;

#line 60
                }
                else
                {

#line 60
                    dx2_1 = dx3_0;

#line 60
                    dx3_1 = _S415;

#line 60
                    dx4_1 = du3_0;

#line 60
                    du2_1 = _S417;

#line 60
                }

#line 37
                Vec4Cart_0 _S418 = dx2_1;

#line 37
                Vec4Cart_0 _S419 = dx3_1;
                Vec4Cart_0 _S420 = dx4_1;

#line 38
                Vec4Cart_0 _S421 = du2_1;

#line 38
                dx2_1 = dx2_0;

#line 38
                dx3_1 = _S418;

#line 38
                dx4_1 = _S419;

#line 38
                du2_1 = du2_0;

#line 38
                du3_1 = _S420;

#line 38
                du4_1 = _S421;

#line 57
            }

#line 57
            dx2_0 = dx2_1;

#line 57
            dx3_0 = dx3_1;

#line 57
            dx4_0 = dx4_1;

#line 57
            du2_0 = du2_1;

#line 57
            du3_0 = du3_1;

#line 57
            du4_0 = du4_1;

#line 54
        }

#line 41
        int stage_1 = stage_0 + int(1);

#line 41
        previousDx_0 = _S415;

#line 41
        previousDu_0 = _S417;

#line 41
        stage_0 = stage_1;

#line 41
    }

#line 41
    thread Vec4Cart_0 _S422 = dx2_0;

#line 41
    Vec4Cart_0 _S423 = ScaleVec4Cart_0(&_S422, 2.0);

#line 41
    thread Vec4Cart_0 _S424 = dx1_0;

#line 41
    thread Vec4Cart_0 _S425 = _S423;

#line 41
    Vec4Cart_0 _S426 = AddVec4Cart_0(&_S424, &_S425);

#line 41
    thread Vec4Cart_0 _S427 = dx3_0;

#line 41
    Vec4Cart_0 _S428 = ScaleVec4Cart_0(&_S427, 2.0);

#line 41
    thread Vec4Cart_0 _S429 = _S428;

#line 41
    thread Vec4Cart_0 _S430 = dx4_0;

#line 41
    Vec4Cart_0 _S431 = AddVec4Cart_0(&_S429, &_S430);

#line 41
    thread Vec4Cart_0 _S432 = _S426;

#line 41
    thread Vec4Cart_0 _S433 = _S431;

#line 41
    Vec4Cart_0 _S434 = AddVec4Cart_0(&_S432, &_S433);

#line 41
    thread Vec4Cart_0 _S435 = _S434;

#line 41
    Vec4Cart_0 _S436 = ScaleVec4Cart_0(&_S435, 0.1666666716337204);

#line 70
    *dx_6 = _S436;

#line 70
    thread Vec4Cart_0 _S437 = du2_0;

#line 70
    Vec4Cart_0 _S438 = ScaleVec4Cart_0(&_S437, 2.0);

#line 70
    thread Vec4Cart_0 _S439 = du1_0;

#line 70
    thread Vec4Cart_0 _S440 = _S438;

#line 70
    Vec4Cart_0 _S441 = AddVec4Cart_0(&_S439, &_S440);

#line 70
    thread Vec4Cart_0 _S442 = du3_0;

#line 70
    Vec4Cart_0 _S443 = ScaleVec4Cart_0(&_S442, 2.0);

#line 70
    thread Vec4Cart_0 _S444 = _S443;

#line 70
    thread Vec4Cart_0 _S445 = du4_0;

#line 70
    Vec4Cart_0 _S446 = AddVec4Cart_0(&_S444, &_S445);

#line 70
    thread Vec4Cart_0 _S447 = _S441;

#line 70
    thread Vec4Cart_0 _S448 = _S446;

#line 70
    Vec4Cart_0 _S449 = AddVec4Cart_0(&_S447, &_S448);

#line 70
    thread Vec4Cart_0 _S450 = _S449;

#line 70
    Vec4Cart_0 _S451 = ScaleVec4Cart_0(&_S450, 0.1666666716337204);

    *du_4 = _S451;

#line 72
    thread Vec4Cart_0 _S452 = *dx_6;

#line 72
    Vec4Cart_0 _S453 = AddVec4Cart_0(&state_8->x_5, &_S452);

#line 72
    thread Vec4Cart_0 _S454 = _S453;

#line 72
    bool _S455 = IsKerrSchildFamilyEventRepresented_0(&_S454, M_4, a_8, Q_4, Lambda_4);


    return _S455;
}

bool IntegrateGeodesicRK4Cart_0(const GeodesicStateCart_0 thread* state_9, float dt_5, float M_5, float a_9, float Q_5, float Lambda_5, GeodesicStateCart_0 thread* result_6)
{

    *result_6 = *state_9;
    thread Vec4Cart_0 dx_7;

#line 82
    thread Vec4Cart_0 du_5;

#line 82
    bool _S456 = IntegrateGeodesicRK4CartDelta_0(state_9, dt_5, M_5, a_9, Q_5, Lambda_5, &dx_7, &du_5);
    if(!_S456)
    {

#line 83
        return false;
    }

#line 83
    thread Vec4Cart_0 _S457 = dx_7;

#line 83
    Vec4Cart_0 _S458 = AddVec4Cart_0(&state_9->x_5, &_S457);
    result_6->x_5 = _S458;

#line 84
    thread Vec4Cart_0 _S459 = du_5;

#line 84
    Vec4Cart_0 _S460 = AddVec4Cart_0(&state_9->u_1, &_S459);
    result_6->u_1 = _S460;
    return true;
}


#line 259
bool AdvanceGeodesicRK4_0(const GeodesicStateCart_0 thread* state_10, float step_0, bool isWarp_0, bool isWormhole_0, float M_6, float a_10, float Q_6, float Lambda_6, float warpVs_0, float warpSigma_0, float warpR_0, float throatB0_0, int wormShape_1, GeodesicStateCart_0 thread* result_7)
{

#line 265
    if(isWarp_0)
    {

#line 265
        bool _S461 = IntegrateGeodesicRK4WarpCart_0(state_10, step_0, warpVs_0, warpSigma_0, warpR_0, 0.0, 0.0, 0.0, result_7);
        return _S461;
    }

    if(isWormhole_0)
    {

#line 269
        bool _S462 = IntegrateGeodesicRK4MorrisThorneCart_0(state_10, step_0, throatB0_0, 0.0, wormShape_1, result_7);
        return _S462;
    }

#line 270
    bool _S463 = IntegrateGeodesicRK4Cart_0(state_10, step_0, M_6, a_10, Q_6, Lambda_6, result_7);


    return _S463;
}


#line 447
void AdvanceCompositeGeodesicRK4Pair_0(const GeodesicStateCart_0 thread* state_11, float step_1, int coarseSubsteps_0, bool isWarp_1, bool isWormhole_1, float M_7, float a_11, float Q_7, float Lambda_7, float warpVs_1, float warpSigma_1, float warpR_1, float throatB0_1, int wormShape_2, GeodesicStateCart_0 thread* coarse_0, GeodesicStateCart_0 thread* fine_0, bool thread* coarseRepresented_0, bool thread* fineRepresented_0)
{



    GeodesicStateCart_0 current_0;
    *coarse_0 = *state_11;
    *fine_0 = *state_11;
    *coarseRepresented_0 = false;
    *fineRepresented_0 = false;

#line 456
    bool represented_0;
    if(coarseSubsteps_0 < int(1))
    {

#line 457
        represented_0 = true;

#line 457
    }
    else
    {

#line 457
        represented_0 = coarseSubsteps_0 > int(32);

#line 457
    }

#line 457
    if(represented_0)
    {

#line 457
        return;
    }

#line 457
    int pass_0 = int(0);
    for(;;)
    {

#line 458
        if(pass_0 < int(2))
        {
        }
        else
        {

#line 458
            break;
        }

#line 459
        bool _S464 = pass_0 == int(0);

#line 459
        int substeps_0;

#line 459
        if(_S464)
        {

#line 459
            substeps_0 = coarseSubsteps_0;

#line 459
        }
        else
        {

#line 459
            substeps_0 = coarseSubsteps_0 * int(2);

#line 459
        }
        float _S465 = step_1 / float(substeps_0);

#line 460
        GeodesicStateCart_0 current_1 = *state_11;

#line 460
        int part_0 = int(0);


        for(;;)
        {

#line 463
            if(part_0 < substeps_0)
            {
            }
            else
            {

#line 463
                current_0 = current_1;

#line 463
                represented_0 = true;

#line 463
                break;
            }

#line 463
            thread GeodesicStateCart_0 _S466 = current_1;
            thread GeodesicStateCart_0 next_0;

#line 464
            bool _S467 = AdvanceGeodesicRK4_0(&_S466, _S465, isWarp_1, isWormhole_1, M_7, a_11, Q_7, Lambda_7, warpVs_1, warpSigma_1, warpR_1, throatB0_1, wormShape_2, &next_0);
            if(!_S467)
            {

#line 465
                current_0 = *state_11;

#line 465
                represented_0 = false;



                break;
            }
            GeodesicStateCart_0 _S468 = next_0;

#line 463
            int part_1 = part_0 + int(1);

#line 463
            current_1 = _S468;

#line 463
            part_0 = part_1;

#line 463
        }

#line 473
        if(_S464)
        {

#line 474
            *coarse_0 = current_0;
            *coarseRepresented_0 = represented_0;

#line 473
        }
        else
        {

            *fine_0 = current_0;
            *fineRepresented_0 = represented_0;

#line 473
        }

#line 458
        pass_0 = pass_0 + int(1);

#line 458
    }

#line 481
    return;
}


#line 276
float RelativeStepDifference_0(const GeodesicStateCart_0 thread* start_0, const GeodesicStateCart_0 thread* first_2, const GeodesicStateCart_0 thread* second_2)
{
    array<float, int(8)> _S469 = { (&start_0->x_5)->t_0, (&start_0->x_5)->x_2, (&start_0->x_5)->y_0, (&start_0->x_5)->z_0, (&start_0->u_1)->t_0, (&start_0->u_1)->x_2, (&start_0->u_1)->y_0, (&start_0->u_1)->z_0 };

    array<float, int(8)> _S470 = { (&first_2->x_5)->t_0, (&first_2->x_5)->x_2, (&first_2->x_5)->y_0, (&first_2->x_5)->z_0, (&first_2->u_1)->t_0, (&first_2->u_1)->x_2, (&first_2->u_1)->y_0, (&first_2->u_1)->z_0 };

    array<float, int(8)> _S471 = { (&second_2->x_5)->t_0, (&second_2->x_5)->x_2, (&second_2->x_5)->y_0, (&second_2->x_5)->z_0, (&second_2->u_1)->t_0, (&second_2->u_1)->x_2, (&second_2->u_1)->y_0, (&second_2->u_1)->z_0 };

#line 282
    float error_0 = 0.0;

#line 282
    int component_0 = int(0);


    for(;;)
    {

#line 285
        if(component_0 < int(8))
        {
        }
        else
        {

#line 285
            break;
        }

        float _S472 = max(error_0, abs(_S471[component_0] - _S470[component_0]) / max(1.0, max(abs(_S469[component_0]), abs(_S471[component_0]))));

#line 285
        int component_1 = component_0 + int(1);

#line 285
        error_0 = _S472;

#line 285
        component_0 = component_1;

#line 285
    }

#line 292
    return error_0 / 15.0;
}

float RelativeNullConstraint_0(const array<array<float, int(4)>, int(4)> thread* metric_5, const Vec4Cart_0 thread* tangent_0)
{

#line 296
    array<float, int(4)> _S473 = { tangent_0->t_0, tangent_0->x_2, tangent_0->y_0, tangent_0->z_0 };

#line 296
    int mu_6 = int(0);

#line 296
    float contraction_0 = 0.0;

#line 296
    float scale_0 = 0.0;


    for(;;)
    {

#line 299
        if(mu_6 < int(4))
        {
        }
        else
        {

#line 299
            break;
        }

#line 299
        int nu_8 = int(0);
        for(;;)
        {

#line 300
            if(nu_8 < int(4))
            {
            }
            else
            {

#line 300
                break;
            }

#line 301
            float term_0 = (*metric_5)[mu_6][nu_8] * _S473[mu_6] * _S473[nu_8];
            float contraction_1 = contraction_0 + term_0;
            float scale_1 = scale_0 + abs(term_0);

#line 300
            nu_8 = nu_8 + int(1);

#line 300
            contraction_0 = contraction_1;

#line 300
            scale_0 = scale_1;

#line 300
        }

#line 299
        mu_6 = mu_6 + int(1);

#line 299
    }

#line 306
    return abs(contraction_0) / max(scale_0, 9.99999968265522539e-21);
}

float DeviceRealEpsilon_0()
{


    return 1.1920928955078125e-07;
}


#line 321
bool ProjectNullTangentWithSelection_0(const array<array<float, int(4)>, int(4)> thread* metric_6, const Vec4Cart_0 thread* tangent_1, Vec4Cart_0 thread* projected_0, int thread* selectedComponent_0)
{

#line 322
    int nu_9;

#line 322
    int j_6;
    *projected_0 = *tangent_1;
    *selectedComponent_0 = int(-1);

#line 324
    float _S474 = tangent_1->t_0;

#line 324
    bool temporalRepresented_0;
    if(!isfinite(tangent_1->t_0))
    {

#line 325
        temporalRepresented_0 = true;

#line 325
    }
    else
    {

#line 325
        temporalRepresented_0 = !isfinite(tangent_1->x_2);

#line 325
    }

#line 325
    if(temporalRepresented_0)
    {

#line 325
        temporalRepresented_0 = true;

#line 325
    }
    else
    {

#line 325
        temporalRepresented_0 = !isfinite(tangent_1->y_0);

#line 325
    }
    if(temporalRepresented_0)
    {

#line 326
        temporalRepresented_0 = true;

#line 326
    }
    else
    {

#line 326
        temporalRepresented_0 = !isfinite(tangent_1->z_0);

#line 326
    }

#line 325
    if(temporalRepresented_0)
    {
        return false;
    }

#line 327
    float _S475 = tangent_1->x_2;

#line 327
    float _S476 = tangent_1->y_0;

#line 327
    float _S477 = tangent_1->z_0;


    array<float, int(4)> components_0 = { _S474, tangent_1->x_2, tangent_1->y_0, tangent_1->z_0 };

#line 330
    float metricScale_0 = 0.0;

#line 330
    int mu_7 = int(0);

    for(;;)
    {

#line 332
        if(mu_7 < int(4))
        {
        }
        else
        {

#line 332
            break;
        }

#line 332
        nu_9 = int(0);
        for(;;)
        {

#line 333
            if(nu_9 < int(4))
            {
            }
            else
            {

#line 333
                break;
            }

#line 333
            float _S478 = (*metric_6)[mu_7][nu_9];
            if(!isfinite((*metric_6)[mu_7][nu_9]))
            {

#line 334
                return false;
            }

#line 335
            float _S479 = max(metricScale_0, abs(_S478));

#line 333
            int nu_10 = nu_9 + int(1);

#line 333
            metricScale_0 = _S479;

#line 333
            nu_9 = nu_10;

#line 333
        }

#line 332
        mu_7 = mu_7 + int(1);

#line 332
    }

#line 338
    if(!(metricScale_0 > 0.0))
    {

#line 338
        return false;
    }
    float _S480 = 64.0 * DeviceRealEpsilon_0();
    thread array<float, int(4)> bestComponents_0;

#line 341
    bestComponents_0[int(0)] = _S474;

#line 341
    bestComponents_0[int(1)] = _S475;

#line 341
    bestComponents_0[int(2)] = _S476;

#line 341
    bestComponents_0[int(3)] = _S477;

#line 341
    temporalRepresented_0 = false;

#line 341
    float bestCorrection_0 = 1.00000001504746622e+30;

#line 341
    bool bestRepresented_0 = false;

#line 341
    nu_9 = int(0);



    for(;;)
    {

#line 345
        if(nu_9 < int(4))
        {
        }
        else
        {

#line 345
            break;
        }

#line 345
        bool _S481;
        if(nu_9 > int(0))
        {

#line 346
            _S481 = temporalRepresented_0;

#line 346
        }
        else
        {

#line 346
            _S481 = false;

#line 346
        }

#line 346
        bool temporalRepresented_1;

#line 346
        bool bestRepresented_1;

#line 346
        float bestCorrection_1;

#line 346
        if(_S481)
        {

#line 346
            temporalRepresented_1 = temporalRepresented_0;

#line 346
            bestCorrection_1 = bestCorrection_0;

#line 346
            bestRepresented_1 = bestRepresented_0;

#line 345
            int component_2 = nu_9 + int(1);

#line 345
            temporalRepresented_0 = temporalRepresented_1;

#line 345
            bestCorrection_0 = bestCorrection_1;

#line 345
            bestRepresented_0 = bestRepresented_1;

#line 345
            nu_9 = component_2;

#line 345
            continue;
        }

#line 345
        float _S482 = (*metric_6)[nu_9][nu_9];

#line 345
        int i_7 = int(0);

#line 345
        float B_0 = 0.0;

#line 345
        float C_0 = 0.0;

#line 350
        for(;;)
        {

#line 350
            if(i_7 < int(4))
            {
            }
            else
            {

#line 350
                break;
            }

#line 351
            if(i_7 == nu_9)
            {

#line 351
                bestCorrection_1 = C_0;

#line 351
                i_7 = i_7 + int(1);

#line 351
                C_0 = bestCorrection_1;

#line 350
                continue;
            }

#line 350
            int _S483 = i_7;

            float B_1 = B_0 + ((*metric_6)[nu_9][i_7] + (*metric_6)[i_7][nu_9]) * components_0[i_7];

#line 352
            j_6 = int(0);

#line 352
            bestCorrection_1 = C_0;
            for(;;)
            {

#line 353
                if(j_6 < int(4))
                {
                }
                else
                {

#line 353
                    break;
                }

#line 354
                if(j_6 == nu_9)
                {

#line 354
                    j_6 = j_6 + int(1);

#line 353
                    continue;
                }

#line 353
                bestCorrection_1 = bestCorrection_1 + (*metric_6)[i_7][j_6] * components_0[_S483] * components_0[j_6];

#line 353
                j_6 = j_6 + int(1);

#line 353
            }

#line 353
            B_0 = B_1;

#line 350
            i_7 = i_7 + int(1);

#line 350
            C_0 = bestCorrection_1;

#line 350
        }

#line 350
        bool _S484;

#line 358
        if(!isfinite(_S482))
        {

#line 358
            _S484 = true;

#line 358
        }
        else
        {

#line 358
            _S484 = !isfinite(B_0);

#line 358
        }

#line 358
        bool _S485;

#line 358
        if(_S484)
        {

#line 358
            _S485 = true;

#line 358
        }
        else
        {

#line 358
            _S485 = !isfinite(C_0);

#line 358
        }

#line 358
        if(_S485)
        {

#line 358
            temporalRepresented_1 = temporalRepresented_0;

#line 358
            bestCorrection_1 = bestCorrection_0;

#line 358
            bestRepresented_1 = bestRepresented_0;

#line 345
            int component_2 = nu_9 + int(1);

#line 345
            temporalRepresented_0 = temporalRepresented_1;

#line 345
            bestCorrection_0 = bestCorrection_1;

#line 345
            bestRepresented_0 = bestRepresented_1;

#line 345
            nu_9 = component_2;

#line 345
            continue;
        }

#line 360
        thread array<float, int(2)> roots_0;

#line 360
        roots_0[int(0)] = 0.0;

#line 360
        roots_0[int(1)] = 0.0;

        if(_S482 == 0.0)
        {

#line 363
            if(B_0 == 0.0)
            {

#line 363
                temporalRepresented_1 = temporalRepresented_0;

#line 363
                bestCorrection_1 = bestCorrection_0;

#line 363
                bestRepresented_1 = bestRepresented_0;

#line 345
                int component_2 = nu_9 + int(1);

#line 345
                temporalRepresented_0 = temporalRepresented_1;

#line 345
                bestCorrection_0 = bestCorrection_1;

#line 345
                bestRepresented_0 = bestRepresented_1;

#line 345
                nu_9 = component_2;

#line 345
                continue;
            }

#line 364
            roots_0[int(0)] = - C_0 / B_0;

#line 364
            j_6 = int(1);

#line 362
        }
        else
        {

            float fourAC_0 = 4.0 * _S482 * C_0;
            float _S486 = B_0 * B_0;

#line 367
            float discriminant_0 = _S486 - fourAC_0;
            float discriminantScale_0 = _S486 + abs(fourAC_0);
            if(!isfinite(discriminant_0))
            {

#line 369
                temporalRepresented_1 = true;

#line 369
            }
            else
            {

#line 369
                temporalRepresented_1 = discriminant_0 < (- _S480 * discriminantScale_0);

#line 369
            }

#line 369
            if(temporalRepresented_1)
            {

#line 369
                temporalRepresented_1 = temporalRepresented_0;

#line 369
                bestCorrection_1 = bestCorrection_0;

#line 369
                bestRepresented_1 = bestRepresented_0;

#line 345
                int component_2 = nu_9 + int(1);

#line 345
                temporalRepresented_0 = temporalRepresented_1;

#line 345
                bestCorrection_0 = bestCorrection_1;

#line 345
                bestRepresented_0 = bestRepresented_1;

#line 345
                nu_9 = component_2;

#line 345
                continue;
            }

#line 374
            float sqrtDiscriminant_0 = sqrt(max(discriminant_0, 0.0));

            if(B_0 >= 0.0)
            {

#line 376
                metricScale_0 = sqrtDiscriminant_0;

#line 376
            }
            else
            {

#line 376
                metricScale_0 = - sqrtDiscriminant_0;

#line 376
            }

#line 375
            float q_1 = -0.5 * (B_0 + metricScale_0);

            if(q_1 == 0.0)
            {

#line 378
                float _S487 = - B_0 / (2.0 * _S482);

#line 378
                roots_0[int(0)] = _S487;
                roots_0[int(1)] = _S487;

#line 377
            }
            else
            {

                roots_0[int(0)] = q_1 / _S482;
                roots_0[int(1)] = C_0 / q_1;

#line 377
            }

#line 377
            j_6 = int(2);

#line 362
        }

#line 362
        bestCorrection_1 = bestCorrection_0;

#line 362
        temporalRepresented_1 = temporalRepresented_0;

#line 362
        bestRepresented_1 = bestRepresented_0;

#line 362
        int rootIndex_0 = int(0);

#line 386
        for(;;)
        {

#line 386
            if(rootIndex_0 < int(2))
            {
            }
            else
            {

#line 386
                break;
            }

#line 387
            if(rootIndex_0 >= j_6)
            {

#line 387
                rootIndex_0 = rootIndex_0 + int(1);

#line 386
                continue;
            }
            float root_0 = roots_0[rootIndex_0];
            if(!isfinite(roots_0[rootIndex_0]))
            {

#line 389
                rootIndex_0 = rootIndex_0 + int(1);

#line 386
                continue;
            }



            float correction_0 = abs(root_0 - components_0[nu_9]) / (1.0 + abs(components_0[nu_9]));

            if(correction_0 >= bestCorrection_1)
            {

#line 393
                rootIndex_0 = rootIndex_0 + int(1);

#line 386
                continue;
            }

#line 386
            int i_8 = int(0);

#line 394
            for(;;)
            {

#line 394
                if(i_8 < int(4))
                {
                }
                else
                {

#line 394
                    break;
                }

#line 394
                bestComponents_0[i_8] = components_0[i_8];

#line 394
                i_8 = i_8 + int(1);

#line 394
            }
            bestComponents_0[nu_9] = root_0;

            *selectedComponent_0 = nu_9;

#line 397
            bool temporalRepresented_2;

            if(nu_9 == int(0))
            {

#line 399
                temporalRepresented_2 = true;

#line 399
            }
            else
            {

#line 399
                temporalRepresented_2 = temporalRepresented_1;

#line 399
            }

#line 399
            bestCorrection_1 = correction_0;

#line 399
            temporalRepresented_1 = temporalRepresented_2;

#line 399
            bestRepresented_1 = true;

#line 386
            rootIndex_0 = rootIndex_0 + int(1);

#line 386
        }

#line 345
        int component_2 = nu_9 + int(1);

#line 345
        temporalRepresented_0 = temporalRepresented_1;

#line 345
        bestCorrection_0 = bestCorrection_1;

#line 345
        bestRepresented_0 = bestRepresented_1;

#line 345
        nu_9 = component_2;

#line 345
    }

#line 402
    if(!bestRepresented_0)
    {

#line 402
        return false;
    }

#line 403
    *projected_0 = MakeVec4Cart_0(bestComponents_0[int(0)], bestComponents_0[int(1)], bestComponents_0[int(2)], bestComponents_0[int(3)]);

    return true;
}

bool ProjectNullTangentPreservingBranch_0(const array<array<float, int(4)>, int(4)> thread* metric_7, const Vec4Cart_0 thread* tangent_2, Vec4Cart_0 thread* projected_1)
{
    thread int selectedComponent_1;

#line 410
    bool _S488 = ProjectNullTangentWithSelection_0(metric_7, tangent_2, projected_1, &selectedComponent_1);
    return _S488;
}


#line 247
struct AdaptiveGeodesicStep_0
{
    GeodesicStateCart_0 candidate_0;
    float errorRatio_0;
    float nullResidual_0;
    float projectedNullResidual_0;
    float stepFactor_0;
    bool accepted_1;
};


#line 499
AdaptiveGeodesicStep_0 AttemptAdaptiveGeodesicStepWithResolution_0(const GeodesicStateCart_0 thread* state_12, float step_2, bool isWarp_2, bool isWormhole_2, float M_8, float a_12, float Q_8, float Lambda_8, float warpVs_2, float warpSigma_2, float warpR_2, float throatB0_2, int wormShape_3, float integrationTolerance_0, float nullTolerance_0, int coarseSubsteps_1)
{

#line 505
    thread AdaptiveGeodesicStep_0 result_8;
    (&result_8)->candidate_0 = *state_12;
    thread GeodesicStateCart_0 fullStep_0 = *state_12;
    thread bool fullRepresented_0;



    thread bool candidateRepresented_0;

#line 512
    AdvanceCompositeGeodesicRK4Pair_0(state_12, step_2, coarseSubsteps_1, isWarp_2, isWormhole_2, M_8, a_12, Q_8, Lambda_8, warpVs_2, warpSigma_2, warpR_2, throatB0_2, wormShape_3, &fullStep_0, &(&result_8)->candidate_0, &fullRepresented_0, &candidateRepresented_0);

#line 512
    bool _S489;

#line 517
    if(!fullRepresented_0)
    {

#line 517
        _S489 = true;

#line 517
    }
    else
    {

#line 517
        _S489 = !candidateRepresented_0;

#line 517
    }

#line 517
    if(_S489)
    {

#line 518
        (&result_8)->candidate_0 = *state_12;
        (&result_8)->errorRatio_0 = 1.00000001504746622e+30;
        (&result_8)->nullResidual_0 = 1.00000001504746622e+30;
        (&result_8)->projectedNullResidual_0 = 1.00000001504746622e+30;
        (&result_8)->stepFactor_0 = 0.5;
        (&result_8)->accepted_1 = false;
        return result_8;
    }

#line 524
    thread GeodesicStateCart_0 _S490 = fullStep_0;

#line 524
    thread GeodesicStateCart_0 _S491 = (&result_8)->candidate_0;

#line 524
    float _S492 = RelativeStepDifference_0(state_12, &_S490, &_S491);

    (&result_8)->errorRatio_0 = _S492 / integrationTolerance_0;


    thread array<array<float, int(4)>, int(4)> candidateMetric_0;
    thread array<array<float, int(4)>, int(4)> candidateInverse_0;
    if(isWarp_2)
    {

#line 531
        thread Vec4Cart_0 _S493 = (&(&result_8)->candidate_0)->x_5;

#line 531
        GetWarpDriveMetric_0(&_S493, warpVs_2, warpSigma_2, warpR_2, 0.0, 0.0, 0.0, &candidateMetric_0, &candidateInverse_0);

#line 531
    }
    else
    {
        if(isWormhole_2)
        {

#line 534
            thread Vec4Cart_0 _S494 = (&(&result_8)->candidate_0)->x_5;

#line 534
            GetMorrisThorneCartesianMetric_0(&_S494, throatB0_2, 0.0, wormShape_3, &candidateMetric_0, &candidateInverse_0);

#line 534
        }
        else
        {

#line 534
            thread Vec4Cart_0 _S495 = (&(&result_8)->candidate_0)->x_5;

#line 534
            GetKerrSchildFamilyMetric_0(&_S495, M_8, a_12, Q_8, Lambda_8, &candidateMetric_0, &candidateInverse_0);

#line 534
        }

#line 531
    }

#line 531
    thread array<array<float, int(4)>, int(4)> _S496 = candidateMetric_0;

#line 531
    thread Vec4Cart_0 _S497 = (&(&result_8)->candidate_0)->u_1;

#line 531
    float _S498 = RelativeNullConstraint_0(&_S496, &_S497);

#line 541
    (&result_8)->nullResidual_0 = _S498;

#line 541
    thread array<array<float, int(4)>, int(4)> _S499 = candidateMetric_0;

#line 541
    thread Vec4Cart_0 _S500 = (&(&result_8)->candidate_0)->u_1;
    thread Vec4Cart_0 projectedTangent_0;

#line 542
    bool _S501 = ProjectNullTangentPreservingBranch_0(&_S499, &_S500, &projectedTangent_0);


    if(_S501)
    {

#line 546
        (&(&result_8)->candidate_0)->u_1 = projectedTangent_0;

#line 546
        thread array<array<float, int(4)>, int(4)> _S502 = candidateMetric_0;

#line 546
        thread Vec4Cart_0 _S503 = projectedTangent_0;

#line 546
        float _S504 = RelativeNullConstraint_0(&_S502, &_S503);
        (&result_8)->projectedNullResidual_0 = _S504;

#line 545
    }
    else
    {


        (&result_8)->projectedNullResidual_0 = 1.00000001504746622e+30;

#line 545
    }

#line 553
    (&result_8)->stepFactor_0 = clamp(0.89999997615814209 * pow(1.0 / max((&result_8)->errorRatio_0, 9.999999960041972e-13), 0.20000000298023224), 0.20000000298023224, 2.0);

    if(_S501)
    {

#line 555
        _S489 = ((&result_8)->errorRatio_0) <= 1.0;

#line 555
    }
    else
    {

#line 555
        _S489 = false;

#line 555
    }

#line 555
    if(_S489)
    {

#line 555
        _S489 = ((&result_8)->nullResidual_0) <= nullTolerance_0;

#line 555
    }
    else
    {

#line 555
        _S489 = false;

#line 555
    }
    if(_S489)
    {

#line 556
        _S489 = ((&result_8)->projectedNullResidual_0) <= nullTolerance_0;

#line 556
    }
    else
    {

#line 556
        _S489 = false;

#line 556
    }

#line 555
    (&result_8)->accepted_1 = _S489;


    if(!_S489)
    {

#line 559
        if(!_S501)
        {

#line 559
            _S489 = true;

#line 559
        }
        else
        {

#line 559
            _S489 = ((&result_8)->nullResidual_0) > nullTolerance_0;

#line 559
        }

#line 559
        if(_S489)
        {

#line 559
            _S489 = true;

#line 559
        }
        else
        {

#line 559
            _S489 = ((&result_8)->projectedNullResidual_0) > nullTolerance_0;

#line 559
        }

#line 558
    }
    else
    {

#line 558
        _S489 = false;

#line 558
    }

#line 558
    if(_S489)
    {

        (&result_8)->stepFactor_0 = min((&result_8)->stepFactor_0, 0.5);

#line 558
    }

#line 563
    return result_8;
}


AdaptiveGeodesicStep_0 AttemptAdaptiveGeodesicStep_0(const GeodesicStateCart_0 thread* state_13, float step_3, bool isWarp_3, bool isWormhole_3, float M_9, float a_13, float Q_9, float Lambda_9, float warpVs_3, float warpSigma_3, float warpR_3, float throatB0_3, int wormShape_4, float integrationTolerance_1, float nullTolerance_1)
{


    AdaptiveGeodesicStep_0 _S505 = AttemptAdaptiveGeodesicStepWithResolution_0(state_13, step_3, isWarp_3, isWormhole_3, M_9, a_13, Q_9, Lambda_9, warpVs_3, warpSigma_3, warpR_3, throatB0_3, wormShape_4, integrationTolerance_1, nullTolerance_1, int(1));
    return _S505;
}


#line 439
bool AdaptiveAttemptAlreadyAtMinimumStep_0(float attemptedStep_0, float minimumStep_0)
{

#line 440
    return attemptedStep_0 <= minimumStep_0;
}


#line 8 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_trace_event.slang"
bool IsFiniteAcceptedSegmentValue_0(const Vec4Cart_0 thread* value_0)
{

#line 8
    bool _S506;
    if(isfinite(value_0->t_0))
    {

#line 9
        _S506 = isfinite(value_0->x_2);

#line 9
    }
    else
    {

#line 9
        _S506 = false;

#line 9
    }

#line 9
    if(_S506)
    {

#line 9
        _S506 = isfinite(value_0->y_0);

#line 9
    }
    else
    {

#line 9
        _S506 = false;

#line 9
    }
    if(_S506)
    {

#line 10
        _S506 = isfinite(value_0->z_0);

#line 10
    }
    else
    {

#line 10
        _S506 = false;

#line 10
    }

#line 9
    return _S506;
}


#line 44
float PolynomialEpsilon_0()
{


    return 1.1920928955078125e-07;
}


#line 207
void SphericalSegmentPolynomial_0(const Vec4Cart_0 thread* start_1, const Vec4Cart_0 thread* startVelocity_0, const Vec4Cart_0 thread* finish_0, const Vec4Cart_0 thread* finishVelocity_0, float affineLength_1, float captureRadius_0, array<float, int(7)> thread* coefficients_0)
{

    int index_0 = int(0);
    for(;;)
    {

#line 211
        if(index_0 < int(7))
        {
        }
        else
        {

#line 211
            break;
        }

#line 211
        (*coefficients_0)[index_0] = 0.0;

#line 211
        index_0 = index_0 + int(1);

#line 211
    }
    (*coefficients_0)[int(0)] = - captureRadius_0 * captureRadius_0;
    array<float, int(3)> _S507 = { start_1->x_2, start_1->y_0, start_1->z_0 };
    array<float, int(3)> _S508 = { finish_0->x_2, finish_0->y_0, finish_0->z_0 };
    array<float, int(3)> _S509 = { startVelocity_0->x_2, startVelocity_0->y_0, startVelocity_0->z_0 };
    array<float, int(3)> _S510 = { finishVelocity_0->x_2, finishVelocity_0->y_0, finishVelocity_0->z_0 };

#line 216
    int axis_0 = int(0);
    for(;;)
    {

#line 217
        if(axis_0 < int(3))
        {
        }
        else
        {

#line 217
            break;
        }

#line 218
        float startDelta_0 = affineLength_1 * _S509[axis_0];
        float finishDelta_0 = affineLength_1 * _S510[axis_0];
        thread array<float, int(4)> component_3;
        component_3[int(0)] = _S507[axis_0];
        component_3[int(1)] = startDelta_0;
        component_3[int(2)] = -3.0 * _S507[axis_0] + 3.0 * _S508[axis_0] - 2.0 * startDelta_0 - finishDelta_0;


        component_3[int(3)] = 2.0 * _S507[axis_0] - 2.0 * _S508[axis_0] + startDelta_0 + finishDelta_0;

#line 226
        int left_0 = int(0);


        for(;;)
        {

#line 229
            if(left_0 <= int(3))
            {
            }
            else
            {

#line 229
                break;
            }

#line 229
            int right_1 = int(0);
            for(;;)
            {

#line 230
                if(right_1 <= int(3))
                {
                }
                else
                {

#line 230
                    break;
                }

#line 231
                (*coefficients_0)[left_0 + right_1] = (*coefficients_0)[left_0 + right_1] + component_3[left_0] * component_3[right_1];

#line 230
                right_1 = right_1 + int(1);

#line 230
            }

#line 229
            left_0 = left_0 + int(1);

#line 229
        }

#line 217
        axis_0 = axis_0 + int(1);

#line 217
    }

#line 233
    return;
}


#line 60
float PolynomialValueTolerance_0(const array<float, int(7)> thread* coefficients_1, int degree_0)
{

#line 60
    float scale_2 = 0.0;

#line 60
    int power_0 = int(0);

    for(;;)
    {

#line 62
        if(power_0 <= degree_0)
        {
        }
        else
        {

#line 62
            break;
        }

#line 63
        float _S511 = max(scale_2, abs((*coefficients_1)[power_0]));

#line 62
        int power_1 = power_0 + int(1);

#line 62
        scale_2 = _S511;

#line 62
        power_0 = power_1;

#line 62
    }


    return 256.0 * PolynomialEpsilon_0() * max(scale_2, 1.00000000317107685e-30);
}

void AppendUnitRoot_0(array<float, int(6)> thread* roots_1, int thread* rootCount_0, float root_1, float fractionTolerance_0)
{
    float _S512 = clamp(root_1, 0.0, 1.0);

#line 70
    bool _S513;
    if((*rootCount_0) > int(0))
    {

#line 71
        _S513 = (abs(_S512 - (*roots_1)[*rootCount_0 - int(1)])) <= fractionTolerance_0;

#line 71
    }
    else
    {

#line 71
        _S513 = false;

#line 71
    }

#line 71
    if(_S513)
    {

#line 71
        return;
    }

#line 72
    if((*rootCount_0) < int(6))
    {

#line 72
        int _S514 = *rootCount_0;

#line 72
        *rootCount_0 = *rootCount_0 + int(1);

#line 72
        (*roots_1)[_S514] = _S512;

#line 72
    }
    return;
}


#line 75
int FindLinearRoots_0(const array<float, int(7)> thread* coefficients_2, array<float, int(6)> thread* roots_2)
{

#line 75
    int index_1 = int(0);
    for(;;)
    {

#line 76
        if(index_1 < int(6))
        {
        }
        else
        {

#line 76
            break;
        }

#line 76
        (*roots_2)[index_1] = 0.0;

#line 76
        index_1 = index_1 + int(1);

#line 76
    }

#line 76
    float _S515 = PolynomialValueTolerance_0(coefficients_2, int(1));

    float fractionTolerance_1 = 512.0 * PolynomialEpsilon_0();

#line 78
    float _S516 = (*coefficients_2)[int(1)];
    if((abs((*coefficients_2)[int(1)])) <= _S515)
    {

#line 79
        return int(0);
    }

#line 80
    float root_2 = - (*coefficients_2)[int(0)] / _S516;
    thread int count_0 = int(0);

#line 81
    bool _S517;
    if(isfinite(root_2))
    {

#line 82
        _S517 = root_2 >= (- fractionTolerance_1);

#line 82
    }
    else
    {

#line 82
        _S517 = false;

#line 82
    }

#line 82
    if(_S517)
    {

#line 82
        _S517 = root_2 <= (1.0 + fractionTolerance_1);

#line 82
    }
    else
    {

#line 82
        _S517 = false;

#line 82
    }

#line 82
    if(_S517)
    {
        AppendUnitRoot_0(roots_2, &count_0, root_2, fractionTolerance_1);

#line 82
    }



    return count_0;
}


#line 52
float EvaluatePolynomial_0(const array<float, int(7)> thread* coefficients_3, int degree_1, float argument_0)
{
    int _S518 = degree_1 - int(1);

#line 54
    float value_1 = (*coefficients_3)[degree_1];

#line 54
    int power_2 = _S518;

#line 54
    for(;;)
    {

#line 54
        if(power_2 >= int(0))
        {
        }
        else
        {

#line 54
            break;
        }

#line 55
        float _S519 = value_1 * argument_0 + (*coefficients_3)[power_2];

#line 54
        int power_3 = power_2 - int(1);

#line 54
        value_1 = _S519;

#line 54
        power_2 = power_3;

#line 54
    }


    return value_1;
}


#line 89
int FindRootsFromCriticalPoints_0(const array<float, int(7)> thread* coefficients_4, int degree_2, const array<float, int(6)> thread* criticalRoots_0, int criticalCount_0, array<float, int(6)> thread* roots_3)
{
    bool _S520;

#line 91
    int index_2 = int(0);
    for(;;)
    {

#line 92
        if(index_2 < int(6))
        {
        }
        else
        {

#line 92
            break;
        }

#line 92
        (*roots_3)[index_2] = 0.0;

#line 92
        index_2 = index_2 + int(1);

#line 92
    }

#line 92
    float _S521 = PolynomialValueTolerance_0(coefficients_4, degree_2);

    float fractionTolerance_2 = 512.0 * PolynomialEpsilon_0();
    thread array<float, int(8)> boundaries_0;

#line 95
    index_2 = int(0);
    for(;;)
    {

#line 96
        if(index_2 < int(8))
        {
        }
        else
        {

#line 96
            break;
        }

#line 96
        boundaries_0[index_2] = 0.0;

#line 96
        index_2 = index_2 + int(1);

#line 96
    }

    boundaries_0[int(0)] = 0.0;

#line 98
    index_2 = int(0);

#line 98
    int boundaryCount_0 = int(1);
    for(;;)
    {

#line 99
        if(index_2 < criticalCount_0)
        {
        }
        else
        {

#line 99
            break;
        }

#line 99
        float _S522 = (*criticalRoots_0)[index_2];

        if(((*criticalRoots_0)[index_2]) > fractionTolerance_2)
        {

#line 101
            _S520 = _S522 < (1.0 - fractionTolerance_2);

#line 101
        }
        else
        {

#line 101
            _S520 = false;

#line 101
        }

#line 101
        if(_S520)
        {

#line 102
            int _S523 = boundaryCount_0 + int(1);

#line 102
            boundaries_0[boundaryCount_0] = _S522;

#line 102
            boundaryCount_0 = _S523;

#line 101
        }

#line 99
        index_2 = index_2 + int(1);

#line 99
    }

#line 105
    int _S524 = boundaryCount_0 + int(1);

#line 105
    boundaries_0[boundaryCount_0] = 1.0;

    thread int rootCount_1 = int(0);

#line 107
    int interval_0 = int(0);
    for(;;)
    {

#line 108
        int _S525 = interval_0 + int(1);

#line 108
        if(_S525 < _S524)
        {
        }
        else
        {

#line 108
            break;
        }

#line 109
        float lower_0 = boundaries_0[interval_0];
        float upper_0 = boundaries_0[_S525];

#line 110
        float _S526 = EvaluatePolynomial_0(coefficients_4, degree_2, boundaries_0[interval_0]);

#line 110
        float _S527 = EvaluatePolynomial_0(coefficients_4, degree_2, boundaries_0[_S525]);


        if((abs(_S526)) <= _S521)
        {

#line 114
            AppendUnitRoot_0(roots_3, &rootCount_1, lower_0, fractionTolerance_2);

#line 113
        }


        float _S528 = - _S521;

#line 116
        if(_S526 < _S528)
        {

#line 116
            _S520 = _S527 > _S521;

#line 116
        }
        else
        {

#line 116
            _S520 = false;

#line 116
        }

#line 116
        bool signChange_0;

#line 116
        if(_S520)
        {

#line 116
            signChange_0 = true;

#line 116
        }
        else
        {

#line 117
            if(_S526 > _S521)
            {

#line 117
                signChange_0 = _S527 < _S528;

#line 117
            }
            else
            {

#line 117
                signChange_0 = false;

#line 117
            }

#line 116
        }

        if(signChange_0)
        {

#line 118
            float rootLower_0 = lower_0;

#line 118
            float rootUpper_0 = upper_0;

#line 118
            float rootLowerValue_0 = _S526;

#line 118
            int iteration_0 = int(0);



            for(;;)
            {

#line 122
                if(iteration_0 < int(48))
                {
                }
                else
                {

#line 122
                    break;
                }

#line 123
                float midpoint_0 = (rootLower_0 + rootUpper_0) * 0.5;

#line 123
                float _S529 = EvaluatePolynomial_0(coefficients_4, degree_2, midpoint_0);

                if((rootLowerValue_0 < 0.0) == (_S529 < 0.0))
                {

#line 125
                    rootLower_0 = midpoint_0;

#line 125
                    rootLowerValue_0 = _S529;

#line 125
                }
                else
                {

#line 125
                    rootUpper_0 = midpoint_0;

#line 125
                }

#line 122
                iteration_0 = iteration_0 + int(1);

#line 122
            }

#line 132
            AppendUnitRoot_0(roots_3, &rootCount_1, (rootLower_0 + rootUpper_0) * 0.5, fractionTolerance_2);

#line 118
        }

#line 108
        interval_0 = _S525;

#line 108
    }

#line 108
    float _S530 = EvaluatePolynomial_0(coefficients_4, degree_2, 1.0);

#line 136
    if((abs(_S530)) <= _S521)
    {

#line 137
        AppendUnitRoot_0(roots_3, &rootCount_1, 1.0, fractionTolerance_2);

#line 136
    }


    return rootCount_1;
}

int FindQuadraticRoots_0(const array<float, int(7)> thread* coefficients_5, array<float, int(6)> thread* roots_4)
{

#line 142
    float _S531 = (*coefficients_5)[int(2)];
    float _S532 = abs((*coefficients_5)[int(2)]);

#line 143
    float _S533 = PolynomialValueTolerance_0(coefficients_5, int(2));

#line 143
    if(_S532 <= _S533)
    {

#line 143
        int _S534 = FindLinearRoots_0(coefficients_5, roots_4);
        return _S534;
    }
    thread array<float, int(7)> derivative_0;

#line 146
    int index_3 = int(0);
    for(;;)
    {

#line 147
        if(index_3 < int(7))
        {
        }
        else
        {

#line 147
            break;
        }

#line 147
        derivative_0[index_3] = 0.0;

#line 147
        index_3 = index_3 + int(1);

#line 147
    }
    derivative_0[int(0)] = (*coefficients_5)[int(1)];
    derivative_0[int(1)] = 2.0 * _S531;

#line 149
    thread array<float, int(7)> _S535 = derivative_0;
    thread array<float, int(6)> critical_0;

#line 150
    int _S536 = FindLinearRoots_0(&_S535, &critical_0);

#line 150
    thread array<float, int(6)> _S537 = critical_0;

#line 150
    int _S538 = FindRootsFromCriticalPoints_0(coefficients_5, int(2), &_S537, _S536, roots_4);

    return _S538;
}

int FindCubicRoots_0(const array<float, int(7)> thread* coefficients_6, array<float, int(6)> thread* roots_5)
{

#line 156
    float _S539 = abs((*coefficients_6)[int(3)]);

#line 156
    float _S540 = PolynomialValueTolerance_0(coefficients_6, int(3));

#line 156
    if(_S539 <= _S540)
    {

#line 156
        int _S541 = FindQuadraticRoots_0(coefficients_6, roots_5);
        return _S541;
    }
    thread array<float, int(7)> derivative_1;

#line 159
    int index_4 = int(0);
    for(;;)
    {

#line 160
        if(index_4 < int(7))
        {
        }
        else
        {

#line 160
            break;
        }

#line 160
        derivative_1[index_4] = 0.0;

#line 160
        index_4 = index_4 + int(1);

#line 160
    }

#line 160
    int power_4 = int(1);
    for(;;)
    {

#line 161
        if(power_4 <= int(3))
        {
        }
        else
        {

#line 161
            break;
        }

#line 162
        derivative_1[power_4 - int(1)] = float(power_4) * (*coefficients_6)[power_4];

#line 161
        power_4 = power_4 + int(1);

#line 161
    }

#line 161
    thread array<float, int(7)> _S542 = derivative_1;

    thread array<float, int(6)> critical_1;

#line 163
    int _S543 = FindQuadraticRoots_0(&_S542, &critical_1);

#line 163
    thread array<float, int(6)> _S544 = critical_1;

#line 163
    int _S545 = FindRootsFromCriticalPoints_0(coefficients_6, int(3), &_S544, _S543, roots_5);

    return _S545;
}

int FindQuarticRoots_0(const array<float, int(7)> thread* coefficients_7, array<float, int(6)> thread* roots_6)
{

#line 169
    float _S546 = abs((*coefficients_7)[int(4)]);

#line 169
    float _S547 = PolynomialValueTolerance_0(coefficients_7, int(4));

#line 169
    if(_S546 <= _S547)
    {

#line 169
        int _S548 = FindCubicRoots_0(coefficients_7, roots_6);
        return _S548;
    }
    thread array<float, int(7)> derivative_2;

#line 172
    int index_5 = int(0);
    for(;;)
    {

#line 173
        if(index_5 < int(7))
        {
        }
        else
        {

#line 173
            break;
        }

#line 173
        derivative_2[index_5] = 0.0;

#line 173
        index_5 = index_5 + int(1);

#line 173
    }

#line 173
    int power_5 = int(1);
    for(;;)
    {

#line 174
        if(power_5 <= int(4))
        {
        }
        else
        {

#line 174
            break;
        }

#line 175
        derivative_2[power_5 - int(1)] = float(power_5) * (*coefficients_7)[power_5];

#line 174
        power_5 = power_5 + int(1);

#line 174
    }

#line 174
    thread array<float, int(7)> _S549 = derivative_2;

    thread array<float, int(6)> critical_2;

#line 176
    int _S550 = FindCubicRoots_0(&_S549, &critical_2);

#line 176
    thread array<float, int(6)> _S551 = critical_2;

#line 176
    int _S552 = FindRootsFromCriticalPoints_0(coefficients_7, int(4), &_S551, _S550, roots_6);

    return _S552;
}

int FindQuinticRoots_0(const array<float, int(7)> thread* coefficients_8, array<float, int(6)> thread* roots_7)
{

#line 182
    float _S553 = abs((*coefficients_8)[int(5)]);

#line 182
    float _S554 = PolynomialValueTolerance_0(coefficients_8, int(5));

#line 182
    if(_S553 <= _S554)
    {

#line 182
        int _S555 = FindQuarticRoots_0(coefficients_8, roots_7);
        return _S555;
    }
    thread array<float, int(7)> derivative_3;

#line 185
    int index_6 = int(0);
    for(;;)
    {

#line 186
        if(index_6 < int(7))
        {
        }
        else
        {

#line 186
            break;
        }

#line 186
        derivative_3[index_6] = 0.0;

#line 186
        index_6 = index_6 + int(1);

#line 186
    }

#line 186
    int power_6 = int(1);
    for(;;)
    {

#line 187
        if(power_6 <= int(5))
        {
        }
        else
        {

#line 187
            break;
        }

#line 188
        derivative_3[power_6 - int(1)] = float(power_6) * (*coefficients_8)[power_6];

#line 187
        power_6 = power_6 + int(1);

#line 187
    }

#line 187
    thread array<float, int(7)> _S556 = derivative_3;

    thread array<float, int(6)> critical_3;

#line 189
    int _S557 = FindQuarticRoots_0(&_S556, &critical_3);

#line 189
    thread array<float, int(6)> _S558 = critical_3;

#line 189
    int _S559 = FindRootsFromCriticalPoints_0(coefficients_8, int(5), &_S558, _S557, roots_7);

    return _S559;
}

int FindSexticRoots_0(const array<float, int(7)> thread* coefficients_9, array<float, int(6)> thread* roots_8)
{

#line 195
    float _S560 = abs((*coefficients_9)[int(6)]);

#line 195
    float _S561 = PolynomialValueTolerance_0(coefficients_9, int(6));

#line 195
    if(_S560 <= _S561)
    {

#line 195
        int _S562 = FindQuinticRoots_0(coefficients_9, roots_8);
        return _S562;
    }
    thread array<float, int(7)> derivative_4;

#line 198
    int index_7 = int(0);
    for(;;)
    {

#line 199
        if(index_7 < int(7))
        {
        }
        else
        {

#line 199
            break;
        }

#line 199
        derivative_4[index_7] = 0.0;

#line 199
        index_7 = index_7 + int(1);

#line 199
    }

#line 199
    int power_7 = int(1);
    for(;;)
    {

#line 200
        if(power_7 <= int(6))
        {
        }
        else
        {

#line 200
            break;
        }

#line 201
        derivative_4[power_7 - int(1)] = float(power_7) * (*coefficients_9)[power_7];

#line 200
        power_7 = power_7 + int(1);

#line 200
    }

#line 200
    thread array<float, int(7)> _S563 = derivative_4;

    thread array<float, int(6)> critical_4;

#line 202
    int _S564 = FindQuinticRoots_0(&_S563, &critical_4);

#line 202
    thread array<float, int(6)> _S565 = critical_4;

#line 202
    int _S566 = FindRootsFromCriticalPoints_0(coefficients_9, int(6), &_S565, _S564, roots_8);

    return _S566;
}


#line 13
Vec4Cart_0 AcceptedSegmentPosition_0(const Vec4Cart_0 thread* start_2, const Vec4Cart_0 thread* startVelocity_1, const Vec4Cart_0 thread* finish_1, const Vec4Cart_0 thread* finishVelocity_1, float affineLength_2, float fraction_0)
{

    float s2_2 = fraction_0 * fraction_0;
    float s3_2 = s2_2 * fraction_0;
    float _S567 = 3.0 * s2_2;
    float h10_0 = s3_2 - 2.0 * s2_2 + fraction_0;
    float h01_0 = -2.0 * s3_2 + _S567;
    float h11_0 = s3_2 - s2_2;

#line 21
    Vec4Cart_0 _S568 = ScaleVec4Cart_0(start_2, 2.0 * s3_2 - _S567 + 1.0);

#line 21
    Vec4Cart_0 _S569 = ScaleVec4Cart_0(startVelocity_1, affineLength_2 * h10_0);

#line 21
    thread Vec4Cart_0 _S570 = _S568;

#line 21
    thread Vec4Cart_0 _S571 = _S569;

#line 21
    Vec4Cart_0 _S572 = AddVec4Cart_0(&_S570, &_S571);

#line 21
    Vec4Cart_0 _S573 = ScaleVec4Cart_0(finish_1, h01_0);

#line 21
    Vec4Cart_0 _S574 = ScaleVec4Cart_0(finishVelocity_1, affineLength_2 * h11_0);

#line 21
    thread Vec4Cart_0 _S575 = _S573;

#line 21
    thread Vec4Cart_0 _S576 = _S574;

#line 21
    Vec4Cart_0 _S577 = AddVec4Cart_0(&_S575, &_S576);

#line 21
    thread Vec4Cart_0 _S578 = _S572;

#line 21
    thread Vec4Cart_0 _S579 = _S577;

#line 21
    Vec4Cart_0 _S580 = AddVec4Cart_0(&_S578, &_S579);
    return _S580;
}


#line 29
Vec4Cart_0 AcceptedSegmentTangent_0(const Vec4Cart_0 thread* start_3, const Vec4Cart_0 thread* startVelocity_2, const Vec4Cart_0 thread* finish_2, const Vec4Cart_0 thread* finishVelocity_2, float affineLength_3, float fraction_1)
{

    float s2_3 = fraction_1 * fraction_1;
    float _S581 = 6.0 * fraction_1;
    float _S582 = 3.0 * s2_3;

#line 34
    float dh10_0 = _S582 - 4.0 * fraction_1 + 1.0;
    float dh01_0 = -6.0 * s2_3 + _S581;
    float dh11_0 = _S582 - 2.0 * fraction_1;

#line 36
    Vec4Cart_0 _S583 = ScaleVec4Cart_0(start_3, (6.0 * s2_3 - _S581) / affineLength_3);

#line 36
    Vec4Cart_0 _S584 = ScaleVec4Cart_0(startVelocity_2, dh10_0);

#line 36
    thread Vec4Cart_0 _S585 = _S583;

#line 36
    thread Vec4Cart_0 _S586 = _S584;

#line 36
    Vec4Cart_0 _S587 = AddVec4Cart_0(&_S585, &_S586);

#line 36
    Vec4Cart_0 _S588 = ScaleVec4Cart_0(finish_2, dh01_0 / affineLength_3);

#line 36
    Vec4Cart_0 _S589 = ScaleVec4Cart_0(finishVelocity_2, dh11_0);

#line 36
    thread Vec4Cart_0 _S590 = _S588;

#line 36
    thread Vec4Cart_0 _S591 = _S589;

#line 36
    Vec4Cart_0 _S592 = AddVec4Cart_0(&_S590, &_S591);

#line 36
    thread Vec4Cart_0 _S593 = _S587;

#line 36
    thread Vec4Cart_0 _S594 = _S592;

#line 36
    Vec4Cart_0 _S595 = AddVec4Cart_0(&_S593, &_S594);
    return _S595;
}


#line 243
bool FindSphericalBoundaryEvent_0(const Vec4Cart_0 thread* start_4, const Vec4Cart_0 thread* startVelocity_3, const Vec4Cart_0 thread* finish_3, const Vec4Cart_0 thread* finishVelocity_3, float affineLength_4, float boundaryRadius_0, int radialSense_0, float thread* fraction_2, Vec4Cart_0 thread* eventPosition_0, Vec4Cart_0 thread* eventTangent_0)
{


    *fraction_2 = 0.0;
    *eventPosition_0 = *start_4;
    *eventTangent_0 = *startVelocity_3;

#line 249
    bool _S596;
    if(isfinite(affineLength_4))
    {

#line 250
        _S596 = affineLength_4 > 0.0;

#line 250
    }
    else
    {

#line 250
        _S596 = false;

#line 250
    }

#line 250
    if(_S596)
    {

#line 250
        _S596 = isfinite(boundaryRadius_0);

#line 250
    }
    else
    {

#line 250
        _S596 = false;

#line 250
    }
    if(_S596)
    {

#line 251
        _S596 = boundaryRadius_0 > 0.0;

#line 251
    }
    else
    {

#line 251
        _S596 = false;

#line 251
    }

#line 251
    if(_S596)
    {

#line 252
        if(radialSense_0 == int(0))
        {

#line 252
            _S596 = true;

#line 252
        }
        else
        {

#line 252
            _S596 = radialSense_0 == int(1);

#line 252
        }
        if(_S596)
        {

#line 253
            _S596 = true;

#line 253
        }
        else
        {

#line 253
            _S596 = radialSense_0 == int(-1);

#line 253
        }

#line 251
    }
    else
    {

#line 251
        _S596 = false;

#line 251
    }


    if(_S596)
    {

#line 254
        bool _S597 = IsFiniteAcceptedSegmentValue_0(start_4);

#line 254
        _S596 = _S597;

#line 254
    }
    else
    {

#line 254
        _S596 = false;

#line 254
    }
    if(_S596)
    {

#line 255
        bool _S598 = IsFiniteAcceptedSegmentValue_0(startVelocity_3);

#line 255
        _S596 = _S598;

#line 255
    }
    else
    {

#line 255
        _S596 = false;

#line 255
    }
    if(_S596)
    {

#line 256
        bool _S599 = IsFiniteAcceptedSegmentValue_0(finish_3);

#line 256
        _S596 = _S599;

#line 256
    }
    else
    {

#line 256
        _S596 = false;

#line 256
    }
    if(_S596)
    {

#line 257
        bool _S600 = IsFiniteAcceptedSegmentValue_0(finishVelocity_3);

#line 257
        _S596 = _S600;

#line 257
    }
    else
    {

#line 257
        _S596 = false;

#line 257
    }

#line 250
    if(!_S596)
    {

#line 258
        return false;
    }

#line 258
    Vec4Cart_0 _S601 = ScaleVec4Cart_0(startVelocity_3, affineLength_4 / 3.0);

#line 258
    thread Vec4Cart_0 _S602 = _S601;

#line 258
    Vec4Cart_0 _S603 = AddVec4Cart_0(start_4, &_S602);

#line 258
    Vec4Cart_0 _S604 = ScaleVec4Cart_0(finishVelocity_3, - affineLength_4 / 3.0);

#line 258
    thread Vec4Cart_0 _S605 = _S604;

#line 258
    Vec4Cart_0 _S606 = AddVec4Cart_0(finish_3, &_S605);

#line 269
    float interiorLimit_0 = boundaryRadius_0 * boundaryRadius_0 * (1.0 - 64.0 * PolynomialEpsilon_0());

    float startSquared_0 = start_4->x_2 * start_4->x_2 + start_4->y_0 * start_4->y_0 + start_4->z_0 * start_4->z_0;
    float finishSquared_0 = finish_3->x_2 * finish_3->x_2 + finish_3->y_0 * finish_3->y_0 + finish_3->z_0 * finish_3->z_0;

#line 272
    float _S607 = _S603.x_2;

#line 272
    float _S608 = _S603.y_0;

#line 272
    float _S609 = _S603.z_0;
    float control1Squared_0 = _S607 * _S607 + _S608 * _S608 + _S609 * _S609;

#line 273
    float _S610 = _S606.x_2;

#line 273
    float _S611 = _S606.y_0;

#line 273
    float _S612 = _S606.z_0;

    float control2Squared_0 = _S610 * _S610 + _S611 * _S611 + _S612 * _S612;

    if(isfinite(interiorLimit_0))
    {

#line 277
        _S596 = startSquared_0 < interiorLimit_0;

#line 277
    }
    else
    {

#line 277
        _S596 = false;

#line 277
    }

#line 277
    if(_S596)
    {

#line 277
        _S596 = finishSquared_0 < interiorLimit_0;

#line 277
    }
    else
    {

#line 277
        _S596 = false;

#line 277
    }
    if(_S596)
    {

#line 278
        _S596 = control1Squared_0 < interiorLimit_0;

#line 278
    }
    else
    {

#line 278
        _S596 = false;

#line 278
    }

#line 278
    if(_S596)
    {

#line 278
        _S596 = control2Squared_0 < interiorLimit_0;

#line 278
    }
    else
    {

#line 278
        _S596 = false;

#line 278
    }

#line 277
    if(_S596)
    {
        return false;
    }
    thread array<float, int(7)> coefficients_10;

#line 281
    SphericalSegmentPolynomial_0(start_4, startVelocity_3, finish_3, finishVelocity_3, affineLength_4, boundaryRadius_0, &coefficients_10);

#line 281
    thread array<float, int(7)> _S613 = coefficients_10;


    thread array<float, int(6)> roots_9;

#line 284
    int _S614 = FindSexticRoots_0(&_S613, &roots_9);

    if(_S614 == int(0))
    {

#line 286
        return false;
    }

#line 286
    int root_3 = int(0);

    for(;;)
    {

#line 288
        if(root_3 < _S614)
        {
        }
        else
        {

#line 288
            break;
        }

#line 289
        float candidateFraction_0 = roots_9[root_3];

#line 289
        Vec4Cart_0 _S615 = AcceptedSegmentPosition_0(start_4, startVelocity_3, finish_3, finishVelocity_3, affineLength_4, roots_9[root_3]);

#line 289
        Vec4Cart_0 _S616 = AcceptedSegmentTangent_0(start_4, startVelocity_3, finish_3, finishVelocity_3, affineLength_4, roots_9[root_3]);

#line 289
        float _S617 = _S615.x_2;

#line 289
        float _S618 = _S615.y_0;

#line 289
        float _S619 = _S615.z_0;

#line 296
        float sampledRadius_0 = sqrt(_S617 * _S617 + _S618 * _S618 + _S619 * _S619);


        if(isfinite(sampledRadius_0))
        {

#line 299
            _S596 = sampledRadius_0 > 0.0;

#line 299
        }
        else
        {

#line 299
            _S596 = false;

#line 299
        }

#line 299
        if(!_S596)
        {

#line 299
            root_3 = root_3 + int(1);

#line 288
            continue;
        }

#line 301
        float radialNumerator_0 = _S615.x_2 * _S616.x_2 + _S615.y_0 * _S616.y_0 + _S615.z_0 * _S616.z_0;

        if(!isfinite(radialNumerator_0))
        {

#line 303
            root_3 = root_3 + int(1);

#line 288
            continue;
        }

#line 288
        bool _S620;

#line 304
        if(radialSense_0 == int(1))
        {

#line 304
            _S620 = !(radialNumerator_0 > 0.0);

#line 304
        }
        else
        {

#line 304
            _S620 = false;

#line 304
        }

#line 304
        if(_S620)
        {

#line 305
            root_3 = root_3 + int(1);

#line 288
            continue;
        }

#line 288
        bool _S621;

#line 306
        if(radialSense_0 == int(-1))
        {

#line 306
            _S621 = !(radialNumerator_0 < 0.0);

#line 306
        }
        else
        {

#line 306
            _S621 = false;

#line 306
        }

#line 306
        if(_S621)
        {

#line 307
            root_3 = root_3 + int(1);

#line 288
            continue;
        }

#line 308
        *fraction_2 = candidateFraction_0;
        *eventPosition_0 = _S615;
        *eventTangent_0 = _S616;
        float surfaceScale_0 = boundaryRadius_0 / sampledRadius_0;
        eventPosition_0->x_2 = eventPosition_0->x_2 * surfaceScale_0;
        eventPosition_0->y_0 = eventPosition_0->y_0 * surfaceScale_0;
        eventPosition_0->z_0 = eventPosition_0->z_0 * surfaceScale_0;
        return true;
    }
    return false;
}


#line 243
bool FindSphericalBoundaryEvent_1(const Vec4Cart_0 thread* start_5, const Vec4Cart_0 thread* startVelocity_4, const Vec4Cart_0 thread* finish_4, const Vec4Cart_0 thread* finishVelocity_4, float affineLength_5, float boundaryRadius_1, int radialSense_1, float thread* fraction_3, Vec4Cart_0 thread* eventPosition_1, Vec4Cart_0 thread* eventTangent_1)
{


    *fraction_3 = 0.0;
    *eventPosition_1 = *start_5;
    *eventTangent_1 = *startVelocity_4;

#line 249
    bool _S622;
    if(isfinite(affineLength_5))
    {

#line 250
        _S622 = affineLength_5 > 0.0;

#line 250
    }
    else
    {

#line 250
        _S622 = false;

#line 250
    }

#line 250
    if(_S622)
    {

#line 250
        _S622 = isfinite(boundaryRadius_1);

#line 250
    }
    else
    {

#line 250
        _S622 = false;

#line 250
    }
    if(_S622)
    {

#line 251
        _S622 = boundaryRadius_1 > 0.0;

#line 251
    }
    else
    {

#line 251
        _S622 = false;

#line 251
    }

#line 251
    if(_S622)
    {

#line 252
        if(radialSense_1 == int(0))
        {

#line 252
            _S622 = true;

#line 252
        }
        else
        {

#line 252
            _S622 = radialSense_1 == int(1);

#line 252
        }
        if(_S622)
        {

#line 253
            _S622 = true;

#line 253
        }
        else
        {

#line 253
            _S622 = radialSense_1 == int(-1);

#line 253
        }

#line 251
    }
    else
    {

#line 251
        _S622 = false;

#line 251
    }


    if(_S622)
    {

#line 254
        bool _S623 = IsFiniteAcceptedSegmentValue_0(start_5);

#line 254
        _S622 = _S623;

#line 254
    }
    else
    {

#line 254
        _S622 = false;

#line 254
    }
    if(_S622)
    {

#line 255
        bool _S624 = IsFiniteAcceptedSegmentValue_0(startVelocity_4);

#line 255
        _S622 = _S624;

#line 255
    }
    else
    {

#line 255
        _S622 = false;

#line 255
    }
    if(_S622)
    {

#line 256
        bool _S625 = IsFiniteAcceptedSegmentValue_0(finish_4);

#line 256
        _S622 = _S625;

#line 256
    }
    else
    {

#line 256
        _S622 = false;

#line 256
    }
    if(_S622)
    {

#line 257
        bool _S626 = IsFiniteAcceptedSegmentValue_0(finishVelocity_4);

#line 257
        _S622 = _S626;

#line 257
    }
    else
    {

#line 257
        _S622 = false;

#line 257
    }

#line 250
    if(!_S622)
    {

#line 258
        return false;
    }

#line 258
    Vec4Cart_0 _S627 = ScaleVec4Cart_0(startVelocity_4, affineLength_5 / 3.0);

#line 258
    thread Vec4Cart_0 _S628 = _S627;

#line 258
    Vec4Cart_0 _S629 = AddVec4Cart_0(start_5, &_S628);

#line 258
    Vec4Cart_0 _S630 = ScaleVec4Cart_0(finishVelocity_4, - affineLength_5 / 3.0);

#line 258
    thread Vec4Cart_0 _S631 = _S630;

#line 258
    Vec4Cart_0 _S632 = AddVec4Cart_0(finish_4, &_S631);

#line 269
    float interiorLimit_1 = boundaryRadius_1 * boundaryRadius_1 * (1.0 - 64.0 * PolynomialEpsilon_0());

    float startSquared_1 = start_5->x_2 * start_5->x_2 + start_5->y_0 * start_5->y_0 + start_5->z_0 * start_5->z_0;
    float finishSquared_1 = finish_4->x_2 * finish_4->x_2 + finish_4->y_0 * finish_4->y_0 + finish_4->z_0 * finish_4->z_0;

#line 272
    float _S633 = _S629.x_2;

#line 272
    float _S634 = _S629.y_0;

#line 272
    float _S635 = _S629.z_0;
    float control1Squared_1 = _S633 * _S633 + _S634 * _S634 + _S635 * _S635;

#line 273
    float _S636 = _S632.x_2;

#line 273
    float _S637 = _S632.y_0;

#line 273
    float _S638 = _S632.z_0;

    float control2Squared_1 = _S636 * _S636 + _S637 * _S637 + _S638 * _S638;

    if(isfinite(interiorLimit_1))
    {

#line 277
        _S622 = startSquared_1 < interiorLimit_1;

#line 277
    }
    else
    {

#line 277
        _S622 = false;

#line 277
    }

#line 277
    if(_S622)
    {

#line 277
        _S622 = finishSquared_1 < interiorLimit_1;

#line 277
    }
    else
    {

#line 277
        _S622 = false;

#line 277
    }
    if(_S622)
    {

#line 278
        _S622 = control1Squared_1 < interiorLimit_1;

#line 278
    }
    else
    {

#line 278
        _S622 = false;

#line 278
    }

#line 278
    if(_S622)
    {

#line 278
        _S622 = control2Squared_1 < interiorLimit_1;

#line 278
    }
    else
    {

#line 278
        _S622 = false;

#line 278
    }

#line 277
    if(_S622)
    {
        return false;
    }
    thread array<float, int(7)> coefficients_11;

#line 281
    SphericalSegmentPolynomial_0(start_5, startVelocity_4, finish_4, finishVelocity_4, affineLength_5, boundaryRadius_1, &coefficients_11);

#line 281
    thread array<float, int(7)> _S639 = coefficients_11;


    thread array<float, int(6)> roots_10;

#line 284
    int _S640 = FindSexticRoots_0(&_S639, &roots_10);

    if(_S640 == int(0))
    {

#line 286
        return false;
    }

#line 286
    int root_4 = int(0);

    for(;;)
    {

#line 288
        if(root_4 < _S640)
        {
        }
        else
        {

#line 288
            break;
        }

#line 289
        float candidateFraction_1 = roots_10[root_4];

#line 289
        Vec4Cart_0 _S641 = AcceptedSegmentPosition_0(start_5, startVelocity_4, finish_4, finishVelocity_4, affineLength_5, roots_10[root_4]);

#line 289
        Vec4Cart_0 _S642 = AcceptedSegmentTangent_0(start_5, startVelocity_4, finish_4, finishVelocity_4, affineLength_5, roots_10[root_4]);

#line 289
        float _S643 = _S641.x_2;

#line 289
        float _S644 = _S641.y_0;

#line 289
        float _S645 = _S641.z_0;

#line 296
        float sampledRadius_1 = sqrt(_S643 * _S643 + _S644 * _S644 + _S645 * _S645);


        if(isfinite(sampledRadius_1))
        {

#line 299
            _S622 = sampledRadius_1 > 0.0;

#line 299
        }
        else
        {

#line 299
            _S622 = false;

#line 299
        }

#line 299
        if(!_S622)
        {

#line 299
            root_4 = root_4 + int(1);

#line 288
            continue;
        }

#line 301
        float radialNumerator_1 = _S641.x_2 * _S642.x_2 + _S641.y_0 * _S642.y_0 + _S641.z_0 * _S642.z_0;

        if(!isfinite(radialNumerator_1))
        {

#line 303
            root_4 = root_4 + int(1);

#line 288
            continue;
        }

#line 288
        bool _S646;

#line 304
        if(radialSense_1 == int(1))
        {

#line 304
            _S646 = !(radialNumerator_1 > 0.0);

#line 304
        }
        else
        {

#line 304
            _S646 = false;

#line 304
        }

#line 304
        if(_S646)
        {

#line 305
            root_4 = root_4 + int(1);

#line 288
            continue;
        }

#line 288
        bool _S647;

#line 306
        if(radialSense_1 == int(-1))
        {

#line 306
            _S647 = !(radialNumerator_1 < 0.0);

#line 306
        }
        else
        {

#line 306
            _S647 = false;

#line 306
        }

#line 306
        if(_S647)
        {

#line 307
            root_4 = root_4 + int(1);

#line 288
            continue;
        }

#line 308
        *fraction_3 = candidateFraction_1;
        *eventPosition_1 = _S641;
        *eventTangent_1 = _S642;
        float surfaceScale_1 = boundaryRadius_1 / sampledRadius_1;
        eventPosition_1->x_2 = eventPosition_1->x_2 * surfaceScale_1;
        eventPosition_1->y_0 = eventPosition_1->y_0 * surfaceScale_1;
        eventPosition_1->z_0 = eventPosition_1->z_0 * surfaceScale_1;
        return true;
    }
    return false;
}


#line 87 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_ellis_topology.slang"
void AttemptEllisTwoSheet_0(EllisTopologyTraceResult_0 thread* result_9, float thread* adaptiveStep_1, Vec4Cart_0 thread* compensationX_1, Vec4Cart_0 thread* compensationU_1, uint thread* accepted_2, float thread* lastStep_1, uint thread* failure_1, float stepScale_1, float minStep_1, float maxStep_1, float escapeRadius_2, float b0_8, int wormShape_5)
{

#line 93
    float oppositeRadius_1 = EllisOppositeEscapeRadius_0(b0_8, escapeRadius_2);
    result_9->attempts_0 = result_9->attempts_0 + int(1);
    float radius_0 = sqrt((&(&result_9->state_0)->x_5)->x_2 * (&(&result_9->state_0)->x_5)->x_2 + (&(&result_9->state_0)->x_5)->y_0 * (&(&result_9->state_0)->x_5)->y_0 + (&(&result_9->state_0)->x_5)->z_0 * (&(&result_9->state_0)->x_5)->z_0);



    float radialNumerator_2 = (&(&result_9->state_0)->x_5)->x_2 * (&(&result_9->state_0)->u_1)->x_2 + (&(&result_9->state_0)->x_5)->y_0 * (&(&result_9->state_0)->u_1)->y_0 + (&(&result_9->state_0)->x_5)->z_0 * (&(&result_9->state_0)->u_1)->z_0;

#line 99
    bool _S648;

    if(radius_0 >= escapeRadius_2)
    {

#line 101
        _S648 = radialNumerator_2 > 0.0;

#line 101
    }
    else
    {

#line 101
        _S648 = false;

#line 101
    }

#line 101
    if(_S648)
    {

#line 102
        result_9->termination_0 = int(1);
        return;
    }
    if(radius_0 <= oppositeRadius_1)
    {

#line 105
        _S648 = radialNumerator_2 < 0.0;

#line 105
    }
    else
    {

#line 105
        _S648 = false;

#line 105
    }

#line 105
    if(_S648)
    {

#line 106
        result_9->termination_0 = int(2);
        return;
    }


    float _S649 = min(clamp(stepScale_1 * max(radius_0, b0_8), minStep_1, maxStep_1), *adaptiveStep_1);



    float _S650 = AdaptiveIntegrationTolerance_0();

#line 115
    float _S651 = AdaptiveNullTolerance_0();

#line 115
    thread GeodesicStateCart_0 _S652 = result_9->state_0;

#line 115
    AdaptiveGeodesicStep_0 _S653 = AttemptAdaptiveGeodesicStep_0(&_S652, _S649, false, true, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, b0_8, wormShape_5, _S650, _S651);
    if(!_S653.accepted_1)
    {

#line 117
        if(AdaptiveAttemptAlreadyAtMinimumStep_0(_S649, minStep_1))
        {

#line 118
            *failure_1 = 3U;
            return;
        }
        *adaptiveStep_1 = max(minStep_1, _S649 * _S653.stepFactor_0);
        return;
    }
    *adaptiveStep_1 = clamp(_S649 * _S653.stepFactor_0, minStep_1, maxStep_1);

#line 124
    thread Vec4Cart_0 _S654 = (&result_9->state_0)->x_5;

#line 124
    thread Vec4Cart_0 _S655 = (&result_9->state_0)->u_1;

#line 124
    thread Vec4Cart_0 _S656 = _S653.candidate_0.x_5;

#line 124
    thread Vec4Cart_0 _S657 = _S653.candidate_0.u_1;

#line 136
    thread float observerFraction_0;
    thread Vec4Cart_0 observerPosition_0;
    thread Vec4Cart_0 observerTangent_0;

#line 138
    bool _S658 = FindSphericalBoundaryEvent_0(&_S654, &_S655, &_S656, &_S657, _S649, escapeRadius_2, int(1), &observerFraction_0, &observerPosition_0, &observerTangent_0);

#line 138
    thread Vec4Cart_0 _S659 = (&result_9->state_0)->x_5;

#line 138
    thread Vec4Cart_0 _S660 = (&result_9->state_0)->u_1;

#line 138
    thread Vec4Cart_0 _S661 = _S653.candidate_0.x_5;

#line 138
    thread Vec4Cart_0 _S662 = _S653.candidate_0.u_1;

#line 143
    thread float oppositeFraction_0;
    thread Vec4Cart_0 oppositePosition_0;
    thread Vec4Cart_0 oppositeTangent_0;

#line 145
    bool _S663 = FindSphericalBoundaryEvent_0(&_S659, &_S660, &_S661, &_S662, _S649, oppositeRadius_1, int(-1), &oppositeFraction_0, &oppositePosition_0, &oppositeTangent_0);

#line 151
    if(_S658)
    {

#line 151
        if(!_S663)
        {

#line 151
            _S648 = true;

#line 151
        }
        else
        {

#line 151
            _S648 = observerFraction_0 <= oppositeFraction_0;

#line 151
        }

#line 151
    }
    else
    {

#line 151
        _S648 = false;

#line 151
    }

#line 151
    if(_S648)
    {

#line 152
        (&result_9->state_0)->x_5 = observerPosition_0;
        (&result_9->state_0)->u_1 = observerTangent_0;
        *lastStep_1 = _S649 * observerFraction_0;
        *accepted_2 = *accepted_2 + 1U;
        result_9->affineLength_0 = result_9->affineLength_0 + _S649 * observerFraction_0;
        result_9->termination_0 = int(1);
        return;
    }
    if(_S663)
    {

#line 161
        (&result_9->state_0)->x_5 = oppositePosition_0;
        (&result_9->state_0)->u_1 = oppositeTangent_0;
        *lastStep_1 = _S649 * oppositeFraction_0;
        *accepted_2 = *accepted_2 + 1U;
        result_9->affineLength_0 = result_9->affineLength_0 + _S649 * oppositeFraction_0;
        result_9->termination_0 = int(2);
        return;
    }

    result_9->state_0 = _S653.candidate_0;
    *lastStep_1 = _S649;
    *accepted_2 = *accepted_2 + 1U;
    result_9->affineLength_0 = result_9->affineLength_0 + _S649;
    return;
}


#line 176
void AdvanceEllisTwoSheet_0(EllisTopologyContinuation_0 thread* continuation_2, int maxAttempts_1, float stepScale_2, float minStep_2, float maxStep_2, float escapeRadius_3, float b0_9, int wormShape_6)
{

    bool _S664;
    if((continuation_2->failure_0) != 0U)
    {

#line 180
        _S664 = true;

#line 180
    }
    else
    {

#line 180
        _S664 = ((&continuation_2->trace_0)->termination_0) != int(0);

#line 180
    }

#line 180
    if(_S664)
    {

#line 181
        return;
    }

#line 182
    if(((&continuation_2->trace_0)->attempts_0) >= maxAttempts_1)
    {

#line 183
        continuation_2->failure_0 = 2U;
        return;
    }
    AttemptEllisTwoSheet_0(&continuation_2->trace_0, &continuation_2->adaptiveStep_0, &continuation_2->compensationX_0, &continuation_2->compensationU_0, &continuation_2->accepted_0, &continuation_2->lastStep_0, &continuation_2->failure_0, stepScale_2, minStep_2, maxStep_2, escapeRadius_3, b0_9, wormShape_6);

#line 191
    if(((&continuation_2->trace_0)->termination_0) == int(0))
    {

#line 191
        _S664 = (continuation_2->failure_0) == 0U;

#line 191
    }
    else
    {

#line 191
        _S664 = false;

#line 191
    }
    if(_S664)
    {

#line 192
        _S664 = ((&continuation_2->trace_0)->attempts_0) >= maxAttempts_1;

#line 192
    }
    else
    {

#line 192
        _S664 = false;

#line 192
    }

#line 191
    if(_S664)
    {

        continuation_2->failure_0 = 2U;

#line 191
    }



    return;
}


#line 487 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_integrator.slang"
int AdaptiveCoarseSubstepsForAttempt_0(float step_4, float minimumStep_1)
{

#line 488
    if(step_4 != minimumStep_1)
    {

#line 488
        return int(1);
    }


    return int(2);
}


#line 581
AdaptiveGeodesicStep_0 AttemptAdaptiveGeodesicStepAtFloor_0(const GeodesicStateCart_0 thread* state_14, float step_5, float minimumStep_2, bool isWarp_4, bool isWormhole_4, float M_10, float a_14, float Q_10, float Lambda_10, float warpVs_4, float warpSigma_4, float warpR_4, float throatB0_4, int wormShape_7, float integrationTolerance_2, float nullTolerance_2)
{


    AdaptiveGeodesicStep_0 _S665 = AttemptAdaptiveGeodesicStepWithResolution_0(state_14, step_5, isWarp_4, isWormhole_4, M_10, a_14, Q_10, Lambda_10, warpVs_4, warpSigma_4, warpR_4, throatB0_4, wormShape_7, integrationTolerance_2, nullTolerance_2, AdaptiveCoarseSubstepsForAttempt_0(step_5, minimumStep_2));

    return _S665;
}


#line 321 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_trace_event.slang"
bool FindSphericalCaptureEvent_0(const Vec4Cart_0 thread* start_6, const Vec4Cart_0 thread* startVelocity_5, const Vec4Cart_0 thread* finish_5, const Vec4Cart_0 thread* finishVelocity_5, float affineLength_6, float captureRadius_1, float thread* fraction_4, Vec4Cart_0 thread* eventPosition_2, Vec4Cart_0 thread* eventTangent_2)
{

    bool _S666 = FindSphericalBoundaryEvent_1(start_6, startVelocity_5, finish_5, finishVelocity_5, affineLength_6, captureRadius_1, int(0), fraction_4, eventPosition_2, eventTangent_2);
    return _S666;
}


#line 286 "/home/astra/.project/Sirius/out/trace-jacobi-shared-segment-rhs-screen/inputs/candidate/trace.slang"
float VolumeScaleHeight_0(float r_3, float innerRadius_0, float hOverR_0, float hPower_0)
{

#line 286
    bool _S667;
    if(r_3 <= 0.0)
    {

#line 287
        _S667 = true;

#line 287
    }
    else
    {

#line 287
        _S667 = innerRadius_0 <= 0.0;

#line 287
    }

#line 287
    if(_S667)
    {

#line 287
        return 0.0;
    }
    return clamp(hOverR_0 * pow(r_3 / innerRadius_0, hPower_0), 0.00999999977648258, 0.5) * r_3;
}

bool IsInsideVolumeDisk_0(float cylindricalR_0, float z_3, float innerRadius_1, float outerRadius_0, float hOverR_1, float hPower_1)
{

#line 293
    bool _S668;
    if(cylindricalR_0 < innerRadius_1)
    {

#line 294
        _S668 = true;

#line 294
    }
    else
    {

#line 294
        _S668 = cylindricalR_0 > outerRadius_0;

#line 294
    }

#line 294
    if(_S668)
    {

#line 294
        return false;
    }
    return (abs(z_3)) <= (3.0 * VolumeScaleHeight_0(cylindricalR_0, innerRadius_1, hOverR_1, hPower_1));
}


#line 39 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_disk.slang"
float TruncatedGaussianOpacityDensity_0(float tauRef_0, float radius_1, float referenceRadius_0, float scaleHeight_0, float height_0)
{
    bool _S669;
    if(!isfinite(tauRef_0))
    {

#line 42
        _S669 = true;

#line 42
    }
    else
    {

#line 42
        _S669 = tauRef_0 < 0.0;

#line 42
    }

#line 42
    if(_S669)
    {

#line 42
        _S669 = true;

#line 42
    }
    else
    {

#line 42
        _S669 = !isfinite(radius_1);

#line 42
    }

#line 42
    if(_S669)
    {

#line 42
        _S669 = true;

#line 42
    }
    else
    {

#line 42
        _S669 = radius_1 <= 0.0;

#line 42
    }

#line 42
    if(_S669)
    {

#line 42
        _S669 = true;

#line 42
    }
    else
    {

#line 42
        _S669 = !isfinite(referenceRadius_0);

#line 42
    }
    if(_S669)
    {

#line 43
        _S669 = true;

#line 43
    }
    else
    {

#line 43
        _S669 = referenceRadius_0 <= 0.0;

#line 43
    }

#line 43
    if(_S669)
    {

#line 43
        _S669 = true;

#line 43
    }
    else
    {

#line 43
        _S669 = !isfinite(scaleHeight_0);

#line 43
    }
    if(_S669)
    {

#line 44
        _S669 = true;

#line 44
    }
    else
    {

#line 44
        _S669 = scaleHeight_0 <= 0.0;

#line 44
    }

#line 44
    if(_S669)
    {

#line 44
        _S669 = true;

#line 44
    }
    else
    {

#line 44
        _S669 = !isfinite(height_0);

#line 44
    }

#line 42
    if(_S669)
    {

        return 0.0;
    }
    float zOverH_0 = height_0 / scaleHeight_0;
    if(!isfinite(zOverH_0))
    {

#line 48
        _S669 = true;

#line 48
    }
    else
    {

#line 48
        _S669 = (abs(zOverH_0)) > 3.0;

#line 48
    }

#line 48
    if(_S669)
    {

#line 49
        return 0.0;
    }

    float density_0 = tauRef_0 * pow(radius_1 / referenceRadius_0, -1.5) * exp(-0.5 * zOverH_0 * zOverH_0) / (2.49986100196838379 * scaleHeight_0);

    if(isfinite(density_0))
    {

#line 54
        _S669 = density_0 >= 0.0;

#line 54
    }
    else
    {

#line 54
        _S669 = false;

#line 54
    }

#line 54
    float _S670;

#line 54
    if(_S669)
    {

#line 54
        _S670 = density_0;

#line 54
    }
    else
    {

#line 54
        _S670 = 0.0;

#line 54
    }

#line 54
    return _S670;
}


#line 243 "/home/astra/.project/Sirius/out/trace-jacobi-shared-segment-rhs-screen/inputs/candidate/trace.slang"
float TurbulenceFade_0(float t_2)
{

#line 244
    return t_2 * t_2 * t_2 * (t_2 * (t_2 * 6.0 - 15.0) + 10.0);
}


#line 221
uint TurbulenceHash_0(uint x_6, uint y_3, uint z_4, uint seed_0)
{


    uint h_0 = ((seed_0 ^ (x_6 * 374761393U)) ^ (y_3 * 668265263U)) ^ (z_4 * 1013904223U);
    uint h_1 = ((h_0 >> 16U) ^ h_0) * 2654435761U;

    return ((h_1 >> 16U) ^ h_1) * 2654435761U;
}

float3 TurbulenceGradient_0(uint h_2)
{

#line 240
    return gradients_0[h_2 % 12U];
}


#line 247
float TurbulencePerlin_0(float3 p_0, uint seed_1)
{

#line 248
    int3 _S671 = int3(floor(p_0));
    float3 f_0 = p_0 - float3(_S671);
    float _S672 = TurbulenceFade_0(f_0.x);

#line 250
    float _S673 = TurbulenceFade_0(f_0.y);

#line 250
    float _S674 = TurbulenceFade_0(f_0.z);
    thread array<float, int(8)> corner_0;

#line 251
    int z_5 = int(0);
    for(;;)
    {

#line 252
        if(z_5 <= int(1))
        {
        }
        else
        {

#line 252
            break;
        }

#line 252
        int y_4 = int(0);
        for(;;)
        {

#line 253
            if(y_4 <= int(1))
            {
            }
            else
            {

#line 253
                break;
            }

#line 253
            int x_7 = int(0);
            for(;;)
            {

#line 254
                if(x_7 <= int(1))
                {
                }
                else
                {

#line 254
                    break;
                }

                corner_0[z_5 * int(4) + y_4 * int(2) + x_7] = dot(TurbulenceGradient_0(TurbulenceHash_0(uint(_S671.x + x_7), uint(_S671.y + y_4), uint(_S671.z + z_5), seed_1)), f_0 - float3(float(x_7), float(y_4), float(z_5)));

#line 254
                x_7 = x_7 + int(1);

#line 254
            }

#line 253
            y_4 = y_4 + int(1);

#line 253
        }

#line 252
        z_5 = z_5 + int(1);

#line 252
    }

#line 268
    return mix(mix(mix(corner_0[int(0)], corner_0[int(1)], _S672), mix(corner_0[int(2)], corner_0[int(3)], _S672), _S673), mix(mix(corner_0[int(4)], corner_0[int(5)], _S672), mix(corner_0[int(6)], corner_0[int(7)], _S672), _S673), _S674);
}

float TurbulenceDensityFactor_0(float3 p_1)
{

#line 271
    uint octave_0 = 0U;

#line 271
    float octaveAmplitude_0 = 1.0;

#line 271
    float frequency_0 = 0.20000000298023224;

#line 271
    float value_2 = 0.0;

#line 271
    float totalAmplitude_0 = 0.0;

#line 276
    for(;;)
    {

#line 276
        if(octave_0 < 6U)
        {
        }
        else
        {

#line 276
            break;
        }

#line 277
        float value_3 = value_2 + octaveAmplitude_0 * TurbulencePerlin_0(p_1 * float3(frequency_0) , 12345U + octave_0);

        float totalAmplitude_1 = totalAmplitude_0 + octaveAmplitude_0;
        float octaveAmplitude_1 = octaveAmplitude_0 * 0.5;
        float frequency_1 = frequency_0 * 2.0;

#line 276
        octave_0 = octave_0 + 1U;

#line 276
        octaveAmplitude_0 = octaveAmplitude_1;

#line 276
        frequency_0 = frequency_1;

#line 276
        value_2 = value_3;

#line 276
        totalAmplitude_0 = totalAmplitude_1;

#line 276
    }

#line 283
    return max(1.0 + 0.30000001192092896 * value_2 / totalAmplitude_0, 0.00999999977648258);
}


#line 231 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_disk.slang"
float ShakuraSunyaevFluxShape_0(float r_4, float rInner_0)
{

#line 231
    bool _S675;
    if(r_4 <= rInner_0)
    {

#line 232
        _S675 = true;

#line 232
    }
    else
    {

#line 232
        _S675 = rInner_0 <= 0.0;

#line 232
    }

#line 232
    if(_S675)
    {

#line 232
        return 0.0;
    }

#line 233
    float ratio_0 = rInner_0 / r_4;
    return ratio_0 * ratio_0 * ratio_0 * max(1.0 - sqrt(ratio_0), 0.0);
}

float ShakuraSunyaevTemperature_0(float temperatureScale_0, float r_5, float rInner_1)
{
    if(temperatureScale_0 <= 0.0)
    {

#line 239
        return 0.0;
    }

#line 240
    float flux_0 = ShakuraSunyaevFluxShape_0(r_5, rInner_1);
    float reference_2 = ShakuraSunyaevFluxShape_0(1.5 * rInner_1, rInner_1);

#line 241
    bool _S676;

    if(flux_0 <= 0.0)
    {

#line 243
        _S676 = true;

#line 243
    }
    else
    {

#line 243
        _S676 = reference_2 <= 0.0;

#line 243
    }

#line 243
    if(_S676)
    {

#line 243
        return 0.0;
    }

#line 244
    return temperatureScale_0 * pow(flux_0 / reference_2, 0.25);
}


#line 158
float PtSpecificEnergy_0(float r_6, float a_15)
{

#line 159
    float sqrtR_0 = sqrt(r_6);
    float _S677 = r_6 * r_6;
    return (_S677 - 2.0 * r_6 + a_15 * sqrtR_0) / (r_6 * sqrt(max(_S677 - 3.0 * r_6 + 2.0 * a_15 * sqrtR_0, 9.99999968265522539e-21)));
}

float PtSpecificAngularMomentum_0(float r_7, float a_16)
{

#line 165
    float sqrtR_1 = sqrt(r_7);
    float _S678 = r_7 * r_7;

#line 166
    float _S679 = 2.0 * a_16 * sqrtR_1;
    return (_S678 - _S679 + a_16 * a_16) / (sqrtR_1 * sqrt(max(_S678 - 3.0 * r_7 + _S679, 9.99999968265522539e-21)));
}

float PtAngularVelocity_0(float r_8, float a_17)
{

#line 171
    return 1.0 / (pow(r_8, 1.5) + a_17);
}

float PtAngularMomentumDerivative_0(float r_9, float a_18)
{

#line 175
    float sqrtR_2 = sqrt(r_9);
    float _S680 = r_9 * r_9;

#line 176
    float _S681 = 2.0 * a_18 * sqrtR_2;

#line 176
    float N_0 = _S680 - _S681 + a_18 * a_18;
    float radicand_0 = _S680 - 3.0 * r_9 + _S681;
    if(radicand_0 <= 9.999999960041972e-13)
    {

#line 178
        return 0.0;
    }

#line 179
    float aOverSqrtR_0 = a_18 / sqrtR_2;
    float _S682 = 2.0 * r_9;



    return (_S682 - aOverSqrtR_0 - N_0 * (0.5 / r_9 + 0.5 * (_S682 - 3.0 + aOverSqrtR_0) / radicand_0)) / (sqrtR_2 * sqrt(radicand_0));
}

float PtAngularVelocityDerivative_0(float r_10, float a_19)
{

#line 188
    float denominator_0 = pow(r_10, 1.5) + a_19;
    return -3.0 * sqrt(r_10) / (2.0 * denominator_0 * denominator_0);
}

float FullPageThorneFluxShape_0(float r_11, float rInner_2, float aStar_0)
{

#line 193
    if(r_11 <= (rInner_2 * 1.0010000467300415))
    {

#line 193
        return 0.0;
    }

#line 194
    float energy_0 = PtSpecificEnergy_0(r_11, aStar_0);
    float angularMomentum_0 = PtSpecificAngularMomentum_0(r_11, aStar_0);
    float omega_0 = PtAngularVelocity_0(r_11, aStar_0);
    float halfWidth_0 = 0.5 * (r_11 - rInner_2);
    float _S683 = 0.5 * (r_11 + rInner_2);

#line 198
    int i_9 = int(0);

#line 198
    float integral_0 = 0.0;

    for(;;)
    {

#line 200
        if(i_9 < int(16))
        {
        }
        else
        {

#line 200
            break;
        }

#line 201
        float ri_0 = _S683 + halfWidth_0 * PT_GL_NODES_0[i_9];



        float integral_1 = integral_0 + PT_GL_WEIGHTS_0[i_9] * (PtSpecificEnergy_0(ri_0, aStar_0) - PtAngularVelocity_0(ri_0, aStar_0) * PtSpecificAngularMomentum_0(ri_0, aStar_0)) * PtAngularMomentumDerivative_0(ri_0, aStar_0);

#line 200
        i_9 = i_9 + int(1);

#line 200
        integral_0 = integral_1;

#line 200
    }

#line 208
    float integral_2 = integral_0 * halfWidth_0;
    float invariant_0 = energy_0 - omega_0 * angularMomentum_0;
    if((abs(invariant_0)) < 9.999999960041972e-13)
    {

#line 210
        return 0.0;
    }

    return max(- PtAngularVelocityDerivative_0(r_11, aStar_0) * integral_2 / (invariant_0 * invariant_0) / r_11, 0.0);
}

float FullPageThorneTemperature_0(float temperatureScale_1, float r_12, float rInner_3, float M_11, float aStar_1)
{

    float innerInM_0 = rInner_3 / M_11;
    float flux_1 = FullPageThorneFluxShape_0(r_12 / M_11, innerInM_0, aStar_1);
    float reference_3 = FullPageThorneFluxShape_0(1.5 * innerInM_0, innerInM_0, aStar_1);

#line 221
    bool _S684;

    if(flux_1 <= 0.0)
    {

#line 223
        _S684 = true;

#line 223
    }
    else
    {

#line 223
        _S684 = reference_3 <= 0.0;

#line 223
    }

#line 223
    if(_S684)
    {

#line 223
        return 0.0;
    }

#line 224
    return temperatureScale_1 * pow(flux_1 / reference_3, 0.25);
}


#line 299 "/home/astra/.project/Sirius/out/trace-jacobi-shared-segment-rhs-screen/inputs/candidate/trace.slang"
float VolumeDiskTemperature_0(float cylindricalR_1, float z_6, float innerRadius_2, float innerTemperature_0, int temperatureModel_0, float M_12, float aStar_2, float hOverR_2, float hPower_2)
{
    float midplane_0;

    if(temperatureModel_0 == int(1))
    {

#line 303
        midplane_0 = ShakuraSunyaevTemperature_0(innerTemperature_0, cylindricalR_1, innerRadius_2);

#line 303
    }
    else
    {

#line 303
        midplane_0 = FullPageThorneTemperature_0(innerTemperature_0, cylindricalR_1, innerRadius_2, M_12, aStar_2);

#line 303
    }



    return max(midplane_0, 0.0);
}


#line 250 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_disk.slang"
bool KerrPhotonKillingQuantitiesCart_0(const Vec4Cart_0 thread* position_0, const Vec4Cart_0 thread* pastVelocity_0, const array<array<float, int(4)>, int(4)> thread* metric_8, float thread* energy_1, float thread* angularMomentum_1)
{

    array<float, int(4)> _S685 = { - pastVelocity_0->t_0, - pastVelocity_0->x_2, - pastVelocity_0->y_0, - pastVelocity_0->z_0 };


    thread array<float, int(4)> p_2;

#line 256
    p_2[int(0)] = 0.0;

#line 256
    p_2[int(1)] = 0.0;

#line 256
    p_2[int(2)] = 0.0;

#line 256
    p_2[int(3)] = 0.0;

#line 256
    int mu_8 = int(0);
    for(;;)
    {

#line 257
        if(mu_8 < int(4))
        {
        }
        else
        {

#line 257
            break;
        }

#line 257
        int nu_11 = int(0);
        for(;;)
        {

#line 258
            if(nu_11 < int(4))
            {
            }
            else
            {

#line 258
                break;
            }

#line 259
            p_2[mu_8] = p_2[mu_8] + (*metric_8)[mu_8][nu_11] * _S685[nu_11];

#line 258
            nu_11 = nu_11 + int(1);

#line 258
        }

#line 257
        mu_8 = mu_8 + int(1);

#line 257
    }

#line 262
    *energy_1 = - p_2[int(0)];
    *angularMomentum_1 = - position_0->y_0 * p_2[int(1)] + position_0->x_2 * p_2[int(2)];

#line 263
    bool _S686;
    if(isfinite(*energy_1))
    {

#line 264
        _S686 = (*energy_1) > 0.0;

#line 264
    }
    else
    {

#line 264
        _S686 = false;

#line 264
    }

#line 264
    if(_S686)
    {

#line 264
        _S686 = isfinite(*angularMomentum_1);

#line 264
    }
    else
    {

#line 264
        _S686 = false;

#line 264
    }

#line 264
    return _S686;
}


#line 272
bool TryKerrStationaryFrameFrequencyTransfer_0(float observerFrequency_0, float energy_2, float angularMomentum_2, float M_13, float a_20, float r_13, float cosTheta_0, float angularVelocity_0, float4 thread* transfer_0)
{

    *transfer_0 = float4(0.0) ;

#line 275
    bool _S687;
    if(!isfinite(observerFrequency_0))
    {

#line 276
        _S687 = true;

#line 276
    }
    else
    {

#line 276
        _S687 = observerFrequency_0 <= 0.0;

#line 276
    }

#line 276
    if(_S687)
    {

#line 276
        _S687 = true;

#line 276
    }
    else
    {

#line 276
        _S687 = !isfinite(energy_2);

#line 276
    }
    if(_S687)
    {

#line 277
        _S687 = true;

#line 277
    }
    else
    {

#line 277
        _S687 = energy_2 <= 0.0;

#line 277
    }

#line 277
    if(_S687)
    {

#line 277
        _S687 = true;

#line 277
    }
    else
    {

#line 277
        _S687 = !isfinite(angularMomentum_2);

#line 277
    }

#line 277
    if(_S687)
    {

#line 277
        _S687 = true;

#line 277
    }
    else
    {

#line 277
        _S687 = !isfinite(M_13);

#line 277
    }
    if(_S687)
    {

#line 278
        _S687 = true;

#line 278
    }
    else
    {

#line 278
        _S687 = M_13 <= 0.0;

#line 278
    }

#line 278
    if(_S687)
    {

#line 278
        _S687 = true;

#line 278
    }
    else
    {

#line 278
        _S687 = !isfinite(a_20);

#line 278
    }

#line 278
    if(_S687)
    {

#line 278
        _S687 = true;

#line 278
    }
    else
    {

#line 278
        _S687 = (abs(a_20)) > M_13;

#line 278
    }

#line 278
    if(_S687)
    {

#line 278
        _S687 = true;

#line 278
    }
    else
    {

#line 278
        _S687 = !isfinite(r_13);

#line 278
    }
    if(_S687)
    {

#line 279
        _S687 = true;

#line 279
    }
    else
    {

#line 279
        _S687 = r_13 <= 0.0;

#line 279
    }

#line 279
    if(_S687)
    {

#line 279
        _S687 = true;

#line 279
    }
    else
    {

#line 279
        _S687 = !isfinite(cosTheta_0);

#line 279
    }

#line 279
    if(_S687)
    {

#line 279
        _S687 = true;

#line 279
    }
    else
    {

#line 279
        _S687 = (abs(cosTheta_0)) >= 1.0;

#line 279
    }

#line 279
    if(_S687)
    {

#line 279
        _S687 = true;

#line 279
    }
    else
    {

#line 279
        _S687 = !isfinite(angularVelocity_0);

#line 279
    }

#line 276
    if(_S687)
    {



        return false;
    }
    float sinThetaSquared_0 = (1.0 - cosTheta_0) * (1.0 + cosTheta_0);
    float _S688 = r_13 * r_13;

#line 284
    float _S689 = a_20 * a_20;

#line 284
    float sigma_9 = _S688 + _S689 * cosTheta_0 * cosTheta_0;
    if(!isfinite(sigma_9))
    {

#line 285
        _S687 = true;

#line 285
    }
    else
    {

#line 285
        _S687 = sigma_9 <= 0.0;

#line 285
    }

#line 285
    if(_S687)
    {

#line 285
        _S687 = true;

#line 285
    }
    else
    {

#line 285
        _S687 = sinThetaSquared_0 <= 0.0;

#line 285
    }

#line 285
    if(_S687)
    {

#line 285
        return false;
    }

#line 286
    float _S690 = 2.0 * M_13;



    float norm_0 = - (1.0 - _S690 * r_13 / sigma_9) + 2.0 * angularVelocity_0 * (-2.0 * M_13 * a_20 * r_13 * sinThetaSquared_0 / sigma_9) + angularVelocity_0 * angularVelocity_0 * (sinThetaSquared_0 * (_S688 + _S689 + _S690 * a_20 * a_20 * r_13 * sinThetaSquared_0 / sigma_9));

    if(!isfinite(norm_0))
    {

#line 292
        _S687 = true;

#line 292
    }
    else
    {

#line 292
        _S687 = norm_0 >= 0.0;

#line 292
    }

#line 292
    if(_S687)
    {

#line 292
        return false;
    }

#line 293
    float timeComponent_0 = 1.0 / sqrt(- norm_0);
    float frameFrequency_0 = timeComponent_0 * (energy_2 - angularVelocity_0 * angularMomentum_2);
    if(!isfinite(frameFrequency_0))
    {

#line 295
        _S687 = true;

#line 295
    }
    else
    {

#line 295
        _S687 = frameFrequency_0 <= 0.0;

#line 295
    }

#line 295
    if(_S687)
    {

#line 295
        return false;
    }

#line 296
    *transfer_0 = float4(observerFrequency_0 / frameFrequency_0, frameFrequency_0, angularVelocity_0, timeComponent_0);

    return true;
}


#line 272
bool TryKerrStationaryFrameFrequencyTransfer_1(float observerFrequency_1, float energy_3, float angularMomentum_3, float M_14, float a_21, float r_14, float cosTheta_1, float angularVelocity_1, float4 thread* transfer_1)
{

    *transfer_1 = float4(0.0) ;

#line 275
    bool _S691;
    if(!isfinite(observerFrequency_1))
    {

#line 276
        _S691 = true;

#line 276
    }
    else
    {

#line 276
        _S691 = observerFrequency_1 <= 0.0;

#line 276
    }

#line 276
    if(_S691)
    {

#line 276
        _S691 = true;

#line 276
    }
    else
    {

#line 276
        _S691 = !isfinite(energy_3);

#line 276
    }
    if(_S691)
    {

#line 277
        _S691 = true;

#line 277
    }
    else
    {

#line 277
        _S691 = energy_3 <= 0.0;

#line 277
    }

#line 277
    if(_S691)
    {

#line 277
        _S691 = true;

#line 277
    }
    else
    {

#line 277
        _S691 = !isfinite(angularMomentum_3);

#line 277
    }

#line 277
    if(_S691)
    {

#line 277
        _S691 = true;

#line 277
    }
    else
    {

#line 277
        _S691 = !isfinite(M_14);

#line 277
    }
    if(_S691)
    {

#line 278
        _S691 = true;

#line 278
    }
    else
    {

#line 278
        _S691 = M_14 <= 0.0;

#line 278
    }

#line 278
    if(_S691)
    {

#line 278
        _S691 = true;

#line 278
    }
    else
    {

#line 278
        _S691 = !isfinite(a_21);

#line 278
    }

#line 278
    if(_S691)
    {

#line 278
        _S691 = true;

#line 278
    }
    else
    {

#line 278
        _S691 = (abs(a_21)) > M_14;

#line 278
    }

#line 278
    if(_S691)
    {

#line 278
        _S691 = true;

#line 278
    }
    else
    {

#line 278
        _S691 = !isfinite(r_14);

#line 278
    }
    if(_S691)
    {

#line 279
        _S691 = true;

#line 279
    }
    else
    {

#line 279
        _S691 = r_14 <= 0.0;

#line 279
    }

#line 279
    if(_S691)
    {

#line 279
        _S691 = true;

#line 279
    }
    else
    {

#line 279
        _S691 = !isfinite(cosTheta_1);

#line 279
    }

#line 279
    if(_S691)
    {

#line 279
        _S691 = true;

#line 279
    }
    else
    {

#line 279
        _S691 = (abs(cosTheta_1)) >= 1.0;

#line 279
    }

#line 279
    if(_S691)
    {

#line 279
        _S691 = true;

#line 279
    }
    else
    {

#line 279
        _S691 = !isfinite(angularVelocity_1);

#line 279
    }

#line 276
    if(_S691)
    {



        return false;
    }
    float sinThetaSquared_1 = (1.0 - cosTheta_1) * (1.0 + cosTheta_1);
    float _S692 = r_14 * r_14;

#line 284
    float _S693 = a_21 * a_21;

#line 284
    float sigma_10 = _S692 + _S693 * cosTheta_1 * cosTheta_1;
    if(!isfinite(sigma_10))
    {

#line 285
        _S691 = true;

#line 285
    }
    else
    {

#line 285
        _S691 = sigma_10 <= 0.0;

#line 285
    }

#line 285
    if(_S691)
    {

#line 285
        _S691 = true;

#line 285
    }
    else
    {

#line 285
        _S691 = sinThetaSquared_1 <= 0.0;

#line 285
    }

#line 285
    if(_S691)
    {

#line 285
        return false;
    }

#line 286
    float _S694 = 2.0 * M_14;



    float norm_1 = - (1.0 - _S694 * r_14 / sigma_10) + 2.0 * angularVelocity_1 * (-2.0 * M_14 * a_21 * r_14 * sinThetaSquared_1 / sigma_10) + angularVelocity_1 * angularVelocity_1 * (sinThetaSquared_1 * (_S692 + _S693 + _S694 * a_21 * a_21 * r_14 * sinThetaSquared_1 / sigma_10));

    if(!isfinite(norm_1))
    {

#line 292
        _S691 = true;

#line 292
    }
    else
    {

#line 292
        _S691 = norm_1 >= 0.0;

#line 292
    }

#line 292
    if(_S691)
    {

#line 292
        return false;
    }

#line 293
    float timeComponent_1 = 1.0 / sqrt(- norm_1);
    float frameFrequency_1 = timeComponent_1 * (energy_3 - angularVelocity_1 * angularMomentum_3);
    if(!isfinite(frameFrequency_1))
    {

#line 295
        _S691 = true;

#line 295
    }
    else
    {

#line 295
        _S691 = frameFrequency_1 <= 0.0;

#line 295
    }

#line 295
    if(_S691)
    {

#line 295
        return false;
    }

#line 296
    *transfer_1 = float4(observerFrequency_1 / frameFrequency_1, frameFrequency_1, angularVelocity_1, timeComponent_1);

    return true;
}




bool KerrZamoFrequencyTransfer_0(float observerFrequency_2, float energy_4, float angularMomentum_4, float M_15, float a_22, float r_15, float cosTheta_2, float4 thread* transfer_2)
{

    *transfer_2 = float4(0.0) ;

#line 307
    bool _S695;
    if(!isfinite(M_15))
    {

#line 308
        _S695 = true;

#line 308
    }
    else
    {

#line 308
        _S695 = M_15 <= 0.0;

#line 308
    }

#line 308
    if(_S695)
    {

#line 308
        _S695 = true;

#line 308
    }
    else
    {

#line 308
        _S695 = !isfinite(a_22);

#line 308
    }

#line 308
    if(_S695)
    {

#line 308
        _S695 = true;

#line 308
    }
    else
    {

#line 308
        _S695 = (abs(a_22)) > M_15;

#line 308
    }

#line 308
    if(_S695)
    {

#line 308
        _S695 = true;

#line 308
    }
    else
    {

#line 308
        _S695 = !isfinite(r_15);

#line 308
    }
    if(_S695)
    {

#line 309
        _S695 = true;

#line 309
    }
    else
    {

#line 309
        _S695 = r_15 <= 0.0;

#line 309
    }

#line 309
    if(_S695)
    {

#line 309
        _S695 = true;

#line 309
    }
    else
    {

#line 309
        _S695 = !isfinite(cosTheta_2);

#line 309
    }

#line 309
    if(_S695)
    {

#line 309
        _S695 = true;

#line 309
    }
    else
    {

#line 309
        _S695 = (abs(cosTheta_2)) >= 1.0;

#line 309
    }

#line 308
    if(_S695)
    {
        return false;
    }
    float sinThetaSquared_2 = (1.0 - cosTheta_2) * (1.0 + cosTheta_2);
    float _S696 = r_15 * r_15;

#line 313
    float _S697 = a_22 * a_22;

#line 313
    float sigma_11 = _S696 + _S697 * cosTheta_2 * cosTheta_2;
    float gtphi_0 = -2.0 * M_15 * a_22 * r_15 * sinThetaSquared_2 / sigma_11;
    float gphiphi_0 = sinThetaSquared_2 * (_S696 + _S697 + 2.0 * M_15 * a_22 * a_22 * r_15 * sinThetaSquared_2 / sigma_11);

    if(!isfinite(gphiphi_0))
    {

#line 317
        _S695 = true;

#line 317
    }
    else
    {

#line 317
        _S695 = gphiphi_0 <= 0.0;

#line 317
    }

#line 317
    if(_S695)
    {

#line 317
        return false;
    }

#line 318
    bool _S698 = TryKerrStationaryFrameFrequencyTransfer_1(observerFrequency_2, energy_4, angularMomentum_4, M_15, a_22, r_15, cosTheta_2, - gtphi_0 / gphiphi_0, transfer_2);

#line 318
    return _S698;
}


#line 14
float ObservedDiskBolometricSource_0(float emittedSource_0, float redshift_0)
{

#line 15
    float g2_0 = redshift_0 * redshift_0;
    return emittedSource_0 * g2_0 * g2_0;
}


#line 59
float CieGaussian_0(float wavelengthNm_0, float centre_0, float leftWidth_0, float rightWidth_0, float amplitude_0)
{

#line 60
    float width_0;
    if(wavelengthNm_0 < centre_0)
    {

#line 61
        width_0 = leftWidth_0;

#line 61
    }
    else
    {

#line 61
        width_0 = rightWidth_0;

#line 61
    }
    float offset_0 = (wavelengthNm_0 - centre_0) * width_0;
    return amplitude_0 * exp(-0.5 * offset_0 * offset_0);
}



float3 Cie1931TwoDegreeFit_0(float wavelengthNm_1)
{

#line 68
    bool _S699;
    if(!isfinite(wavelengthNm_1))
    {

#line 69
        _S699 = true;

#line 69
    }
    else
    {

#line 69
        _S699 = wavelengthNm_1 < 380.0;

#line 69
    }

#line 69
    if(_S699)
    {

#line 69
        _S699 = true;

#line 69
    }
    else
    {

#line 69
        _S699 = wavelengthNm_1 > 780.0;

#line 69
    }

#line 69
    if(_S699)
    {

#line 70
        return float3(0.0) ;
    }

#line 79
    return float3(CieGaussian_0(wavelengthNm_1, 442.0, 0.06239999830722809, 0.03739999979734421, 0.3619999885559082) + CieGaussian_0(wavelengthNm_1, 599.79998779296875, 0.02639999985694885, 0.03229999914765358, 1.0559999942779541) - CieGaussian_0(wavelengthNm_1, 501.100006103515625, 0.04899999871850014, 0.03819999843835831, 0.06499999761581421), CieGaussian_0(wavelengthNm_1, 568.79998779296875, 0.02129999920725822, 0.02470000088214874, 0.82099997997283936) + CieGaussian_0(wavelengthNm_1, 530.9000244140625, 0.06129999831318855, 0.03220000118017197, 0.28600001335144043), CieGaussian_0(wavelengthNm_1, 437.0, 0.08449999988079071, 0.02779999934136868, 1.21700000762939453) + CieGaussian_0(wavelengthNm_1, 459.0, 0.03849999979138374, 0.07249999791383743, 0.6809999942779541));
}




float3 XyzD65ToLinearSrgb_0(float3 xyz_0)
{

#line 86
    thread float3 rgb_0;
    float _S700 = xyz_0.x;

#line 87
    float _S701 = xyz_0.y;
    float _S702 = xyz_0.z;

#line 87
    rgb_0.x = 3.24096989631652832 * _S700 - 1.53738319873809814 * _S701 - 0.49861076474189758 * _S702;

    rgb_0.y = -0.96924364566802979 * _S700 + 1.87596750259399414 * _S701 + 0.04155505821108818 * _S702;

    rgb_0.z = 0.05563008040189743 * _S700 - 0.20397695899009705 * _S701 + 1.05697154998779297 * _S702;

    return rgb_0;
}



float3 BlackbodyColor_0(float T_0)
{

#line 99
    if(T_0 <= 0.0)
    {

#line 99
        return float3(0.0) ;
    }

#line 100
    float3 _S703 = float3(0.0) ;

#line 100
    int sample_1 = int(0);

#line 100
    float3 xyz_1 = _S703;

#line 105
    for(;;)
    {

#line 105
        if(sample_1 < int(32))
        {
        }
        else
        {

#line 105
            break;
        }

#line 106
        float wavelengthNm_2 = 380.0 + (float(sample_1) + 0.5) * 12.5;



        float3 xyz_2 = xyz_1 + Cie1931TwoDegreeFit_0(wavelengthNm_2) * float3((1.0 / (pow(wavelengthNm_2, 5.0) * (exp(1.4387769e+07 / (wavelengthNm_2 * T_0)) - 1.0)))) ;

#line 105
        sample_1 = sample_1 + int(1);

#line 105
        xyz_1 = xyz_2;

#line 105
    }

#line 114
    float3 rgb_1 = max(XyzD65ToLinearSrgb_0(xyz_1), _S703);
    return rgb_1 / float3(max(max(rgb_1.x, rgb_1.y), max(rgb_1.z, 9.99999968265522539e-21))) ;
}


#line 22
float GreyLayerAbsorbedFraction_0(float deltaTau_0)
{

#line 22
    bool _S704;
    if(!isfinite(deltaTau_0))
    {

#line 23
        _S704 = true;

#line 23
    }
    else
    {

#line 23
        _S704 = deltaTau_0 < 0.0;

#line 23
    }

#line 23
    if(_S704)
    {

#line 23
        return 0.0;
    }

#line 24
    if(deltaTau_0 <= 0.00999999977648258)
    {

#line 25
        return deltaTau_0 * (1.0 - deltaTau_0 * (0.5 - deltaTau_0 / 6.0));
    }
    return 1.0 - exp(- deltaTau_0);
}


#line 319 "/home/astra/.project/Sirius/out/trace-jacobi-shared-segment-rhs-screen/inputs/candidate/trace.slang"
void AccumulateVolumeSegment_0(const Vec4Cart_0 thread* start_7, const Vec4Cart_0 thread* finish_6, const Vec4Cart_0 thread* startVelocity_6, const Vec4Cart_0 thread* finishVelocity_6, float affineLength_7, int sampleCount_0, bool volumeEnabled_0, bool turbulenceEnabled_0, bool dopplerBeaming_0, float innerRadius_3, float outerRadius_1, float innerTemperature_1, int temperatureModel_1, float M_16, float a_23, float aStar_3, float hOverR_3, float hPower_3, float tauMidplane_0, float observerFrequency_3, float thread* opticalDepth_0, float3 thread* emission_0)
{

#line 328
    bool _S705;
    if(affineLength_7 <= 0.0)
    {

#line 329
        _S705 = true;

#line 329
    }
    else
    {

#line 329
        _S705 = (*opticalDepth_0) >= 10.0;

#line 329
    }

#line 329
    if(_S705)
    {

#line 329
        return;
    }

#line 330
    float _S706 = float(sampleCount_0);

#line 330
    float _S707 = affineLength_7 / _S706;

#line 330
    int sample_2 = int(0);
    for(;;)
    {

#line 331
        if(sample_2 < sampleCount_0)
        {
        }
        else
        {

#line 331
            break;
        }

#line 332
        float fraction_5 = (float(sample_2) + 0.5) / _S706;

#line 332
        Vec4Cart_0 _S708 = AcceptedSegmentPosition_0(start_7, startVelocity_6, finish_6, finishVelocity_6, affineLength_7, fraction_5);

#line 332
        Vec4Cart_0 _S709 = AcceptedSegmentTangent_0(start_7, startVelocity_6, finish_6, finishVelocity_6, affineLength_7, fraction_5);

#line 339
        float diskR_0 = ComputeKerrRadiusCart_0(_S708.x_2, _S708.y_0, _S708.z_0, a_23);


        float3 _S710 = float3(0.0) ;

        if(volumeEnabled_0)
        {

#line 344
            _S705 = IsInsideVolumeDisk_0(diskR_0, _S708.z_0, innerRadius_3, outerRadius_1, hOverR_3, hPower_3);

#line 344
        }
        else
        {

#line 344
            _S705 = false;

#line 344
        }

#line 344
        float diskDeltaTau_0;

#line 344
        float3 diskSource_0;

#line 344
        if(_S705)
        {


            float diskOpacity_0 = TruncatedGaussianOpacityDensity_0(tauMidplane_0, diskR_0, innerRadius_3, VolumeScaleHeight_0(diskR_0, innerRadius_3, hOverR_3, hPower_3), _S708.z_0);

#line 348
            float diskOpacity_1;

            if(turbulenceEnabled_0)
            {


                float diskPhi_0 = atan2(_S708.y_0, _S708.x_2);

#line 354
                diskOpacity_1 = diskOpacity_0 * TurbulenceDensityFactor_0(float3(diskR_0 * cos(diskPhi_0), diskR_0 * sin(diskPhi_0), _S708.z_0));

#line 350
            }
            else
            {

#line 350
                diskOpacity_1 = diskOpacity_0;

#line 350
            }

#line 360
            float temperature_0 = VolumeDiskTemperature_0(diskR_0, _S708.z_0, innerRadius_3, innerTemperature_1, temperatureModel_1, M_16, aStar_3, hOverR_3, hPower_3);

#line 360
            bool _S711;


            if(isfinite(diskOpacity_1))
            {

#line 363
                _S711 = diskOpacity_1 > 0.0;

#line 363
            }
            else
            {

#line 363
                _S711 = false;

#line 363
            }

#line 363
            if(_S711)
            {

#line 363
                thread Vec4Cart_0 _S712 = _S708;


                thread array<array<float, int(4)>, int(4)> localMetric_0;
                thread array<array<float, int(4)>, int(4)> localInverse_0;

#line 367
                GetKerrSchildFamilyMetric_0(&_S712, M_16, a_23, 0.0, 0.0, &localMetric_0, &localInverse_0);

#line 367
                thread Vec4Cart_0 _S713 = _S708;

#line 367
                thread Vec4Cart_0 _S714 = _S709;

#line 367
                thread array<array<float, int(4)>, int(4)> _S715 = localMetric_0;


                thread float energy_5;
                thread float angularMomentum_5;

#line 371
                bool _S716 = KerrPhotonKillingQuantitiesCart_0(&_S713, &_S714, &_S715, &energy_5, &angularMomentum_5);
                if(!_S716)
                {
                    sample_2 = sample_2 + int(1);

#line 331
                    continue;
                }

#line 376
                float sqrtM_0 = sqrt(M_16);

                float cosTheta_3 = _S708.z_0 / diskR_0;

                thread float4 zamoTransfer_0;

#line 379
                thread float4 emitterTransfer_0;

                bool _S717 = TryKerrStationaryFrameFrequencyTransfer_0(observerFrequency_3, energy_5, angularMomentum_5, M_16, a_23, diskR_0, cosTheta_3, sqrtM_0 / (pow(diskR_0, 1.5) + a_23 * sqrtM_0), &emitterTransfer_0);

#line 381
                bool _S718;

                if(!_S717)
                {

#line 383
                    _S718 = true;

#line 383
                }
                else
                {

#line 384
                    bool _S719 = KerrZamoFrequencyTransfer_0(observerFrequency_3, energy_5, angularMomentum_5, M_16, a_23, diskR_0, cosTheta_3, &zamoTransfer_0);

#line 384
                    _S718 = !_S719;

#line 383
                }

#line 381
                if(_S718)
                {

#line 387
                    sample_2 = sample_2 + int(1);

#line 331
                    continue;
                }

#line 393
                float emittedFrequency_0 = emitterTransfer_0.y;
                float _S720 = diskOpacity_1 * emittedFrequency_0 * _S707;

#line 394
                bool _S721;


                if(isfinite(temperature_0))
                {

#line 397
                    _S721 = temperature_0 > 0.0;

#line 397
                }
                else
                {

#line 397
                    _S721 = false;

#line 397
                }

#line 397
                if(_S721)
                {

#line 398
                    float emittedSource_1 = pow(temperature_0 / innerTemperature_1, 4.0);
                    if(dopplerBeaming_0)
                    {

#line 399
                        diskDeltaTau_0 = emittedFrequency_0;

#line 399
                    }
                    else
                    {

#line 399
                        diskDeltaTau_0 = zamoTransfer_0.y;

#line 399
                    }
                    float gFactor_0 = observerFrequency_3 / diskDeltaTau_0;

#line 400
                    diskSource_0 = BlackbodyColor_0(temperature_0 * gFactor_0) * float3(ObservedDiskBolometricSource_0(emittedSource_1, gFactor_0)) ;

#line 397
                }
                else
                {

#line 397
                    diskSource_0 = _S710;

#line 397
                }

#line 397
                diskDeltaTau_0 = _S720;

#line 363
            }
            else
            {

#line 363
                diskDeltaTau_0 = 0.0;

#line 363
                diskSource_0 = _S710;

#line 363
            }

#line 344
        }
        else
        {

#line 344
            diskDeltaTau_0 = 0.0;

#line 344
            diskSource_0 = _S710;

#line 344
        }

#line 408
        if(diskDeltaTau_0 <= 0.0)
        {

#line 408
            sample_2 = sample_2 + int(1);

#line 331
            continue;
        }

#line 409
        float _S722 = min(diskDeltaTau_0, 10.0 - *opticalDepth_0);


        *emission_0 = *emission_0 + float3(exp(- *opticalDepth_0))  * diskSource_0 * float3(GreyLayerAbsorbedFraction_0(_S722)) ;
        float _S723 = *opticalDepth_0 + _S722;

#line 413
        *opticalDepth_0 = _S723;
        if(_S723 >= 10.0)
        {

#line 414
            break;
        }

#line 331
        sample_2 = sample_2 + int(1);

#line 331
    }

#line 416
    return;
}


#line 310
float EvaluateCubic_0(float coefficientA_0, float coefficientB_0, float coefficientC_0, float coefficientD_0, float fraction_6)
{
    return ((coefficientA_0 * fraction_6 + coefficientB_0) * fraction_6 + coefficientC_0) * fraction_6 + coefficientD_0;
}


#line 326 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_disk.slang"
float2 ComputeKerrDiskTransferCart_0(float r_16, float phi_0, const Vec4Cart_0 thread* photonVelocity_0, float observerFrequency_4, float M_17, float a_24)
{

#line 327
    bool _S724;
    if(M_17 <= 0.0)
    {

#line 328
        _S724 = true;

#line 328
    }
    else
    {

#line 328
        _S724 = r_16 <= 0.0;

#line 328
    }

#line 328
    if(_S724)
    {

#line 328
        _S724 = true;

#line 328
    }
    else
    {

#line 328
        _S724 = observerFrequency_4 <= 0.0;

#line 328
    }

#line 328
    if(_S724)
    {

#line 328
        _S724 = true;

#line 328
    }
    else
    {

#line 328
        _S724 = (abs(a_24)) > M_17;

#line 328
    }

#line 328
    if(_S724)
    {

#line 329
        return float2(0.0) ;
    }

    float rho_3 = sqrt(r_16 * r_16 + a_24 * a_24);

    Vec4Cart_0 position_1 = MakeVec4Cart_0(0.0, rho_3 * cos(phi_0), rho_3 * sin(phi_0), 0.0);

#line 334
    thread Vec4Cart_0 _S725 = position_1;
    thread array<array<float, int(4)>, int(4)> metric_9;
    thread array<array<float, int(4)>, int(4)> metricInverse_0;

#line 336
    GetKerrSchildFamilyMetric_0(&_S725, M_17, a_24, 0.0, 0.0, &metric_9, &metricInverse_0);

#line 336
    thread Vec4Cart_0 _S726 = position_1;

#line 336
    thread array<array<float, int(4)>, int(4)> _S727 = metric_9;


    thread float energy_6;
    thread float angularMomentum_6;

#line 340
    bool _S728 = KerrPhotonKillingQuantitiesCart_0(&_S726, photonVelocity_0, &_S727, &energy_6, &angularMomentum_6);
    if(!_S728)
    {
        return float2(0.0) ;
    }
    float sqrtM_1 = sqrt(M_17);


    thread float4 zamoTransfer_1;

#line 347
    thread float4 emitterTransfer_1;

    bool _S729 = TryKerrStationaryFrameFrequencyTransfer_0(observerFrequency_4, energy_6, angularMomentum_6, M_17, a_24, r_16, 0.0, sqrtM_1 / (pow(r_16, 1.5) + a_24 * sqrtM_1), &emitterTransfer_1);

    if(!_S729)
    {

#line 351
        _S724 = true;

#line 351
    }
    else
    {

#line 352
        bool _S730 = KerrZamoFrequencyTransfer_0(observerFrequency_4, energy_6, angularMomentum_6, M_17, a_24, r_16, 0.0, &zamoTransfer_1);

#line 352
        _S724 = !_S730;

#line 351
    }

#line 349
    if(_S724)
    {



        return float2(0.0) ;
    }
    return float2(emitterTransfer_1.x, zamoTransfer_1.x);
}


#line 27 "/home/astra/.project/Sirius/out/trace-jacobi-shared-segment-rhs-screen/inputs/candidate/gr_deviation.slang"
void GetRiemannTensorCart_0(const Vec4Cart_0 thread* x_8, float M_18, float a_25, float Q_11, float Lambda_11, array<array<array<array<float, int(4)>, int(4)>, int(4)>, int(4)> thread* R_8)
{
    int rho_4;

#line 29
    int nu_12;

#line 29
    int sig_0;

#line 29
    int distance_0;

#line 29
    int mu_9 = int(0);

    for(;;)
    {

#line 31
        if(mu_9 < int(4))
        {
        }
        else
        {

#line 31
            break;
        }

#line 31
        nu_12 = int(0);
        for(;;)
        {

#line 32
            if(nu_12 < int(4))
            {
            }
            else
            {

#line 32
                break;
            }

#line 32
            rho_4 = int(0);
            for(;;)
            {

#line 33
                if(rho_4 < int(4))
                {
                }
                else
                {

#line 33
                    break;
                }

#line 33
                sig_0 = int(0);
                for(;;)
                {

#line 34
                    if(sig_0 < int(4))
                    {
                    }
                    else
                    {

#line 34
                        break;
                    }

#line 35
                    (*R_8)[mu_9][nu_12][rho_4][sig_0] = 0.0;

#line 34
                    sig_0 = sig_0 + int(1);

#line 34
                }

#line 33
                rho_4 = rho_4 + int(1);

#line 33
            }

#line 32
            nu_12 = nu_12 + int(1);

#line 32
        }

#line 31
        mu_9 = mu_9 + int(1);

#line 31
    }

#line 37
    thread array<array<array<float, int(4)>, int(4)>, int(4)> Gamma_6;

#line 37
    GetKerrSchildFamilyChristoffel_0(x_8, M_18, a_25, Q_11, Lambda_11, &Gamma_6);

#line 51
    float _S731 = 0.03999999910593033 * max(1.0, ComputeKerrRadiusCart_0(x_8->x_2, x_8->y_0, x_8->z_0, a_25));



    thread array<array<array<array<float, int(4)>, int(4)>, int(4)>, int(4)> dGamma_0;

#line 55
    mu_9 = int(0);
    for(;;)
    {

#line 56
        if(mu_9 < int(4))
        {
        }
        else
        {

#line 56
            break;
        }

#line 56
        nu_12 = int(0);
        for(;;)
        {

#line 57
            if(nu_12 < int(4))
            {
            }
            else
            {

#line 57
                break;
            }

#line 57
            sig_0 = int(0);
            for(;;)
            {

#line 58
                if(sig_0 < int(4))
                {
                }
                else
                {

#line 58
                    break;
                }

#line 59
                dGamma_0[int(0)][mu_9][nu_12][sig_0] = 0.0;

#line 58
                sig_0 = sig_0 + int(1);

#line 58
            }

#line 57
            nu_12 = nu_12 + int(1);

#line 57
        }

#line 56
        mu_9 = mu_9 + int(1);

#line 56
    }

#line 56
    rho_4 = int(1);

#line 61
    for(;;)
    {

#line 61
        if(rho_4 < int(4))
        {
        }
        else
        {

#line 61
            break;
        }
        thread array<array<array<float, int(4)>, int(4)>, int(4)> Gamma_plus_0;
        thread array<array<array<float, int(4)>, int(4)>, int(4)> Gamma_minus_0;

#line 64
        distance_0 = int(1);



        for(;;)
        {

#line 68
            if(distance_0 <= int(3))
            {
            }
            else
            {

#line 68
                break;
            }
            thread Vec4Cart_0 xp_0 = *x_8;
            thread Vec4Cart_0 xm_0 = *x_8;
            if(rho_4 == int(1))
            {

#line 73
                if(distance_0 == int(1))
                {

#line 74
                    (&xp_0)->x_2 = (&xp_0)->x_2 + _S731;
                    (&xm_0)->x_2 = (&xm_0)->x_2 - _S731;

#line 73
                }
                else
                {
                    if(distance_0 == int(2))
                    {

#line 77
                        float _S732 = 2.0 * _S731;

#line 77
                        (&xp_0)->x_2 = (&xp_0)->x_2 + _S732;
                        (&xm_0)->x_2 = (&xm_0)->x_2 - _S732;

#line 76
                    }
                    else
                    {


                        float _S733 = 3.0 * _S731;

#line 81
                        (&xp_0)->x_2 = (&xp_0)->x_2 + _S733;
                        (&xm_0)->x_2 = (&xm_0)->x_2 - _S733;

#line 76
                    }

#line 73
                }

#line 72
            }
            else
            {

#line 85
                if(rho_4 == int(2))
                {

#line 86
                    if(distance_0 == int(1))
                    {

#line 87
                        (&xp_0)->y_0 = (&xp_0)->y_0 + _S731;
                        (&xm_0)->y_0 = (&xm_0)->y_0 - _S731;

#line 86
                    }
                    else
                    {
                        if(distance_0 == int(2))
                        {

#line 90
                            float _S734 = 2.0 * _S731;

#line 90
                            (&xp_0)->y_0 = (&xp_0)->y_0 + _S734;
                            (&xm_0)->y_0 = (&xm_0)->y_0 - _S734;

#line 89
                        }
                        else
                        {


                            float _S735 = 3.0 * _S731;

#line 94
                            (&xp_0)->y_0 = (&xp_0)->y_0 + _S735;
                            (&xm_0)->y_0 = (&xm_0)->y_0 - _S735;

#line 89
                        }

#line 86
                    }

#line 85
                }
                else
                {

#line 99
                    if(distance_0 == int(1))
                    {

#line 100
                        (&xp_0)->z_0 = (&xp_0)->z_0 + _S731;
                        (&xm_0)->z_0 = (&xm_0)->z_0 - _S731;

#line 99
                    }
                    else
                    {
                        if(distance_0 == int(2))
                        {

#line 103
                            float _S736 = 2.0 * _S731;

#line 103
                            (&xp_0)->z_0 = (&xp_0)->z_0 + _S736;
                            (&xm_0)->z_0 = (&xm_0)->z_0 - _S736;

#line 102
                        }
                        else
                        {


                            float _S737 = 3.0 * _S731;

#line 107
                            (&xp_0)->z_0 = (&xp_0)->z_0 + _S737;
                            (&xm_0)->z_0 = (&xm_0)->z_0 - _S737;

#line 102
                        }

#line 99
                    }

#line 85
                }

#line 72
            }

#line 72
            thread Vec4Cart_0 _S738 = xp_0;

#line 72
            GetKerrSchildFamilyChristoffel_0(&_S738, M_18, a_25, Q_11, Lambda_11, &Gamma_plus_0);

#line 72
            thread Vec4Cart_0 _S739 = xm_0;

#line 72
            GetKerrSchildFamilyChristoffel_0(&_S739, M_18, a_25, Q_11, Lambda_11, &Gamma_minus_0);

#line 116
            if(distance_0 == int(1))
            {

#line 116
                mu_9 = int(0);
                for(;;)
                {

#line 117
                    if(mu_9 < int(4))
                    {
                    }
                    else
                    {

#line 117
                        break;
                    }

#line 117
                    nu_12 = int(0);
                    for(;;)
                    {

#line 118
                        if(nu_12 < int(4))
                        {
                        }
                        else
                        {

#line 118
                            break;
                        }

#line 118
                        sig_0 = int(0);
                        for(;;)
                        {

#line 119
                            if(sig_0 < int(4))
                            {
                            }
                            else
                            {

#line 119
                                break;
                            }

#line 120
                            dGamma_0[rho_4][mu_9][nu_12][sig_0] = 45.0 * (Gamma_plus_0[mu_9][nu_12][sig_0] - Gamma_minus_0[mu_9][nu_12][sig_0]);

#line 119
                            sig_0 = sig_0 + int(1);

#line 119
                        }

#line 118
                        nu_12 = nu_12 + int(1);

#line 118
                    }

#line 117
                    mu_9 = mu_9 + int(1);

#line 117
                }

#line 116
            }
            else
            {

#line 126
                if(distance_0 == int(2))
                {

#line 126
                    mu_9 = int(0);
                    for(;;)
                    {

#line 127
                        if(mu_9 < int(4))
                        {
                        }
                        else
                        {

#line 127
                            break;
                        }

#line 127
                        nu_12 = int(0);
                        for(;;)
                        {

#line 128
                            if(nu_12 < int(4))
                            {
                            }
                            else
                            {

#line 128
                                break;
                            }

#line 128
                            sig_0 = int(0);
                            for(;;)
                            {

#line 129
                                if(sig_0 < int(4))
                                {
                                }
                                else
                                {

#line 129
                                    break;
                                }

#line 136
                                dGamma_0[rho_4][mu_9][nu_12][sig_0] = dGamma_0[rho_4][mu_9][nu_12][sig_0] - 9.0 * (Gamma_plus_0[mu_9][nu_12][sig_0] - Gamma_minus_0[mu_9][nu_12][sig_0]);

#line 129
                                sig_0 = sig_0 + int(1);

#line 129
                            }

#line 128
                            nu_12 = nu_12 + int(1);

#line 128
                        }

#line 127
                        mu_9 = mu_9 + int(1);

#line 127
                    }

#line 126
                }
                else
                {

#line 126
                    mu_9 = int(0);

#line 141
                    for(;;)
                    {

#line 141
                        if(mu_9 < int(4))
                        {
                        }
                        else
                        {

#line 141
                            break;
                        }

#line 141
                        nu_12 = int(0);
                        for(;;)
                        {

#line 142
                            if(nu_12 < int(4))
                            {
                            }
                            else
                            {

#line 142
                                break;
                            }

#line 142
                            sig_0 = int(0);
                            for(;;)
                            {

#line 143
                                if(sig_0 < int(4))
                                {
                                }
                                else
                                {

#line 143
                                    break;
                                }

#line 144
                                dGamma_0[rho_4][mu_9][nu_12][sig_0] = (dGamma_0[rho_4][mu_9][nu_12][sig_0] + (Gamma_plus_0[mu_9][nu_12][sig_0] - Gamma_minus_0[mu_9][nu_12][sig_0])) / (60.0 * _S731);

#line 143
                                sig_0 = sig_0 + int(1);

#line 143
                            }

#line 142
                            nu_12 = nu_12 + int(1);

#line 142
                        }

#line 141
                        mu_9 = mu_9 + int(1);

#line 141
                    }

#line 126
                }

#line 116
            }

#line 68
            distance_0 = distance_0 + int(1);

#line 68
        }

#line 61
        rho_4 = rho_4 + int(1);

#line 61
    }

#line 61
    mu_9 = int(0);

#line 156
    for(;;)
    {

#line 156
        if(mu_9 < int(4))
        {
        }
        else
        {

#line 156
            break;
        }

#line 156
        distance_0 = int(0);
        for(;;)
        {

#line 157
            if(distance_0 < int(4))
            {
            }
            else
            {

#line 157
                break;
            }

#line 157
            rho_4 = int(0);
            for(;;)
            {

#line 158
                if(rho_4 < int(4))
                {
                }
                else
                {

#line 158
                    break;
                }

#line 158
                nu_12 = int(0);
                for(;;)
                {

#line 159
                    if(nu_12 < int(4))
                    {
                    }
                    else
                    {

#line 159
                        break;
                    }

#line 160
                    float _S740 = dGamma_0[rho_4][mu_9][distance_0][nu_12] - dGamma_0[nu_12][mu_9][distance_0][rho_4];

#line 160
                    sig_0 = int(0);

#line 160
                    float val_0 = _S740;
                    for(;;)
                    {

#line 161
                        if(sig_0 < int(4))
                        {
                        }
                        else
                        {

#line 161
                            break;
                        }

#line 162
                        float val_1 = val_0 + (Gamma_6[mu_9][sig_0][rho_4] * Gamma_6[sig_0][distance_0][nu_12] - Gamma_6[mu_9][sig_0][nu_12] * Gamma_6[sig_0][distance_0][rho_4]);

#line 161
                        sig_0 = sig_0 + int(1);

#line 161
                        val_0 = val_1;

#line 161
                    }


                    (*R_8)[mu_9][distance_0][rho_4][nu_12] = val_0;

#line 159
                    nu_12 = nu_12 + int(1);

#line 159
                }

#line 158
                rho_4 = rho_4 + int(1);

#line 158
            }

#line 157
            distance_0 = distance_0 + int(1);

#line 157
        }

#line 156
        mu_9 = mu_9 + int(1);

#line 156
    }

#line 166
    return;
}


#line 190
void GeodesicDeviationRhs_0(const array<array<array<float, int(4)>, int(4)>, int(4)> thread* Gamma_7, const array<array<array<array<float, int(4)>, int(4)>, int(4)>, int(4)> thread* R_9, const Vec4Cart_0 thread* k_3, float4 xi_0, float4 V_1, float4 thread* dXi_0, float4 thread* dV_0)
{
    int beta_0;
    array<float, int(4)> _S741 = { k_3->t_0, k_3->x_2, k_3->y_0, k_3->z_0 };
    array<float, int(4)> _S742 = { xi_0.x, xi_0.y, xi_0.z, xi_0.w };
    array<float, int(4)> _S743 = { V_1.x, V_1.y, V_1.z, V_1.w };
    thread array<float, int(4)> d_xi_0;

#line 196
    d_xi_0[int(0)] = 0.0;

#line 196
    d_xi_0[int(1)] = 0.0;

#line 196
    d_xi_0[int(2)] = 0.0;

#line 196
    d_xi_0[int(3)] = 0.0;
    thread array<float, int(4)> d_v_0;

#line 197
    d_v_0[int(0)] = 0.0;

#line 197
    d_v_0[int(1)] = 0.0;

#line 197
    d_v_0[int(2)] = 0.0;

#line 197
    d_v_0[int(3)] = 0.0;

#line 197
    int mu_10 = int(0);

    for(;;)
    {

#line 199
        if(mu_10 < int(4))
        {
        }
        else
        {

#line 199
            break;
        }

#line 200
        d_xi_0[mu_10] = _S743[mu_10];

#line 200
        int alpha_0 = int(0);
        for(;;)
        {

#line 201
            if(alpha_0 < int(4))
            {
            }
            else
            {

#line 201
                break;
            }

#line 201
            beta_0 = int(0);
            for(;;)
            {

#line 202
                if(beta_0 < int(4))
                {
                }
                else
                {

#line 202
                    break;
                }

#line 203
                float _S744 = (*Gamma_7)[mu_10][alpha_0][beta_0] * _S741[alpha_0];

#line 203
                d_xi_0[mu_10] = d_xi_0[mu_10] - _S744 * _S742[beta_0];
                d_v_0[mu_10] = d_v_0[mu_10] - _S744 * _S743[beta_0];

#line 202
                beta_0 = beta_0 + int(1);

#line 202
            }

#line 201
            alpha_0 = alpha_0 + int(1);

#line 201
        }

#line 201
        beta_0 = int(0);

#line 207
        for(;;)
        {

#line 207
            if(beta_0 < int(4))
            {
            }
            else
            {

#line 207
                break;
            }

#line 207
            int rho_5 = int(0);
            for(;;)
            {

#line 208
                if(rho_5 < int(4))
                {
                }
                else
                {

#line 208
                    break;
                }

#line 208
                int sigma_12 = int(0);
                for(;;)
                {

#line 209
                    if(sigma_12 < int(4))
                    {
                    }
                    else
                    {

#line 209
                        break;
                    }

#line 210
                    d_v_0[mu_10] = d_v_0[mu_10] - (*R_9)[mu_10][beta_0][rho_5][sigma_12] * _S741[beta_0] * _S742[rho_5] * _S741[sigma_12];

#line 209
                    sigma_12 = sigma_12 + int(1);

#line 209
                }

#line 208
                rho_5 = rho_5 + int(1);

#line 208
            }

#line 207
            beta_0 = beta_0 + int(1);

#line 207
        }

#line 199
        mu_10 = mu_10 + int(1);

#line 199
    }

#line 212
    *dXi_0 = float4(d_xi_0[int(0)], d_xi_0[int(1)], d_xi_0[int(2)], d_xi_0[int(3)]);
    *dV_0 = float4(d_v_0[int(0)], d_v_0[int(1)], d_v_0[int(2)], d_v_0[int(3)]);
    return;
}




void AdvanceJacobiPairAcceptedSegment_0(const Vec4Cart_0 thread* startPosition_0, const Vec4Cart_0 thread* startTangent_0, const Vec4Cart_0 thread* endPosition_0, const Vec4Cart_0 thread* endTangent_0, float dLambda_0, float M_19, float a_26, float Q_12, float Lambda_12, float4 thread* xiA_0, float4 thread* VA_0, float4 thread* xiB_0, float4 thread* VB_0)
{



    thread array<array<array<array<float, int(4)>, int(4)>, int(4)>, int(4)> R_10;
    thread array<array<array<float, int(4)>, int(4)>, int(4)> Gamma_8;

#line 226
    float4 k1XiA_0;

#line 226
    float4 k1VA_0;

#line 226
    float4 k2XiA_0;

#line 226
    float4 k2VA_0;

#line 226
    float4 k1XiB_0;

#line 226
    float4 k1VB_0;

#line 226
    float4 k2XiB_0;

#line 226
    float4 k2VB_0;

#line 226
    float4 k3XiA_0;

#line 226
    float4 k3VA_0;

#line 226
    float4 k3XiB_0;

#line 226
    float4 k3VB_0;

#line 226
    float4 k4XiA_0;

#line 226
    float4 k4VA_0;

#line 226
    float4 k4XiB_0;

#line 226
    float4 k4VB_0;

#line 226
    Vec4Cart_0 phaseTangent_0;

#line 226
    int step_6 = int(0);

#line 231
    for(;;)
    {

#line 231
        if(step_6 < int(8))
        {
        }
        else
        {

#line 231
            break;
        }

#line 232
        bool _S745 = step_6 == int(0);

#line 232
        bool _S746;

#line 232
        if(_S745)
        {

#line 232
            _S746 = true;

#line 232
        }
        else
        {

#line 232
            _S746 = step_6 == int(2);

#line 232
        }

#line 232
        bool _S747;

#line 232
        if(_S746)
        {

#line 232
            _S747 = true;

#line 232
        }
        else
        {

#line 232
            _S747 = step_6 == int(6);

#line 232
        }

#line 232
        Vec4Cart_0 phaseTangent_1;

#line 232
        if(_S747)
        {

#line 232
            Vec4Cart_0 phasePosition_0;

            if(_S745)
            {

#line 234
                phasePosition_0 = *startPosition_0;

#line 234
                phaseTangent_1 = *startTangent_0;

#line 234
            }
            else
            {
                if(step_6 == int(2))
                {

#line 237
                    Vec4Cart_0 _S748 = AddVec4Cart_0(startPosition_0, endPosition_0);

#line 237
                    thread Vec4Cart_0 _S749 = _S748;

#line 237
                    Vec4Cart_0 _S750 = ScaleVec4Cart_0(&_S749, 0.5);

#line 237
                    Vec4Cart_0 _S751 = ScaleVec4Cart_0(endTangent_0, -1.0);

#line 237
                    thread Vec4Cart_0 _S752 = _S751;

#line 237
                    Vec4Cart_0 _S753 = AddVec4Cart_0(startTangent_0, &_S752);

#line 242
                    float _S754 = dLambda_0 / 8.0;

#line 242
                    thread Vec4Cart_0 _S755 = _S753;

#line 242
                    Vec4Cart_0 _S756 = ScaleVec4Cart_0(&_S755, _S754);

#line 242
                    thread Vec4Cart_0 _S757 = _S750;

#line 242
                    thread Vec4Cart_0 _S758 = _S756;

#line 242
                    Vec4Cart_0 _S759 = AddVec4Cart_0(&_S757, &_S758);

#line 242
                    Vec4Cart_0 _S760 = ScaleVec4Cart_0(startPosition_0, -1.0);

#line 242
                    thread Vec4Cart_0 _S761 = _S760;

#line 242
                    Vec4Cart_0 _S762 = AddVec4Cart_0(endPosition_0, &_S761);



                    float _S763 = 1.5 / dLambda_0;

#line 246
                    thread Vec4Cart_0 _S764 = _S762;

#line 246
                    Vec4Cart_0 _S765 = ScaleVec4Cart_0(&_S764, _S763);

#line 246
                    Vec4Cart_0 _S766 = AddVec4Cart_0(startTangent_0, endTangent_0);

#line 246
                    thread Vec4Cart_0 _S767 = _S766;

#line 246
                    Vec4Cart_0 _S768 = ScaleVec4Cart_0(&_S767, -0.25);

#line 246
                    thread Vec4Cart_0 _S769 = _S765;

#line 246
                    thread Vec4Cart_0 _S770 = _S768;

#line 246
                    Vec4Cart_0 _S771 = AddVec4Cart_0(&_S769, &_S770);

#line 246
                    phasePosition_0 = _S759;

#line 246
                    phaseTangent_1 = _S771;

#line 237
                }
                else
                {

#line 237
                    phasePosition_0 = *endPosition_0;

#line 237
                    phaseTangent_1 = *endTangent_0;

#line 237
                }

#line 234
            }

#line 234
            thread Vec4Cart_0 _S772 = phasePosition_0;

#line 234
            GetRiemannTensorCart_0(&_S772, M_19, a_26, Q_12, Lambda_12, &R_10);

#line 234
            thread Vec4Cart_0 _S773 = phasePosition_0;

#line 234
            GetKerrSchildFamilyChristoffel_0(&_S773, M_19, a_26, Q_12, Lambda_12, &Gamma_8);

#line 232
        }
        else
        {

#line 232
            phaseTangent_1 = phaseTangent_0;

#line 232
        }

#line 232
        float4 stepXi_0;

#line 232
        float4 stepV_0;

#line 259
        if(_S745)
        {

#line 259
            stepXi_0 = *xiA_0;

#line 259
            stepV_0 = *VA_0;

#line 259
        }
        else
        {
            if(step_6 == int(1))
            {

#line 262
                stepXi_0 = *xiB_0;

#line 262
                stepV_0 = *VB_0;

#line 262
            }
            else
            {
                if(step_6 == int(2))
                {

#line 265
                    float4 _S774 = float4((0.5 * dLambda_0)) ;

                    float4 _S775 = *VA_0 + k1VA_0 * _S774;

#line 267
                    stepXi_0 = *xiA_0 + k1XiA_0 * _S774;

#line 267
                    stepV_0 = _S775;

#line 265
                }
                else
                {
                    if(step_6 == int(3))
                    {

#line 268
                        float4 _S776 = float4((0.5 * dLambda_0)) ;

                        float4 _S777 = *VA_0 + k2VA_0 * _S776;

#line 270
                        stepXi_0 = *xiA_0 + k2XiA_0 * _S776;

#line 270
                        stepV_0 = _S777;

#line 268
                    }
                    else
                    {
                        if(step_6 == int(4))
                        {

#line 271
                            float4 _S778 = float4((0.5 * dLambda_0)) ;

                            float4 _S779 = *VB_0 + k1VB_0 * _S778;

#line 273
                            stepXi_0 = *xiB_0 + k1XiB_0 * _S778;

#line 273
                            stepV_0 = _S779;

#line 271
                        }
                        else
                        {
                            if(step_6 == int(5))
                            {

#line 274
                                float4 _S780 = float4((0.5 * dLambda_0)) ;

                                float4 _S781 = *VB_0 + k2VB_0 * _S780;

#line 276
                                stepXi_0 = *xiB_0 + k2XiB_0 * _S780;

#line 276
                                stepV_0 = _S781;

#line 274
                            }
                            else
                            {
                                if(step_6 == int(6))
                                {

#line 277
                                    float4 _S782 = float4(dLambda_0) ;

                                    float4 _S783 = *VA_0 + k3VA_0 * _S782;

#line 279
                                    stepXi_0 = *xiA_0 + k3XiA_0 * _S782;

#line 279
                                    stepV_0 = _S783;

#line 277
                                }
                                else
                                {

#line 277
                                    float4 _S784 = float4(dLambda_0) ;

#line 282
                                    float4 _S785 = *VB_0 + k3VB_0 * _S784;

#line 282
                                    stepXi_0 = *xiB_0 + k3XiB_0 * _S784;

#line 282
                                    stepV_0 = _S785;

#line 277
                                }

#line 274
                            }

#line 271
                        }

#line 268
                    }

#line 265
                }

#line 262
            }

#line 259
        }

#line 259
        thread array<array<array<float, int(4)>, int(4)>, int(4)> _S786 = Gamma_8;

#line 259
        thread array<array<array<array<float, int(4)>, int(4)>, int(4)>, int(4)> _S787 = R_10;

#line 259
        thread Vec4Cart_0 _S788 = phaseTangent_1;

#line 285
        thread float4 stepDXi_0;

#line 285
        thread float4 stepDV_0;

#line 285
        GeodesicDeviationRhs_0(&_S786, &_S787, &_S788, stepXi_0, stepV_0, &stepDXi_0, &stepDV_0);

#line 285
        float4 k2XiA_1;

#line 285
        float4 k2VA_1;

        if(_S745)
        {

#line 287
            k1XiA_0 = stepDXi_0;

#line 287
            k1VA_0 = stepDV_0;

#line 287
            k2XiA_1 = k2XiA_0;

#line 287
            k2VA_1 = k2VA_0;

#line 287
        }
        else
        {

#line 287
            float4 k1XiB_1;

#line 287
            float4 k1VB_1;

#line 287
            float4 k2XiB_1;

#line 287
            float4 k2VB_1;

#line 287
            float4 k3XiA_1;

#line 287
            float4 k3VA_1;

#line 287
            float4 k3XiB_1;

#line 287
            float4 k3VB_1;

#line 287
            float4 k4XiA_1;

#line 287
            float4 k4VA_1;

#line 287
            float4 k4XiB_1;

#line 287
            float4 k4VB_1;


            if(step_6 == int(1))
            {

#line 290
                k2XiA_1 = k2XiA_0;

#line 290
                k2VA_1 = k2VA_0;

#line 290
                k1XiB_1 = stepDXi_0;

#line 290
                k1VB_1 = stepDV_0;

#line 290
                k2XiB_1 = k2XiB_0;

#line 290
                k2VB_1 = k2VB_0;

#line 290
                k3XiA_1 = k3XiA_0;

#line 290
                k3VA_1 = k3VA_0;

#line 290
                k3XiB_1 = k3XiB_0;

#line 290
                k3VB_1 = k3VB_0;

#line 290
                k4XiA_1 = k4XiA_0;

#line 290
                k4VA_1 = k4VA_0;

#line 290
                k4XiB_1 = k4XiB_0;

#line 290
                k4VB_1 = k4VB_0;

#line 290
            }
            else
            {
                if(step_6 == int(2))
                {

#line 293
                    k2XiA_1 = stepDXi_0;

#line 293
                    k2VA_1 = stepDV_0;

#line 293
                    k1XiB_1 = k2XiB_0;

#line 293
                    k1VB_1 = k2VB_0;

#line 293
                    k2XiB_1 = k3XiA_0;

#line 293
                    k2VB_1 = k3VA_0;

#line 293
                    k3XiA_1 = k3XiB_0;

#line 293
                    k3VA_1 = k3VB_0;

#line 293
                    k3XiB_1 = k4XiA_0;

#line 293
                    k3VB_1 = k4VA_0;

#line 293
                    k4XiA_1 = k4XiB_0;

#line 293
                    k4VA_1 = k4VB_0;

#line 293
                }
                else
                {
                    if(step_6 == int(3))
                    {

#line 296
                        k2XiA_1 = k2XiB_0;

#line 296
                        k2VA_1 = k2VB_0;

#line 296
                        k1XiB_1 = stepDXi_0;

#line 296
                        k1VB_1 = stepDV_0;

#line 296
                        k2XiB_1 = k3XiB_0;

#line 296
                        k2VB_1 = k3VB_0;

#line 296
                        k3XiA_1 = k4XiA_0;

#line 296
                        k3VA_1 = k4VA_0;

#line 296
                        k3XiB_1 = k4XiB_0;

#line 296
                        k3VB_1 = k4VB_0;

#line 296
                    }
                    else
                    {
                        if(step_6 == int(4))
                        {

#line 299
                            k2XiA_1 = stepDXi_0;

#line 299
                            k2VA_1 = stepDV_0;

#line 299
                            k1XiB_1 = k3XiB_0;

#line 299
                            k1VB_1 = k3VB_0;

#line 299
                            k2XiB_1 = k4XiA_0;

#line 299
                            k2VB_1 = k4VA_0;

#line 299
                            k3XiA_1 = k4XiB_0;

#line 299
                            k3VA_1 = k4VB_0;

#line 299
                        }
                        else
                        {
                            if(step_6 == int(5))
                            {

#line 302
                                k2XiA_1 = stepDXi_0;

#line 302
                                k2VA_1 = stepDV_0;

#line 302
                                k1XiB_1 = k4XiA_0;

#line 302
                                k1VB_1 = k4VA_0;

#line 302
                                k2XiB_1 = k4XiB_0;

#line 302
                                k2VB_1 = k4VB_0;

#line 302
                            }
                            else
                            {
                                if(step_6 == int(6))
                                {

#line 305
                                    k2XiA_1 = stepDXi_0;

#line 305
                                    k2VA_1 = stepDV_0;

#line 305
                                    k1XiB_1 = k4XiB_0;

#line 305
                                    k1VB_1 = k4VB_0;

#line 305
                                }
                                else
                                {

#line 305
                                    k2XiA_1 = k4XiA_0;

#line 305
                                    k2VA_1 = k4VA_0;

#line 305
                                    k1XiB_1 = stepDXi_0;

#line 305
                                    k1VB_1 = stepDV_0;

#line 305
                                }

#line 227
                                float4 _S789 = k2XiA_1;

#line 227
                                float4 _S790 = k2VA_1;
                                float4 _S791 = k1XiB_1;

#line 228
                                float4 _S792 = k1VB_1;

#line 228
                                k2XiA_1 = k3XiB_0;

#line 228
                                k2VA_1 = k3VB_0;

#line 228
                                k1XiB_1 = _S789;

#line 228
                                k1VB_1 = _S790;

#line 228
                                k2XiB_1 = _S791;

#line 228
                                k2VB_1 = _S792;

#line 302
                            }

#line 228
                            float4 _S793 = k2XiA_1;

#line 228
                            float4 _S794 = k2VA_1;

#line 227
                            float4 _S795 = k1XiB_1;

#line 227
                            float4 _S796 = k1VB_1;
                            float4 _S797 = k2XiB_1;

#line 228
                            float4 _S798 = k2VB_1;

#line 228
                            k2XiA_1 = k2XiB_0;

#line 228
                            k2VA_1 = k2VB_0;

#line 228
                            k1XiB_1 = _S793;

#line 228
                            k1VB_1 = _S794;

#line 228
                            k2XiB_1 = _S795;

#line 228
                            k2VB_1 = _S796;

#line 228
                            k3XiA_1 = _S797;

#line 228
                            k3VA_1 = _S798;

#line 299
                        }

#line 228
                        float4 _S799 = k1XiB_1;

#line 228
                        float4 _S800 = k1VB_1;

#line 227
                        float4 _S801 = k2XiB_1;

#line 227
                        float4 _S802 = k2VB_1;
                        float4 _S803 = k3XiA_1;

#line 228
                        float4 _S804 = k3VA_1;

#line 228
                        k1XiB_1 = k3XiA_0;

#line 228
                        k1VB_1 = k3VA_0;

#line 228
                        k2XiB_1 = _S799;

#line 228
                        k2VB_1 = _S800;

#line 228
                        k3XiA_1 = _S801;

#line 228
                        k3VA_1 = _S802;

#line 228
                        k3XiB_1 = _S803;

#line 228
                        k3VB_1 = _S804;

#line 296
                    }

#line 228
                    float4 _S805 = k2XiA_1;

#line 228
                    float4 _S806 = k2VA_1;

#line 227
                    float4 _S807 = k1XiB_1;

#line 227
                    float4 _S808 = k1VB_1;
                    float4 _S809 = k2XiB_1;

#line 228
                    float4 _S810 = k2VB_1;

#line 227
                    float4 _S811 = k3XiA_1;

#line 227
                    float4 _S812 = k3VA_1;
                    float4 _S813 = k3XiB_1;

#line 228
                    float4 _S814 = k3VB_1;

#line 228
                    k2XiA_1 = k2XiA_0;

#line 228
                    k2VA_1 = k2VA_0;

#line 228
                    k1XiB_1 = _S805;

#line 228
                    k1VB_1 = _S806;

#line 228
                    k2XiB_1 = _S807;

#line 228
                    k2VB_1 = _S808;

#line 228
                    k3XiA_1 = _S809;

#line 228
                    k3VA_1 = _S810;

#line 228
                    k3XiB_1 = _S811;

#line 228
                    k3VB_1 = _S812;

#line 228
                    k4XiA_1 = _S813;

#line 228
                    k4VA_1 = _S814;

#line 293
                }

#line 228
                float4 _S815 = k1XiB_1;

#line 228
                float4 _S816 = k1VB_1;

#line 227
                float4 _S817 = k2XiB_1;

#line 227
                float4 _S818 = k2VB_1;
                float4 _S819 = k3XiA_1;

#line 228
                float4 _S820 = k3VA_1;

#line 227
                float4 _S821 = k3XiB_1;

#line 227
                float4 _S822 = k3VB_1;
                float4 _S823 = k4XiA_1;

#line 228
                float4 _S824 = k4VA_1;

#line 228
                k1XiB_1 = k1XiB_0;

#line 228
                k1VB_1 = k1VB_0;

#line 228
                k2XiB_1 = _S815;

#line 228
                k2VB_1 = _S816;

#line 228
                k3XiA_1 = _S817;

#line 228
                k3VA_1 = _S818;

#line 228
                k3XiB_1 = _S819;

#line 228
                k3VB_1 = _S820;

#line 228
                k4XiA_1 = _S821;

#line 228
                k4VA_1 = _S822;

#line 228
                k4XiB_1 = _S823;

#line 228
                k4VB_1 = _S824;

#line 290
            }

#line 290
            k1XiB_0 = k1XiB_1;

#line 290
            k1VB_0 = k1VB_1;

#line 290
            k2XiB_0 = k2XiB_1;

#line 290
            k2VB_0 = k2VB_1;

#line 290
            k3XiA_0 = k3XiA_1;

#line 290
            k3VA_0 = k3VA_1;

#line 290
            k3XiB_0 = k3XiB_1;

#line 290
            k3VB_0 = k3VB_1;

#line 290
            k4XiA_0 = k4XiA_1;

#line 290
            k4VA_0 = k4VA_1;

#line 290
            k4XiB_0 = k4XiB_1;

#line 290
            k4VB_0 = k4VB_1;

#line 287
        }

#line 231
        int step_7 = step_6 + int(1);

#line 231
        k2XiA_0 = k2XiA_1;

#line 231
        k2VA_0 = k2VA_1;

#line 231
        phaseTangent_0 = phaseTangent_1;

#line 231
        step_6 = step_7;

#line 231
    }

#line 231
    float4 _S825 = float4(2.0) ;

#line 231
    float4 _S826 = float4((dLambda_0 / 6.0)) ;

#line 315
    *xiA_0 = *xiA_0 + (k1XiA_0 + k2XiA_0 * _S825 + k3XiA_0 * _S825 + k4XiA_0) * _S826;
    *VA_0 = *VA_0 + (k1VA_0 + k2VA_0 * _S825 + k3VA_0 * _S825 + k4VA_0) * _S826;
    *xiB_0 = *xiB_0 + (k1XiB_0 + k2XiB_0 * _S825 + k3XiB_0 * _S825 + k4XiB_0) * _S826;
    *VB_0 = *VB_0 + (k1VB_0 + k2VB_0 * _S825 + k3VB_0 * _S825 + k4VB_0) * _S826;
    return;
}


#line 31 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_ellis_topology.slang"
bool MapEllisSecondSheetSkyDirection_0(const Vec4Cart_0 thread* position_2, float3 chartDirection_0, float3 thread* mappedDirection_0)
{
    *mappedDirection_0 = float3(0.0, 0.0, 1.0);

#line 33
    float _S827 = position_2->x_2;

#line 33
    float _S828 = position_2->y_0;

#line 33
    float _S829 = position_2->z_0;
    float radius_2 = sqrt(position_2->x_2 * position_2->x_2 + position_2->y_0 * position_2->y_0 + position_2->z_0 * position_2->z_0);

#line 34
    bool _S830;

    if(isfinite(radius_2))
    {

#line 36
        _S830 = radius_2 > 0.0;

#line 36
    }
    else
    {

#line 36
        _S830 = false;

#line 36
    }

#line 36
    if(!_S830)
    {

#line 36
        _S830 = true;

#line 36
    }
    else
    {

#line 36
        _S830 = !all(isfinite(chartDirection_0));

#line 36
    }

#line 36
    if(_S830)
    {

#line 37
        return false;
    }

#line 38
    float3 normal_0 = float3(_S827, _S828, _S829) / float3(radius_2) ;
    float3 direction_4 = float3(chartDirection_0.x, chartDirection_0.y, chartDirection_0.z);

    float norm_2 = sqrt(dot(direction_4, direction_4));
    if(isfinite(norm_2))
    {

#line 42
        _S830 = norm_2 > 0.0;

#line 42
    }
    else
    {

#line 42
        _S830 = false;

#line 42
    }

#line 42
    if(!_S830)
    {

#line 42
        return false;
    }

#line 43
    float3 reflected_0 = direction_4 - float3((2.0 * dot(normal_0, direction_4)))  * normal_0;
    float3 _S831 = float3(reflected_0.x, reflected_0.y, reflected_0.z);

#line 44
    *mappedDirection_0 = _S831;
    return all(isfinite(_S831));
}


#line 377 "/home/astra/.project/Sirius/out/trace-jacobi-shared-segment-rhs-screen/inputs/candidate/gr_deviation.slang"
float BeamScaleExponent_0(float value_4, int exponent_0)
{

#line 378
    if(value_4 == 0.0)
    {

#line 378
        return value_4;
    }

#line 378
    float _S832 = value_4;

#line 378
    int _S833 = exponent_0;
    for(;;)
    {

#line 379
        if(_S833 > int(120))
        {
        }
        else
        {

#line 379
            break;
        }
        int _S834 = _S833 - int(120);

#line 381
        _S832 = ldexp_0(_S832, 120.0);

#line 381
        _S833 = _S834;

#line 379
    }



    for(;;)
    {

#line 383
        if(_S833 < int(-120))
        {
        }
        else
        {

#line 383
            break;
        }
        int _S835 = _S833 + int(120);

#line 385
        _S832 = ldexp_0(_S832, -120.0);

#line 385
        _S833 = _S835;

#line 383
    }



    return ldexp_0(_S832, float(_S833));
}


#line 362
float BeamProductResidual_0(float a_27, float b_1, float roundedProduct_0)
{

#line 363
    precise float sa_0 = 4097.0 * a_27;

#line 363
    precise float sb_0 = 4097.0 * b_1;
    precise float ah_0 = sa_0 - (sa_0 - a_27);

#line 364
    precise float bh_0 = sb_0 - (sb_0 - b_1);
    precise float al_0 = a_27 - ah_0;

#line 365
    precise float bl_0 = b_1 - bh_0;
    precise float e1_0 = roundedProduct_0 - ah_0 * bh_0;
    precise float e2_0 = e1_0 - al_0 * bh_0;
    precise float e3_0 = e2_0 - ah_0 * bl_0;
    precise float residual_0 = al_0 * bl_0 - e3_0;
    return residual_0;
}


#line 390
void BeamEllipseFromTransverseMap_0(float beamA_0, float beamB_0, float beamC_0, float beamD_0, float thread* semiMajor_0, float thread* semiMinor_0, float thread* orientation_0, float thread* determinant_0)
{

    float _S836 = max(max(abs(beamA_0), abs(beamB_0)), max(abs(beamC_0), abs(beamD_0)));
    if(_S836 == 0.0)
    {

#line 395
        *determinant_0 = 0.0;

#line 395
        *orientation_0 = 0.0;

#line 395
        *semiMinor_0 = 0.0;

#line 395
        *semiMajor_0 = 0.0;
        return;
    }
    thread int exponent_1;
    float _S837 = (frexp((_S836), *((&exponent_1))));
    float a_28 = BeamScaleExponent_0(beamA_0, - exponent_1);

#line 400
    float b_2 = BeamScaleExponent_0(beamB_0, - exponent_1);
    float c_0 = BeamScaleExponent_0(beamC_0, - exponent_1);

#line 401
    float d_0 = BeamScaleExponent_0(beamD_0, - exponent_1);


    thread int ae_0;
    float af_0 = (frexp((beamA_0), *((&ae_0))));

#line 404
    thread int be_0;
    float bf_0 = (frexp((beamB_0), *((&be_0))));

#line 404
    thread int ce_0;

    float cf_0 = (frexp((beamC_0), *((&ce_0))));

#line 404
    thread int de_0;

    float df_0 = (frexp((beamD_0), *((&de_0))));

#line 406
    bool hasAD_0;
    if(beamA_0 != 0.0)
    {

#line 407
        hasAD_0 = beamD_0 != 0.0;

#line 407
    }
    else
    {

#line 407
        hasAD_0 = false;

#line 407
    }

#line 407
    bool hasBC_0;
    if(beamB_0 != 0.0)
    {

#line 408
        hasBC_0 = beamC_0 != 0.0;

#line 408
    }
    else
    {

#line 408
        hasBC_0 = false;

#line 408
    }

#line 408
    bool _S838;
    if(hasAD_0)
    {

#line 409
        _S838 = hasBC_0;

#line 409
    }
    else
    {

#line 409
        _S838 = false;

#line 409
    }

#line 409
    int determinantExponent_0;

#line 409
    if(_S838)
    {

#line 409
        determinantExponent_0 = max(ae_0 + de_0, be_0 + ce_0);

#line 409
    }
    else
    {

#line 410
        if(hasAD_0)
        {

#line 410
            determinantExponent_0 = ae_0 + de_0;

#line 410
        }
        else
        {

#line 410
            if(hasBC_0)
            {

#line 410
                determinantExponent_0 = be_0 + ce_0;

#line 410
            }
            else
            {

#line 410
                determinantExponent_0 = int(0);

#line 410
            }

#line 410
        }

#line 409
    }

    precise float ad_0 = af_0 * df_0;

#line 411
    precise float bc_0 = bf_0 * cf_0;


    precise float normalizedDeterminant_0 = BeamScaleExponent_0(ad_0, ae_0 + de_0 - determinantExponent_0) - BeamScaleExponent_0(bc_0, be_0 + ce_0 - determinantExponent_0) + (BeamScaleExponent_0(BeamProductResidual_0(af_0, df_0, ad_0), ae_0 + de_0 - determinantExponent_0) - BeamScaleExponent_0(BeamProductResidual_0(bf_0, cf_0, bc_0), be_0 + ce_0 - determinantExponent_0));


    float normalizedMajor_0 = 0.5 * (length(float2(a_28 + d_0, b_2 - c_0)) + length(float2(a_28 - d_0, b_2 + c_0)));

    *semiMajor_0 = BeamScaleExponent_0(normalizedMajor_0, exponent_1);


    *semiMinor_0 = BeamScaleExponent_0(abs(normalizedDeterminant_0) / normalizedMajor_0, determinantExponent_0 - exponent_1);

    *orientation_0 = 0.5 * atan2(2.0 * (a_28 * c_0 + b_2 * d_0), a_28 * a_28 + b_2 * b_2 - c_0 * c_0 - d_0 * d_0);
    *determinant_0 = BeamScaleExponent_0(normalizedDeterminant_0, determinantExponent_0);
    return;
}


#line 6 "/home/astra/.project/Sirius/src/sirius/backend/../kernels/gr_point_star.slang"
float PointStarAngularWeight_0(float3 direction_5, float3 star_0, float3 tangentA_0, float3 tangentB_0, float cutoff_0, float cosOrientation_0, float sinOrientation_0, float invMajorSquared_0, float invMinorSquared_0)
{


    float normSquared_0 = dot(direction_5, direction_5);

#line 10
    bool _S839;
    if(normSquared_0 < 0.25)
    {

#line 11
        _S839 = true;

#line 11
    }
    else
    {

#line 11
        _S839 = normSquared_0 > 4.0;

#line 11
    }

#line 11
    float3 _S840;

#line 11
    if(_S839)
    {

#line 11
        _S840 = direction_5 / float3(sqrt(normSquared_0)) ;

#line 11
    }
    else
    {

#line 11
        _S840 = direction_5;

#line 11
    }


    float3 normal_1 = cross(_S840, star_0 - _S840);
    float sine_0 = length(normal_1);
    float angle_0 = atan2(sine_0, dot(_S840, star_0));
    if(angle_0 > cutoff_0)
    {

#line 17
        return 0.0;
    }

#line 18
    thread float2 offset_1 = float2(0.0) ;
    if(sine_0 > 0.0)
    {

#line 20
        float3 tangent_3 = cross(normal_1, normalize(_S840)) / float3(sine_0) ;
        offset_1 = float2(angle_0)  * float2(dot(tangent_3, tangentA_0), dot(tangent_3, tangentB_0));

#line 19
    }
    else
    {
        if(angle_0 > 0.0)
        {

            if(invMajorSquared_0 != invMinorSquared_0)
            {

#line 25
                return 0.0;
            }

#line 26
            offset_1.x = angle_0;

#line 22
        }

#line 19
    }

#line 28
    float alongMajor_0 = cosOrientation_0 * offset_1.x + sinOrientation_0 * offset_1.y;
    float alongMinor_0 = - sinOrientation_0 * offset_1.x + cosOrientation_0 * offset_1.y;
    return exp(-0.5 * (alongMajor_0 * alongMajor_0 * invMajorSquared_0 + alongMinor_0 * alongMinor_0 * invMinorSquared_0));
}


#line 127 "/home/astra/.project/Sirius/out/trace-jacobi-shared-segment-rhs-screen/inputs/candidate/trace.slang"
float3 PointStarColor_0(uint starIndex_0, KernelContext_0 thread* kernelContext_1)
{

    return BlackbodyColor_0(kernelContext_1->pointStars_0[starIndex_0 * 8U + 6U]);
}

void AccumulatePointStarCell_0(int thetaBin_0, int phiBin_0, float3 dir_0, float3 tangentA_1, float3 tangentB_1, float cutoff_1, float cosOrientation_1, float sinOrientation_1, float invMajorSquared_1, float invMinorSquared_1, float brightnessScale_0, float3 thread* color_0, KernelContext_0 thread* kernelContext_2)
{

#line 139
    uint cell_0 = uint(thetaBin_0 * int(512) + phiBin_0);

#line 139
    uint at_0 = kernelContext_2->pointOffsets_0[cell_0];
    for(;;)
    {

#line 140
        if(at_0 < (kernelContext_2->pointOffsets_0)[cell_0 + 1U])
        {
        }
        else
        {

#line 140
            break;
        }

#line 141
        uint starIndex_1 = kernelContext_2->pointIndices_0[at_0];
        uint base_1 = starIndex_1 * 8U;


        float weight_0 = PointStarAngularWeight_0(dir_0, float3(kernelContext_2->pointStars_0[base_1], kernelContext_2->pointStars_0[base_1 + 1U], kernelContext_2->pointStars_0[base_1 + 2U]), tangentA_1, tangentB_1, cutoff_1, cosOrientation_1, sinOrientation_1, invMajorSquared_1, invMinorSquared_1);


        if(weight_0 == 0.0)
        {

#line 148
            at_0 = at_0 + 1U;

#line 140
            continue;
        }

#line 150
        float intensity_0 = pow(10.0, -0.40000000596046448 * kernelContext_2->pointStars_0[base_1 + 4U]) * weight_0 * brightnessScale_0;

#line 150
        float3 _S841 = PointStarColor_0(starIndex_1, kernelContext_2);
        *color_0 = *color_0 + _S841 * float3(intensity_0) ;

#line 140
        at_0 = at_0 + 1U;

#line 140
    }

#line 153
    return;
}




float3 SamplePointStarfield_0(float3 direction_6, float sigmaMajor_0, float sigmaMinor_0, float orientation_1, float brightnessScale_1, KernelContext_0 thread* kernelContext_3)
{


    float lengthDirection_0 = length(direction_6);

#line 163
    bool _S842;
    if(lengthDirection_0 < 9.99999968265522539e-21)
    {

#line 164
        _S842 = true;

#line 164
    }
    else
    {

#line 164
        _S842 = sigmaMajor_0 <= 0.0;

#line 164
    }

#line 164
    if(_S842)
    {

#line 164
        _S842 = true;

#line 164
    }
    else
    {

#line 164
        _S842 = sigmaMinor_0 <= 0.0;

#line 164
    }

#line 164
    if(_S842)
    {

#line 165
        return float3(0.0) ;
    }

#line 166
    float3 dir_1 = direction_6 / float3(lengthDirection_0) ;
    thread float3 tangentA_2;
    thread float3 tangentB_2;
    MakeCelestialTangentBasis_0(direction_6, &tangentA_2, &tangentB_2);
    float theta_0 = atan2(length(dir_1.xy), dir_1.z);
    float phi_1 = atan2(dir_1.y, dir_1.x);

#line 171
    float phi_2;
    if(phi_1 < 0.0)
    {

#line 172
        phi_2 = phi_1 + 6.28318548202514648;

#line 172
    }
    else
    {

#line 172
        phi_2 = phi_1;

#line 172
    }
    float _S843 = max(sigmaMajor_0, sigmaMinor_0);
    float _S844 = min(sigmaMajor_0, sigmaMinor_0);
    float _S845 = min(4.0 * _S843, 3.14159274101257324);
    int firstTheta_0 = clamp(int(floor(max(theta_0 - _S845, 0.0) / 3.14159274101257324 * 256.0)) - int(1), int(0), int(255));

    int _S846 = clamp(int(floor(min(theta_0 + _S845, 3.14159274101257324) / 3.14159274101257324 * 256.0)) + int(1), int(0), int(255));



    float sinTheta_0 = sin(theta_0);
    if(theta_0 > _S845)
    {

#line 183
        _S842 = (3.14159274101257324 - theta_0) > _S845;

#line 183
    }
    else
    {

#line 183
        _S842 = false;

#line 183
    }

#line 183
    if(_S842)
    {

#line 183
        _S842 = sinTheta_0 > 9.99999993922529029e-09;

#line 183
    }
    else
    {

#line 183
        _S842 = false;

#line 183
    }

#line 183
    float phiHalfWidth_0;

#line 183
    if(_S842)
    {

#line 183
        phiHalfWidth_0 = asin(clamp(sin(_S845) / sinTheta_0, 0.0, 1.0));

#line 183
    }
    else
    {

#line 183
        phiHalfWidth_0 = 3.14159274101257324;

#line 183
    }



    int _S847 = min(int(floor(phi_2 / 6.28318548202514648 * 512.0)), int(511));
    int _S848 = min(int(ceil(phiHalfWidth_0 / 0.01227184664458036)) + int(1), int(512));
    bool _S849 = (int(2) * _S848 + int(1)) >= int(512);
    float _S850 = 1.0 / (_S843 * _S843);
    float _S851 = 1.0 / (_S844 * _S844);
    float _S852 = cos(orientation_1);
    float _S853 = sin(orientation_1);
    thread float3 color_1 = float3(0.0) ;

#line 194
    int thetaBin_1 = firstTheta_0;

    for(;;)
    {

#line 196
        if(thetaBin_1 <= _S846)
        {
        }
        else
        {

#line 196
            break;
        }

#line 196
        int phiBin_1;
        if(_S849)
        {

#line 197
            phiBin_1 = int(0);
            for(;;)
            {

#line 198
                if(phiBin_1 < int(512))
                {
                }
                else
                {

#line 198
                    break;
                }

#line 198
                AccumulatePointStarCell_0(thetaBin_1, phiBin_1, direction_6, tangentA_2, tangentB_2, _S845, _S852, _S853, _S850, _S851, brightnessScale_1, &color_1, kernelContext_3);

#line 198
                phiBin_1 = phiBin_1 + int(1);

#line 198
            }

#line 197
        }
        else
        {

#line 197
            int delta_0 = - _S848;

#line 205
            for(;;)
            {

#line 205
                if(delta_0 <= _S848)
                {
                }
                else
                {

#line 205
                    break;
                }

#line 206
                int phiBin_2 = (_S847 + delta_0) % int(512);
                if(phiBin_2 < int(0))
                {

#line 207
                    phiBin_1 = phiBin_2 + int(512);

#line 207
                }
                else
                {

#line 207
                    phiBin_1 = phiBin_2;

#line 207
                }

#line 207
                AccumulatePointStarCell_0(thetaBin_1, phiBin_1, direction_6, tangentA_2, tangentB_2, _S845, _S852, _S853, _S850, _S851, brightnessScale_1, &color_1, kernelContext_3);

#line 205
                delta_0 = delta_0 + int(1);

#line 205
            }

#line 197
        }

#line 196
        thetaBin_1 = thetaBin_1 + int(1);

#line 196
    }

#line 215
    return color_1;
}


#line 77
float3 StarfieldTexel_0(int x_9, int y_5, int width_1, KernelContext_0 thread* kernelContext_4)
{

#line 78
    uint packed_0 = kernelContext_4->starfield_0[y_5 * width_1 + x_9];



    return float3(float(packed_0 & 255U) / 255.0, float((packed_0 >> 8U) & 255U) / 255.0, float((packed_0 >> 16U) & 255U) / 255.0);
}


#line 89
float3 SampleStarfield_0(float3 dir_2, int width_2, int height_1, KernelContext_0 thread* kernelContext_5)
{

#line 90
    float len_0 = length(dir_2);
    if(len_0 < 1.00000001335143196e-10)
    {

#line 92
        return float3(0.00100000004749745, 0.00100000004749745, 0.00100000004749745);
    }
    float3 d_1 = dir_2 / float3(len_0) ;

    float theta_1 = acos(clamp(d_1.z, -1.0, 1.0));
    float phi_3 = atan2(d_1.y, d_1.x);

#line 97
    float phi_4;
    if(phi_3 < 0.0)
    {

#line 98
        phi_4 = phi_3 + 6.28318548202514648;

#line 98
    }
    else
    {

#line 98
        phi_4 = phi_3;

#line 98
    }

#line 103
    int _S854 = width_2 - int(1);

#line 103
    float px_0 = phi_4 / 6.28318548202514648 * float(_S854);
    int _S855 = height_1 - int(1);

#line 104
    float py_0 = theta_1 / 3.14159274101257324 * float(_S855);

    int x0_0 = int(floor(px_0));
    int y0_0 = int(floor(py_0));
    int _S856 = min(x0_0 + int(1), _S854);
    int _S857 = min(y0_0 + int(1), _S855);

    float fx_0 = px_0 - float(x0_0);
    float fy_0 = py_0 - float(y0_0);

#line 112
    float3 _S858 = StarfieldTexel_0(x0_0, y0_0, width_2, kernelContext_5);

#line 112
    float3 _S859 = StarfieldTexel_0(_S856, y0_0, width_2, kernelContext_5);

#line 112
    float3 _S860 = StarfieldTexel_0(x0_0, _S857, width_2, kernelContext_5);

#line 112
    float3 _S861 = StarfieldTexel_0(_S856, _S857, width_2, kernelContext_5);

#line 119
    float _S862 = 1.0 - fx_0;

#line 119
    float _S863 = 1.0 - fy_0;

#line 124
    return _S858 * float3((_S862 * _S863))  + _S859 * float3((fx_0 * _S863))  + _S860 * float3((_S862 * fy_0))  + _S861 * float3((fx_0 * fy_0)) ;
}


#line 71
float3 SampleBackgroundGradient_0(float3 dir_3)
{

#line 72
    return float3(0.5, 0.5, 0.5) + float3(0.5)  * dir_3;
}


#line 420
[[kernel]] void ComputeMain(uint3 tid_0 [[thread_position_in_grid]], float device* params_1 [[buffer(1)]], TraceContinuation_natural_0 device* continuation_3 [[buffer(6)]], float device* radiance_1 [[buffer(0)]], uint device* pointOffsets_1 [[buffer(4)]], uint device* pointIndices_1 [[buffer(5)]], float device* pointStars_1 [[buffer(3)]], uint device* starfield_1 [[buffer(2)]])
{

#line 420
    int right_2;

#line 420
    int rootCount_2;

#line 420
    thread KernelContext_0 kernelContext_6;

#line 420
    (&kernelContext_6)->params_0 = params_1;

#line 420
    (&kernelContext_6)->continuation_0 = continuation_3;

#line 420
    (&kernelContext_6)->radiance_0 = radiance_1;

#line 420
    (&kernelContext_6)->pointOffsets_0 = pointOffsets_1;

#line 420
    (&kernelContext_6)->pointIndices_0 = pointIndices_1;

#line 420
    (&kernelContext_6)->pointStars_0 = pointStars_1;

#line 420
    (&kernelContext_6)->starfield_0 = starfield_1;
    int imageWidth_0 = int(params_1[int(0)]);
    int imageHeight_0 = int(params_1[int(1)]);
    uint metricId_0 = uint(params_1[int(2)]);

    int tileOriginX_0 = int(params_1[int(31)]);
    int tileOriginY_0 = int(params_1[int(32)]);
    int tileWidth_0 = int(params_1[int(33)]);
    int tileHeight_0 = int(params_1[int(34)]);

    int tx_0 = int(tid_0.x);
    int ty_0 = int(tid_0.y);

#line 431
    bool escaped_0;
    if(tx_0 >= tileWidth_0)
    {

#line 432
        escaped_0 = true;

#line 432
    }
    else
    {

#line 432
        escaped_0 = ty_0 >= tileHeight_0;

#line 432
    }

#line 432
    if(escaped_0)
    {

#line 432
        return;
    }
    int ix_0 = tileOriginX_0 + tx_0;
    int iy_0 = tileOriginY_0 + ty_0;
    if(ix_0 >= imageWidth_0)
    {

#line 436
        escaped_0 = true;

#line 436
    }
    else
    {

#line 436
        escaped_0 = iy_0 >= imageHeight_0;

#line 436
    }

#line 436
    if(escaped_0)
    {

#line 436
        return;
    }
    float M_20 = (&kernelContext_6)->params_0[int(3)];
    float a_29 = (&kernelContext_6)->params_0[int(4)];
    float Q_13 = (&kernelContext_6)->params_0[int(5)];
    float Lambda_13 = (&kernelContext_6)->params_0[int(6)];

#line 441
    float3 _S864 = P3_0(7U, &kernelContext_6);

#line 441
    float3 _S865 = P3_0(10U, &kernelContext_6);

#line 441
    float3 _S866 = P3_0(13U, &kernelContext_6);

#line 441
    float3 _S867 = P3_0(16U, &kernelContext_6);

#line 447
    float fov_0 = (&kernelContext_6)->params_0[int(19)];
    float aspect_0 = (&kernelContext_6)->params_0[int(20)];

    int maxSteps_0 = int((&kernelContext_6)->params_0[int(21)]);
    float stepScale_3 = (&kernelContext_6)->params_0[int(22)];
    float minStep_3 = (&kernelContext_6)->params_0[int(23)];
    float maxStep_3 = (&kernelContext_6)->params_0[int(24)];
    float escapeRadius_4 = (&kernelContext_6)->params_0[int(25)];
    float captureFactor_0 = (&kernelContext_6)->params_0[int(26)];

    bool _S868 = (&kernelContext_6)->params_0[int(27)] > 0.5;
    float _S869 = (&kernelContext_6)->params_0[int(28)];
    float _S870 = (&kernelContext_6)->params_0[int(29)];
    float _S871 = (&kernelContext_6)->params_0[int(30)];

    bool starfieldEnabled_0 = (&kernelContext_6)->params_0[int(35)] > 0.5;
    int starfieldWidth_0 = int((&kernelContext_6)->params_0[int(36)]);
    int starfieldHeight_0 = int((&kernelContext_6)->params_0[int(37)]);

    float warpVs_5 = (&kernelContext_6)->params_0[int(38)];
    float warpSigma_5 = (&kernelContext_6)->params_0[int(39)];
    float warpR_5 = (&kernelContext_6)->params_0[int(40)];

#line 473
    float throatB0_5 = (&kernelContext_6)->params_0[int(42)];
    int wormShape_8 = int((&kernelContext_6)->params_0[int(43)]);
    float sampleU_1 = (&kernelContext_6)->params_0[int(44)];
    float sampleV_1 = (&kernelContext_6)->params_0[int(45)];
    float sampleIndex_0 = (&kernelContext_6)->params_0[int(46)];

#line 477
    float3 _S872 = P3_0(47U, &kernelContext_6);

    int _S873 = int((&kernelContext_6)->params_0[int(50)]);
    bool _S874 = (&kernelContext_6)->params_0[int(51)] > 0.5;
    int lensModel_0 = int((&kernelContext_6)->params_0[int(52)]);
    float focalLength_1 = (&kernelContext_6)->params_0[int(53)];
    float aperture_1 = (&kernelContext_6)->params_0[int(54)];
    float focusDistance_1 = (&kernelContext_6)->params_0[int(55)];
    bool pointStarfieldEnabled_0 = (&kernelContext_6)->params_0[int(56)] > 0.5;
    float pixelAngularSize_0 = (&kernelContext_6)->params_0[int(57)];
    float pointBrightnessScale_0 = (&kernelContext_6)->params_0[int(58)];
    bool _S875 = (&kernelContext_6)->params_0[int(59)] > 0.5;
    float _S876 = (&kernelContext_6)->params_0[int(60)];
    float _S877 = (&kernelContext_6)->params_0[int(61)];
    float _S878 = (&kernelContext_6)->params_0[int(62)];
    int _S879 = int((&kernelContext_6)->params_0[int(63)]);
    bool _S880 = (&kernelContext_6)->params_0[int(64)] > 0.5;
    bool twoSheetWormhole_0 = (&kernelContext_6)->params_0[int(65)] > 0.5;

#line 501
    bool beamsEnabled_0 = (&kernelContext_6)->params_0[int(41)] > 0.5;

    int _S881 = ty_0 * tileWidth_0 + tx_0;

#line 503
    uint slot_0 = uint(_S881);
    uint action_0 = uint((&kernelContext_6)->params_0[int(68)]);
    uint4 identity_1 = uint4(uint(ix_0), uint(iy_0), uint(sampleIndex_0), (as_type<uint>(((&kernelContext_6)->params_0[int(69)]))));
    thread TraceContinuation_0 record_1;
    bool _S882 = action_0 == 0U;

#line 507
    if(_S882)
    {

#line 507
        float4 _S883 = float4(0.0) ;

#line 507
        array<float4, int(4)> _S884 = { _S883, _S883, _S883, _S883 };

#line 507
        uint4 _S885 = uint4(0U) ;

#line 507
        (&record_1)->x_1 = _S883;

#line 507
        (&record_1)->u_0 = _S883;

#line 507
        (&record_1)->X_0 = _S884;

#line 507
        (&record_1)->V_0 = _S884;

#line 507
        (&record_1)->compX_0 = _S883;

#line 507
        (&record_1)->compU_0 = _S883;

#line 507
        (&record_1)->integration_0 = _S883;

#line 507
        (&record_1)->source_0 = _S883;

#line 507
        (&record_1)->volume_0 = _S883;

#line 507
        (&record_1)->control_0 = _S885;

#line 507
        (&record_1)->identity_0 = _S885;

#line 507
    }
    else
    {

#line 507
        TraceContinuation_natural_0 device* _S886 = (&kernelContext_6)->continuation_0+slot_0;

#line 507
        float4 _S887 = float4((*_S886).u_0) ;

#line 507
        array<float4, int(4)> _S888 = { float4((*_S886).X_0.data_0[int(0)]) , float4((*_S886).X_0.data_0[int(1)]) , float4((*_S886).X_0.data_0[int(2)]) , float4((*_S886).X_0.data_0[int(3)])  };

#line 507
        array<float4, int(4)> _S889 = { float4((*_S886).V_0.data_0[int(0)]) , float4((*_S886).V_0.data_0[int(1)]) , float4((*_S886).V_0.data_0[int(2)]) , float4((*_S886).V_0.data_0[int(3)])  };

#line 507
        float4 _S890 = float4((*_S886).compX_0) ;

#line 507
        float4 _S891 = float4((*_S886).compU_0) ;

#line 507
        float4 _S892 = float4((*_S886).integration_0) ;

#line 507
        float4 _S893 = float4((*_S886).source_0) ;

#line 507
        float4 _S894 = float4((*_S886).volume_0) ;

#line 507
        uint4 _S895 = uint4((*_S886).control_0) ;

#line 507
        uint4 _S896 = uint4((*_S886).identity_0) ;

#line 507
        (&record_1)->x_1 = float4((*_S886).x_1) ;

#line 507
        (&record_1)->u_0 = _S887;

#line 507
        (&record_1)->X_0 = _S888;

#line 507
        (&record_1)->V_0 = _S889;

#line 507
        (&record_1)->compX_0 = _S890;

#line 507
        (&record_1)->compU_0 = _S891;

#line 507
        (&record_1)->integration_0 = _S892;

#line 507
        (&record_1)->source_0 = _S893;

#line 507
        (&record_1)->volume_0 = _S894;

#line 507
        (&record_1)->control_0 = _S895;

#line 507
        (&record_1)->identity_0 = _S896;

#line 507
    }

    if((&kernelContext_6)->params_0[int(68)] == 0.0)
    {

#line 509
        escaped_0 = true;

#line 509
    }
    else
    {

#line 509
        escaped_0 = (&kernelContext_6)->params_0[int(68)] == 1.0;

#line 509
    }

#line 509
    if(escaped_0)
    {

#line 509
        escaped_0 = true;

#line 509
    }
    else
    {

#line 509
        escaped_0 = (&kernelContext_6)->params_0[int(68)] == 2.0;

#line 509
    }

#line 509
    if(!escaped_0)
    {

#line 510
        (&record_1)->control_0.x = 2U;
        (&record_1)->control_0.y = 7U;

#line 511
        TraceContinuation_natural_0 device* _S897 = (&kernelContext_6)->continuation_0+slot_0;

#line 511
        _S897->x_1 = packed_float4(record_1.x_1) ;

#line 511
        _S897->u_0 = packed_float4(record_1.u_0) ;

#line 511
        (&_S897->X_0)->data_0[int(0)] = packed_float4(record_1.X_0[int(0)]) ;

#line 511
        (&_S897->X_0)->data_0[int(1)] = packed_float4(record_1.X_0[int(1)]) ;

#line 511
        (&_S897->X_0)->data_0[int(2)] = packed_float4(record_1.X_0[int(2)]) ;

#line 511
        (&_S897->X_0)->data_0[int(3)] = packed_float4(record_1.X_0[int(3)]) ;

#line 511
        (&_S897->V_0)->data_0[int(0)] = packed_float4(record_1.V_0[int(0)]) ;

#line 511
        (&_S897->V_0)->data_0[int(1)] = packed_float4(record_1.V_0[int(1)]) ;

#line 511
        (&_S897->V_0)->data_0[int(2)] = packed_float4(record_1.V_0[int(2)]) ;

#line 511
        (&_S897->V_0)->data_0[int(3)] = packed_float4(record_1.V_0[int(3)]) ;

#line 511
        _S897->compX_0 = packed_float4(record_1.compX_0) ;

#line 511
        _S897->compU_0 = packed_float4(record_1.compU_0) ;

#line 511
        _S897->integration_0 = packed_float4(record_1.integration_0) ;

#line 511
        _S897->source_0 = packed_float4(record_1.source_0) ;

#line 511
        _S897->volume_0 = packed_float4(record_1.volume_0) ;

#line 511
        _S897->control_0 = packed_uint4(record_1.control_0) ;

#line 511
        _S897->identity_0 = packed_uint4(record_1.identity_0) ;

        return;
    }
    if(action_0 != 0U)
    {

#line 515
        if(any(((&record_1)->identity_0) != identity_1))
        {

#line 515
            escaped_0 = true;

#line 515
        }
        else
        {

#line 515
            escaped_0 = ((&record_1)->control_0.x) < 1U;

#line 515
        }
        if(escaped_0)
        {

#line 516
            escaped_0 = true;

#line 516
        }
        else
        {

#line 516
            escaped_0 = ((&record_1)->control_0.x) > 3U;

#line 516
        }

#line 515
    }
    else
    {

#line 515
        escaped_0 = false;

#line 515
    }

#line 515
    if(escaped_0)
    {

        (&record_1)->control_0.x = 2U;
        (&record_1)->control_0.y = 8U;

#line 519
        TraceContinuation_natural_0 device* _S898 = (&kernelContext_6)->continuation_0+slot_0;

#line 519
        _S898->x_1 = packed_float4(record_1.x_1) ;

#line 519
        _S898->u_0 = packed_float4(record_1.u_0) ;

#line 519
        (&_S898->X_0)->data_0[int(0)] = packed_float4(record_1.X_0[int(0)]) ;

#line 519
        (&_S898->X_0)->data_0[int(1)] = packed_float4(record_1.X_0[int(1)]) ;

#line 519
        (&_S898->X_0)->data_0[int(2)] = packed_float4(record_1.X_0[int(2)]) ;

#line 519
        (&_S898->X_0)->data_0[int(3)] = packed_float4(record_1.X_0[int(3)]) ;

#line 519
        (&_S898->V_0)->data_0[int(0)] = packed_float4(record_1.V_0[int(0)]) ;

#line 519
        (&_S898->V_0)->data_0[int(1)] = packed_float4(record_1.V_0[int(1)]) ;

#line 519
        (&_S898->V_0)->data_0[int(2)] = packed_float4(record_1.V_0[int(2)]) ;

#line 519
        (&_S898->V_0)->data_0[int(3)] = packed_float4(record_1.V_0[int(3)]) ;

#line 519
        _S898->compX_0 = packed_float4(record_1.compX_0) ;

#line 519
        _S898->compU_0 = packed_float4(record_1.compU_0) ;

#line 519
        _S898->integration_0 = packed_float4(record_1.integration_0) ;

#line 519
        _S898->source_0 = packed_float4(record_1.source_0) ;

#line 519
        _S898->volume_0 = packed_float4(record_1.volume_0) ;

#line 519
        _S898->control_0 = packed_uint4(record_1.control_0) ;

#line 519
        _S898->identity_0 = packed_uint4(record_1.identity_0) ;

        return;
    }


    bool _S899 = action_0 == 1U;

#line 525
    if(_S899)
    {

#line 525
        escaped_0 = ((&record_1)->control_0.x) != 1U;

#line 525
    }
    else
    {

#line 525
        escaped_0 = false;

#line 525
    }

#line 525
    if(escaped_0)
    {

#line 525
        return;
    }

#line 526
    bool _S900 = action_0 == 2U;

#line 526
    if(_S900)
    {

#line 526
        escaped_0 = ((&record_1)->control_0.x) == 3U;

#line 526
    }
    else
    {

#line 526
        escaped_0 = false;

#line 526
    }

#line 526
    if(escaped_0)
    {

#line 526
        return;
    }

#line 527
    if(_S900)
    {

#line 527
        if(((&record_1)->control_0.x) != 2U)
        {

#line 527
            escaped_0 = true;

#line 527
        }
        else
        {

#line 527
            escaped_0 = !TraceReasonFinalizable_0((&record_1)->control_0.y);

#line 527
        }

#line 527
    }
    else
    {

#line 527
        escaped_0 = false;

#line 527
    }

#line 527
    if(escaped_0)
    {
        if(((&record_1)->control_0.x) == 1U)
        {

#line 530
            (&record_1)->control_0.x = 2U;
            (&record_1)->control_0.y = 8U;

#line 531
            TraceContinuation_natural_0 device* _S901 = (&kernelContext_6)->continuation_0+slot_0;

#line 531
            _S901->x_1 = packed_float4(record_1.x_1) ;

#line 531
            _S901->u_0 = packed_float4(record_1.u_0) ;

#line 531
            (&_S901->X_0)->data_0[int(0)] = packed_float4(record_1.X_0[int(0)]) ;

#line 531
            (&_S901->X_0)->data_0[int(1)] = packed_float4(record_1.X_0[int(1)]) ;

#line 531
            (&_S901->X_0)->data_0[int(2)] = packed_float4(record_1.X_0[int(2)]) ;

#line 531
            (&_S901->X_0)->data_0[int(3)] = packed_float4(record_1.X_0[int(3)]) ;

#line 531
            (&_S901->V_0)->data_0[int(0)] = packed_float4(record_1.V_0[int(0)]) ;

#line 531
            (&_S901->V_0)->data_0[int(1)] = packed_float4(record_1.V_0[int(1)]) ;

#line 531
            (&_S901->V_0)->data_0[int(2)] = packed_float4(record_1.V_0[int(2)]) ;

#line 531
            (&_S901->V_0)->data_0[int(3)] = packed_float4(record_1.V_0[int(3)]) ;

#line 531
            _S901->compX_0 = packed_float4(record_1.compX_0) ;

#line 531
            _S901->compU_0 = packed_float4(record_1.compU_0) ;

#line 531
            _S901->integration_0 = packed_float4(record_1.integration_0) ;

#line 531
            _S901->source_0 = packed_float4(record_1.source_0) ;

#line 531
            _S901->volume_0 = packed_float4(record_1.volume_0) ;

#line 531
            _S901->control_0 = packed_uint4(record_1.control_0) ;

#line 531
            _S901->identity_0 = packed_uint4(record_1.identity_0) ;

#line 529
        }

#line 534
        return;
    }
    bool isWarp_5 = metricId_0 == 2U;
    bool isWormhole_5 = metricId_0 == 3U;

#line 542
    float disc_0 = M_20 * M_20 - a_29 * a_29 - Q_13 * Q_13;

#line 542
    float horizonRadius_0;
    if(disc_0 > 0.0)
    {

#line 543
        horizonRadius_0 = M_20 + sqrt(disc_0);

#line 543
    }
    else
    {

#line 543
        horizonRadius_0 = M_20;

#line 543
    }
    if(isWormhole_5)
    {

#line 544
        horizonRadius_0 = 0.5 * throatB0_5;

#line 544
    }
    float horizonRadius_1 = horizonRadius_0 * captureFactor_0;

    thread GeodesicStateCart_0 state_15;



    thread float3 volumeEmission_0;
    thread float volumeOpticalDepth_0;

    thread array<float4, int(4)> X_1;

#line 554
    thread array<float4, int(4)> V_2;
    Vec4Cart_0 compX_1 = MakeVec4Cart_0(0.0, 0.0, 0.0, 0.0);

#line 555
    float pupilUp_1;

#line 555
    float pupilRight_1;

#line 555
    int column_0;

#line 555
    float3 screenDirection_1;

#line 555
    uint attempts_1;

#line 555
    uint reason_1;

#line 555
    Vec4Cart_0 pos0_0;



    if(_S882)
    {


        float imageX_0 = (2.0 * (float(ix_0) + sampleU_1) / float(imageWidth_0) - 1.0) * aspect_0;
        float imageY_0 = 1.0 - 2.0 * (float(iy_0) + sampleV_1) / float(imageHeight_0);
        float tanHalfFov_1 = tan(fov_0 * 0.5);



        bool _S902 = lensModel_0 == int(1);

#line 569
        if(_S902)
        {

            ThinLensProjectionSample_0 pupil_0 = ProjectThinLensSample_0(imageX_0, imageY_0, tanHalfFov_1, focalLength_1, aperture_1, focusDistance_1, (&kernelContext_6)->params_0[int(66)], (&kernelContext_6)->params_0[int(67)]);

#line 572
            pupilUp_1 = pupil_0.pupilUp_0;

#line 572
            pupilRight_1 = pupil_0.pupilRight_0;

#line 572
            screenDirection_1 = pupil_0.direction_1;

#line 569
        }
        else
        {

#line 579
            float3 _S903 = normalize(float3(1.0, imageY_0 * tanHalfFov_1, imageX_0 * tanHalfFov_1));

#line 579
            pupilUp_1 = 0.0;

#line 579
            pupilRight_1 = 0.0;

#line 579
            screenDirection_1 = _S903;

#line 569
        }

#line 585
        thread array<array<float, int(4)>, int(4)> g_3;
        thread array<array<float, int(4)>, int(4)> g_inv_3;
        Vec4Cart_0 pos0_1 = MakeVec4Cart_0(0.0, _S864.x, _S864.y, _S864.z);
        if(isWarp_5)
        {

#line 588
            thread Vec4Cart_0 _S904 = pos0_1;

#line 588
            GetWarpDriveMetric_0(&_S904, warpVs_5, warpSigma_5, warpR_5, 0.0, 0.0, 0.0, &g_3, &g_inv_3);

#line 588
        }
        else
        {

#line 590
            if(isWormhole_5)
            {

#line 590
                thread Vec4Cart_0 _S905 = pos0_1;

#line 590
                GetMorrisThorneCartesianMetric_0(&_S905, throatB0_5, 0.0, wormShape_8, &g_3, &g_inv_3);

#line 590
            }
            else
            {

#line 590
                thread Vec4Cart_0 _S906 = pos0_1;

#line 590
                GetKerrSchildFamilyMetric_0(&_S906, M_20, a_29, Q_13, Lambda_13, &g_3, &g_inv_3);

#line 590
            }

#line 588
        }

#line 588
        thread array<array<float, int(4)>, int(4)> _S907 = g_3;

#line 588
        thread array<array<float, int(4)>, int(4)> _S908 = g_inv_3;

#line 588
        CameraObserverFrame_0 _S909 = BuildCameraObserverFrame_0(&_S907, &_S908, _S865, _S867, _S866, _S872);

#line 588
        CameraObserverFrame_0 cameraFrame_0;

#line 597
        if(_S902)
        {

#line 597
            thread Vec4Cart_0 _S910 = _S909.up_0;

#line 597
            Vec4Cart_0 _S911 = ScaleVec4Cart_0(&_S910, pupilUp_1);

#line 597
            thread Vec4Cart_0 _S912 = pos0_1;

#line 597
            thread Vec4Cart_0 _S913 = _S911;

#line 597
            Vec4Cart_0 _S914 = AddVec4Cart_0(&_S912, &_S913);

#line 597
            thread Vec4Cart_0 _S915 = _S909.right_0;

#line 597
            Vec4Cart_0 _S916 = ScaleVec4Cart_0(&_S915, pupilRight_1);

#line 597
            thread Vec4Cart_0 _S917 = _S914;

#line 597
            thread Vec4Cart_0 _S918 = _S916;

#line 597
            Vec4Cart_0 _S919 = AddVec4Cart_0(&_S917, &_S918);

#line 604
            if(isWarp_5)
            {

#line 604
                thread Vec4Cart_0 _S920 = _S919;

#line 604
                GetWarpDriveMetric_0(&_S920, warpVs_5, warpSigma_5, warpR_5, 0.0, 0.0, 0.0, &g_3, &g_inv_3);

#line 604
            }
            else
            {

#line 606
                if(isWormhole_5)
                {

#line 606
                    thread Vec4Cart_0 _S921 = _S919;

#line 606
                    GetMorrisThorneCartesianMetric_0(&_S921, throatB0_5, 0.0, wormShape_8, &g_3, &g_inv_3);

#line 606
                }
                else
                {

#line 606
                    thread Vec4Cart_0 _S922 = _S919;

#line 606
                    GetKerrSchildFamilyMetric_0(&_S922, M_20, a_29, Q_13, Lambda_13, &g_3, &g_inv_3);

#line 606
                }

#line 604
            }

#line 604
            thread array<array<float, int(4)>, int(4)> _S923 = g_3;

#line 604
            thread array<array<float, int(4)>, int(4)> _S924 = g_inv_3;

#line 604
            CameraObserverFrame_0 _S925 = BuildCameraObserverFrame_0(&_S923, &_S924, _S865, _S867, _S866, _S872);

#line 604
            pos0_0 = _S919;

#line 604
            cameraFrame_0 = _S925;

#line 597
        }
        else
        {

#line 597
            pos0_0 = pos0_1;

#line 597
            cameraFrame_0 = _S909;

#line 597
        }

#line 614
        (&state_15)->x_5 = pos0_0;

#line 614
        thread CameraObserverFrame_0 _S926 = cameraFrame_0;

#line 614
        Vec4Cart_0 _S927 = PastDirectedCameraRay_0(&_S926, screenDirection_1);
        (&state_15)->u_1 = _S927;

        float3 color_2 = float3(0.0, 0.0, 0.0);

        volumeOpticalDepth_0 = 0.0;
        volumeEmission_0 = float3(0.0) ;

#line 620
        thread CameraObserverFrame_0 _S928 = cameraFrame_0;

#line 625
        thread Vec4Cart_0 launchScreenA_0;
        thread Vec4Cart_0 launchScreenB_0;

#line 626
        SachsScreenBasis_0(&_S928, screenDirection_1, &launchScreenA_0, &launchScreenB_0);

        float4 _S929 = float4(0.0) ;

#line 628
        X_1[int(0)] = _S929;

#line 628
        X_1[int(1)] = _S929;

#line 628
        thread Vec4Cart_0 _S930 = launchScreenA_0;

#line 628
        float4 _S931 = TracePack_0(&_S930);

#line 628
        float4 _S932 = float4(0.00100000004749745) ;
        float4 _S933 = _S931 * _S932;

#line 629
        X_1[int(2)] = _S933;

#line 629
        thread Vec4Cart_0 _S934 = launchScreenB_0;

#line 629
        float4 _S935 = TracePack_0(&_S934);
        float4 _S936 = _S935 * _S932;

#line 630
        X_1[int(3)] = _S936;
        V_2[int(0)] = _S933;

#line 631
        V_2[int(1)] = _S936;
        V_2[int(2)] = _S929;

#line 632
        V_2[int(3)] = _S929;

#line 645
        if(maxSteps_0 <= int(0))
        {

#line 645
            reason_1 = 4U;

#line 645
        }
        else
        {

#line 645
            reason_1 = 0U;

#line 645
        }
        if(isWormhole_5)
        {

#line 646
            escaped_0 = twoSheetWormhole_0;

#line 646
        }
        else
        {

#line 646
            escaped_0 = false;

#line 646
        }

#line 646
        if(escaped_0)
        {

#line 646
            thread GeodesicStateCart_0 _S937 = state_15;

#line 646
            EllisTopologyContinuation_0 _S938 = InitEllisTwoSheet_0(&_S937, maxSteps_0, stepScale_3, minStep_3, maxStep_3, escapeRadius_4, throatB0_5, wormShape_8);



            if((_S938.failure_0) != 0U)
            {

#line 650
                reason_1 = 7U;

#line 650
            }

#line 646
        }

#line 652
        if(sampleIndex_0 < 0.5)
        {

#line 653
            uint _S939 = slot_0 * 4U;

#line 653
            *((&kernelContext_6)->radiance_0+_S939) = 0.0;

#line 653
            *((&kernelContext_6)->radiance_0+(_S939 + 1U)) = 0.0;
            *((&kernelContext_6)->radiance_0+(_S939 + 2U)) = 0.0;

#line 654
            *((&kernelContext_6)->radiance_0+(_S939 + 3U)) = 0.0;

#line 652
        }

#line 652
        thread Vec4Cart_0 _S940 = (&state_15)->x_5;

#line 652
        float4 _S941 = TracePack_0(&_S940);

#line 657
        (&record_1)->x_1 = _S941;

#line 657
        thread Vec4Cart_0 _S942 = (&state_15)->u_1;

#line 657
        float4 _S943 = TracePack_0(&_S942);

#line 657
        (&record_1)->u_0 = _S943;

#line 657
        column_0 = int(0);
        for(;;)
        {

#line 658
            if(column_0 < int(4))
            {
            }
            else
            {

#line 658
                break;
            }

#line 659
            (&record_1)->X_0[column_0] = X_1[column_0];

#line 659
            (&record_1)->V_0[column_0] = V_2[column_0];

#line 658
            column_0 = column_0 + int(1);

#line 658
        }

#line 658
        thread Vec4Cart_0 _S944 = compX_1;

#line 658
        float4 _S945 = TracePack_0(&_S944);


        (&record_1)->compX_0 = _S945;

#line 661
        (&record_1)->compU_0 = _S945;
        (&record_1)->integration_0 = float4(maxStep_3, 0.0, 0.0, 0.00100000004749745);
        (&record_1)->source_0 = float4(color_2, 1.0);
        (&record_1)->volume_0 = float4(volumeEmission_0, volumeOpticalDepth_0);
        if(reason_1 == 0U)
        {

#line 665
            attempts_1 = 1U;

#line 665
        }
        else
        {

#line 665
            attempts_1 = 2U;

#line 665
        }

#line 665
        (&record_1)->control_0 = uint4(attempts_1, reason_1, 0U, 0U);

        (&record_1)->identity_0 = identity_1;

#line 667
        thread TraceContinuation_0 _S946 = record_1;

#line 667
        bool _S947 = TraceStateFinite_0(&_S946);
        if(!_S947)
        {

#line 669
            (&record_1)->control_0.x = 2U;
            (&record_1)->control_0.y = 9U;

#line 668
        }

#line 668
        TraceContinuation_natural_0 device* _S948 = (&kernelContext_6)->continuation_0+slot_0;

#line 668
        _S948->x_1 = packed_float4(record_1.x_1) ;

#line 668
        _S948->u_0 = packed_float4(record_1.u_0) ;

#line 668
        (&_S948->X_0)->data_0[int(0)] = packed_float4(record_1.X_0[int(0)]) ;

#line 668
        (&_S948->X_0)->data_0[int(1)] = packed_float4(record_1.X_0[int(1)]) ;

#line 668
        (&_S948->X_0)->data_0[int(2)] = packed_float4(record_1.X_0[int(2)]) ;

#line 668
        (&_S948->X_0)->data_0[int(3)] = packed_float4(record_1.X_0[int(3)]) ;

#line 668
        (&_S948->V_0)->data_0[int(0)] = packed_float4(record_1.V_0[int(0)]) ;

#line 668
        (&_S948->V_0)->data_0[int(1)] = packed_float4(record_1.V_0[int(1)]) ;

#line 668
        (&_S948->V_0)->data_0[int(2)] = packed_float4(record_1.V_0[int(2)]) ;

#line 668
        (&_S948->V_0)->data_0[int(3)] = packed_float4(record_1.V_0[int(3)]) ;

#line 668
        _S948->compX_0 = packed_float4(record_1.compX_0) ;

#line 668
        _S948->compU_0 = packed_float4(record_1.compU_0) ;

#line 668
        _S948->integration_0 = packed_float4(record_1.integration_0) ;

#line 668
        _S948->source_0 = packed_float4(record_1.source_0) ;

#line 668
        _S948->volume_0 = packed_float4(record_1.volume_0) ;

#line 668
        _S948->control_0 = packed_uint4(record_1.control_0) ;

#line 668
        _S948->identity_0 = packed_uint4(record_1.identity_0) ;

#line 673
        return;
    }

#line 673
    thread TraceContinuation_0 _S949 = record_1;

#line 673
    bool _S950 = TraceStateFinite_0(&_S949);


    if(!_S950)
    {

#line 677
        (&record_1)->control_0.x = 2U;

#line 677
        (&record_1)->control_0.y = 9U;

#line 677
        TraceContinuation_natural_0 device* _S951 = (&kernelContext_6)->continuation_0+slot_0;

#line 677
        _S951->x_1 = packed_float4(record_1.x_1) ;

#line 677
        _S951->u_0 = packed_float4(record_1.u_0) ;

#line 677
        (&_S951->X_0)->data_0[int(0)] = packed_float4(record_1.X_0[int(0)]) ;

#line 677
        (&_S951->X_0)->data_0[int(1)] = packed_float4(record_1.X_0[int(1)]) ;

#line 677
        (&_S951->X_0)->data_0[int(2)] = packed_float4(record_1.X_0[int(2)]) ;

#line 677
        (&_S951->X_0)->data_0[int(3)] = packed_float4(record_1.X_0[int(3)]) ;

#line 677
        (&_S951->V_0)->data_0[int(0)] = packed_float4(record_1.V_0[int(0)]) ;

#line 677
        (&_S951->V_0)->data_0[int(1)] = packed_float4(record_1.V_0[int(1)]) ;

#line 677
        (&_S951->V_0)->data_0[int(2)] = packed_float4(record_1.V_0[int(2)]) ;

#line 677
        (&_S951->V_0)->data_0[int(3)] = packed_float4(record_1.V_0[int(3)]) ;

#line 677
        _S951->compX_0 = packed_float4(record_1.compX_0) ;

#line 677
        _S951->compU_0 = packed_float4(record_1.compU_0) ;

#line 677
        _S951->integration_0 = packed_float4(record_1.integration_0) ;

#line 677
        _S951->source_0 = packed_float4(record_1.source_0) ;

#line 677
        _S951->volume_0 = packed_float4(record_1.volume_0) ;

#line 677
        _S951->control_0 = packed_uint4(record_1.control_0) ;

#line 677
        _S951->identity_0 = packed_uint4(record_1.identity_0) ;
        return;
    }
    (&state_15)->x_5 = TraceUnpack_0((&record_1)->x_1);

#line 680
    (&state_15)->u_1 = TraceUnpack_0((&record_1)->u_0);

#line 680
    column_0 = int(0);
    for(;;)
    {

#line 681
        if(column_0 < int(4))
        {
        }
        else
        {

#line 681
            break;
        }

#line 682
        X_1[column_0] = (&record_1)->X_0[column_0];

#line 682
        V_2[column_0] = (&record_1)->V_0[column_0];

#line 681
        column_0 = column_0 + int(1);

#line 681
    }


    Vec4Cart_0 compX_2 = TraceUnpack_0((&record_1)->compX_0);

#line 684
    Vec4Cart_0 compU_1 = TraceUnpack_0((&record_1)->compU_0);
    float adaptiveStep_2 = (&record_1)->integration_0.x;

#line 685
    float affineLength_8 = (&record_1)->integration_0.y;
    float lastStep_2 = (&record_1)->integration_0.z;
    float3 color_3 = (&record_1)->source_0.xyz;

#line 687
    float transmittance_0 = (&record_1)->source_0.w;
    volumeEmission_0 = (&record_1)->volume_0.xyz;

#line 688
    volumeOpticalDepth_0 = (&record_1)->volume_0.w;
    uint attempts_2 = (&record_1)->control_0.z;

#line 689
    uint accepted_3 = (&record_1)->control_0.w;
    uint reason_2 = (&record_1)->control_0.y;

#line 690
    bool terminalEscapeSurface_0;

#line 690
    float lastStep_3;

#line 690
    float affineLength_9;

#line 690
    float exoticScale_0;

#line 690
    int step_8;

#line 690
    float3 color_4;

    if(_S899)
    {

#line 693
        if(maxSteps_0 <= int(0))
        {

#line 693
            escaped_0 = true;

#line 693
        }
        else
        {

#line 693
            escaped_0 = attempts_2 >= uint(maxSteps_0);

#line 693
        }

#line 693
        if(escaped_0)
        {

#line 694
            (&record_1)->control_0.x = 2U;

#line 694
            (&record_1)->control_0.y = 4U;

#line 694
            TraceContinuation_natural_0 device* _S952 = (&kernelContext_6)->continuation_0+slot_0;

#line 694
            _S952->x_1 = packed_float4(record_1.x_1) ;

#line 694
            _S952->u_0 = packed_float4(record_1.u_0) ;

#line 694
            (&_S952->X_0)->data_0[int(0)] = packed_float4(record_1.X_0[int(0)]) ;

#line 694
            (&_S952->X_0)->data_0[int(1)] = packed_float4(record_1.X_0[int(1)]) ;

#line 694
            (&_S952->X_0)->data_0[int(2)] = packed_float4(record_1.X_0[int(2)]) ;

#line 694
            (&_S952->X_0)->data_0[int(3)] = packed_float4(record_1.X_0[int(3)]) ;

#line 694
            (&_S952->V_0)->data_0[int(0)] = packed_float4(record_1.V_0[int(0)]) ;

#line 694
            (&_S952->V_0)->data_0[int(1)] = packed_float4(record_1.V_0[int(1)]) ;

#line 694
            (&_S952->V_0)->data_0[int(2)] = packed_float4(record_1.V_0[int(2)]) ;

#line 694
            (&_S952->V_0)->data_0[int(3)] = packed_float4(record_1.V_0[int(3)]) ;

#line 694
            _S952->compX_0 = packed_float4(record_1.compX_0) ;

#line 694
            _S952->compU_0 = packed_float4(record_1.compU_0) ;

#line 694
            _S952->integration_0 = packed_float4(record_1.integration_0) ;

#line 694
            _S952->source_0 = packed_float4(record_1.source_0) ;

#line 694
            _S952->volume_0 = packed_float4(record_1.volume_0) ;

#line 694
            _S952->control_0 = packed_uint4(record_1.control_0) ;

#line 694
            _S952->identity_0 = packed_uint4(record_1.identity_0) ;
            return;
        }

        if(isWormhole_5)
        {

#line 698
            escaped_0 = twoSheetWormhole_0;

#line 698
        }
        else
        {

#line 698
            escaped_0 = false;

#line 698
        }

#line 698
        float transmittance_1;

#line 698
        uint accepted_4;

#line 698
        Vec4Cart_0 compU_2;

#line 698
        if(escaped_0)
        {

#line 699
            Vec4Cart_0 _S953 = { 0.0, 0.0, 0.0, 0.0 };

#line 699
            GeodesicStateCart_0 _S954 = { _S953, _S953 };

#line 699
            EllisTopologyTraceResult_0 _S955 = { _S954, int(0), int(0), 0.0 };

#line 699
            thread EllisTopologyContinuation_0 ellis_0;

#line 699
            (&ellis_0)->trace_0 = _S955;

#line 699
            (&ellis_0)->adaptiveStep_0 = 0.0;

#line 699
            (&ellis_0)->compensationX_0 = _S953;

#line 699
            (&ellis_0)->compensationU_0 = _S953;

#line 699
            (&ellis_0)->accepted_0 = 0U;

#line 699
            (&ellis_0)->lastStep_0 = 0.0;

#line 699
            (&ellis_0)->failure_0 = 0U;
            (&(&ellis_0)->trace_0)->state_0 = state_15;
            (&(&ellis_0)->trace_0)->attempts_0 = int(attempts_2);
            (&(&ellis_0)->trace_0)->affineLength_0 = affineLength_8;
            (&ellis_0)->adaptiveStep_0 = adaptiveStep_2;
            (&ellis_0)->compensationX_0 = compX_2;

#line 704
            (&ellis_0)->compensationU_0 = compU_1;
            (&ellis_0)->accepted_0 = accepted_3;

#line 705
            (&ellis_0)->lastStep_0 = lastStep_2;
            AdvanceEllisTwoSheet_0(&ellis_0, maxSteps_0, stepScale_3, minStep_3, maxStep_3, escapeRadius_4, throatB0_5, wormShape_8);

            state_15 = (&(&ellis_0)->trace_0)->state_0;

#line 708
            uint _S956 = uint((&(&ellis_0)->trace_0)->attempts_0);
            float _S957 = (&(&ellis_0)->trace_0)->affineLength_0;

#line 709
            float _S958 = (&ellis_0)->adaptiveStep_0;
            Vec4Cart_0 _S959 = (&ellis_0)->compensationX_0;

#line 710
            Vec4Cart_0 _S960 = (&ellis_0)->compensationU_0;
            uint _S961 = (&ellis_0)->accepted_0;

#line 711
            float _S962 = (&ellis_0)->lastStep_0;
            if(((&(&ellis_0)->trace_0)->termination_0) == int(1))
            {

#line 712
                reason_1 = 3U;

#line 712
            }
            else
            {

#line 713
                if(((&(&ellis_0)->trace_0)->termination_0) == int(2))
                {

#line 713
                    reason_1 = 10U;

#line 713
                }
                else
                {

#line 715
                    if(((&ellis_0)->failure_0) == 3U)
                    {

#line 715
                        reason_1 = 5U;

#line 715
                    }
                    else
                    {

#line 716
                        if(((&ellis_0)->failure_0) == 2U)
                        {

#line 716
                            reason_1 = 4U;

#line 716
                        }
                        else
                        {

#line 717
                            if(((&ellis_0)->failure_0) != 0U)
                            {

#line 717
                                reason_1 = 7U;

#line 717
                            }
                            else
                            {

#line 717
                                reason_1 = reason_2;

#line 717
                            }

#line 716
                        }

#line 715
                    }

#line 713
                }

#line 712
            }

#line 712
            pos0_0 = _S959;

#line 712
            compU_2 = _S960;

#line 712
            pupilUp_1 = _S958;

#line 712
            lastStep_3 = _S962;

#line 712
            color_4 = color_3;

#line 712
            transmittance_1 = transmittance_0;

#line 712
            attempts_1 = _S956;

#line 712
            affineLength_9 = _S957;

#line 712
            accepted_4 = _S961;

#line 698
        }
        else
        {

#line 698
            pupilUp_1 = adaptiveStep_2;

#line 698
            pupilRight_1 = transmittance_0;

#line 698
            lastStep_3 = lastStep_2;

#line 698
            screenDirection_1 = color_3;

#line 698
            step_8 = int(0);

#line 698
            attempts_1 = attempts_2;

#line 698
            affineLength_9 = affineLength_8;

#line 698
            accepted_4 = accepted_3;

#line 719
            for(;;)
            {

#line 719
                if(step_8 < int(1))
                {
                }
                else
                {

#line 719
                    reason_1 = reason_2;

#line 719
                    color_4 = screenDirection_1;

#line 719
                    transmittance_1 = pupilRight_1;

#line 719
                    break;
                }

#line 720
                float px_1 = (&(&state_15)->x_5)->x_2;
                float py_1 = (&(&state_15)->x_5)->y_0;
                float pz_0 = (&(&state_15)->x_5)->z_0;
                float rBig_0 = sqrt((&(&state_15)->x_5)->x_2 * (&(&state_15)->x_5)->x_2 + (&(&state_15)->x_5)->y_0 * (&(&state_15)->x_5)->y_0 + (&(&state_15)->x_5)->z_0 * (&(&state_15)->x_5)->z_0);

                if(isWormhole_5)
                {

#line 726
                    if(rBig_0 <= horizonRadius_1)
                    {

#line 726
                        reason_1 = 1U;

#line 726
                        color_4 = screenDirection_1;

#line 726
                        transmittance_1 = pupilRight_1;

                        break;
                    }

#line 725
                }
                else
                {


                    if(!isWarp_5)
                    {
                        if((ComputeKerrRadiusCart_0(px_1, py_1, pz_0, a_29)) <= horizonRadius_1)
                        {

#line 732
                            reason_1 = 1U;

#line 732
                            color_4 = screenDirection_1;

#line 732
                            transmittance_1 = pupilRight_1;

                            break;
                        }

#line 730
                    }

#line 725
                }

#line 742
                if(isWarp_5)
                {

#line 742
                    exoticScale_0 = max(warpR_5, 1.0 / warpSigma_5);

#line 742
                }
                else
                {

#line 742
                    exoticScale_0 = throatB0_5;

#line 742
                }

                if(isWarp_5)
                {

#line 744
                    escaped_0 = true;

#line 744
                }
                else
                {

#line 744
                    escaped_0 = isWormhole_5;

#line 744
                }

#line 744
                float rStep_0;
                if(escaped_0)
                {

#line 745
                    rStep_0 = max(rBig_0, exoticScale_0);

#line 745
                }
                else
                {

#line 745
                    rStep_0 = ComputeKerrRadiusCart_0(px_1, py_1, pz_0, a_29);

#line 745
                }


                float _S963 = min(clamp(stepScale_3 * rStep_0, minStep_3, maxStep_3), pupilUp_1);

                uint attempts_3 = attempts_1 + 1U;



                float _S964 = AdaptiveIntegrationTolerance_0();

#line 754
                float _S965 = AdaptiveNullTolerance_0();

#line 754
                thread GeodesicStateCart_0 _S966 = state_15;

#line 754
                AdaptiveGeodesicStep_0 _S967 = AttemptAdaptiveGeodesicStepAtFloor_0(&_S966, _S963, minStep_3, isWarp_5, isWormhole_5, M_20, a_29, Q_13, Lambda_13, warpVs_5, warpSigma_5, warpR_5, throatB0_5, wormShape_8, _S964, _S965);
                if(!_S967.accepted_1)
                {

#line 756
                    if(AdaptiveAttemptAlreadyAtMinimumStep_0(_S963, minStep_3))
                    {

#line 756
                        reason_1 = 5U;

#line 756
                        color_4 = screenDirection_1;

#line 756
                        transmittance_1 = pupilRight_1;

#line 756
                        attempts_1 = attempts_3;

#line 756
                        break;
                    }

#line 756
                    pupilUp_1 = max(minStep_3, _S963 * _S967.stepFactor_0);

#line 756
                    transmittance_1 = pupilRight_1;

#line 756
                    color_4 = screenDirection_1;

#line 719
                    int _S968 = step_8 + int(1);

#line 719
                    pupilRight_1 = transmittance_1;

#line 719
                    screenDirection_1 = color_4;

#line 719
                    step_8 = _S968;

#line 719
                    attempts_1 = attempts_3;

#line 719
                    continue;
                }

#line 760
                float _S969 = clamp(_S963 * _S967.stepFactor_0, minStep_3, maxStep_3);

                thread GeodesicStateCart_0 next_1 = _S967.candidate_0;

#line 762
                thread Vec4Cart_0 _S970 = (&state_15)->x_5;

#line 762
                thread Vec4Cart_0 _S971 = (&state_15)->u_1;

#line 762
                thread Vec4Cart_0 _S972 = (&next_1)->x_5;

#line 762
                thread Vec4Cart_0 _S973 = (&next_1)->u_1;

#line 776
                thread float escapeFraction_0;
                thread Vec4Cart_0 escapePosition_0;
                thread Vec4Cart_0 escapeTangent_0;

#line 778
                bool _S974 = FindSphericalBoundaryEvent_0(&_S970, &_S971, &_S972, &_S973, _S963, escapeRadius_4, int(1), &escapeFraction_0, &escapePosition_0, &escapeTangent_0);

#line 778
                float dt_6;
                if(_S974)
                {


                    (&next_1)->x_5 = escapePosition_0;
                    (&next_1)->u_1 = escapeTangent_0;
                    float dt_7 = _S963 * escapeFraction_0;

#line 785
                    terminalEscapeSurface_0 = true;

#line 785
                    dt_6 = dt_7;

#line 779
                }
                else
                {

#line 790
                    float nextRadialNumerator_0 = (&(&next_1)->x_5)->x_2 * (&(&next_1)->u_1)->x_2 + (&(&next_1)->x_5)->y_0 * (&(&next_1)->u_1)->y_0 + (&(&next_1)->x_5)->z_0 * (&(&next_1)->u_1)->z_0;

#line 795
                    if(((&(&next_1)->x_5)->x_2 * (&(&next_1)->x_5)->x_2 + (&(&next_1)->x_5)->y_0 * (&(&next_1)->x_5)->y_0 + (&(&next_1)->x_5)->z_0 * (&(&next_1)->x_5)->z_0) > (escapeRadius_4 * escapeRadius_4))
                    {

#line 795
                        terminalEscapeSurface_0 = nextRadialNumerator_0 > 0.0;

#line 795
                    }
                    else
                    {

#line 795
                        terminalEscapeSurface_0 = false;

#line 795
                    }

#line 795
                    if(terminalEscapeSurface_0)
                    {
                        (&record_1)->control_0.x = 2U;
                        (&record_1)->control_0.y = 6U;
                        (&record_1)->control_0.z = attempts_3;

#line 799
                        TraceContinuation_natural_0 device* _S975 = (&kernelContext_6)->continuation_0+slot_0;

#line 799
                        _S975->x_1 = packed_float4(record_1.x_1) ;

#line 799
                        _S975->u_0 = packed_float4(record_1.u_0) ;

#line 799
                        (&_S975->X_0)->data_0[int(0)] = packed_float4(record_1.X_0[int(0)]) ;

#line 799
                        (&_S975->X_0)->data_0[int(1)] = packed_float4(record_1.X_0[int(1)]) ;

#line 799
                        (&_S975->X_0)->data_0[int(2)] = packed_float4(record_1.X_0[int(2)]) ;

#line 799
                        (&_S975->X_0)->data_0[int(3)] = packed_float4(record_1.X_0[int(3)]) ;

#line 799
                        (&_S975->V_0)->data_0[int(0)] = packed_float4(record_1.V_0[int(0)]) ;

#line 799
                        (&_S975->V_0)->data_0[int(1)] = packed_float4(record_1.V_0[int(1)]) ;

#line 799
                        (&_S975->V_0)->data_0[int(2)] = packed_float4(record_1.V_0[int(2)]) ;

#line 799
                        (&_S975->V_0)->data_0[int(3)] = packed_float4(record_1.V_0[int(3)]) ;

#line 799
                        _S975->compX_0 = packed_float4(record_1.compX_0) ;

#line 799
                        _S975->compU_0 = packed_float4(record_1.compU_0) ;

#line 799
                        _S975->integration_0 = packed_float4(record_1.integration_0) ;

#line 799
                        _S975->source_0 = packed_float4(record_1.source_0) ;

#line 799
                        _S975->volume_0 = packed_float4(record_1.volume_0) ;

#line 799
                        _S975->control_0 = packed_uint4(record_1.control_0) ;

#line 799
                        _S975->identity_0 = packed_uint4(record_1.identity_0) ;

                        return;
                    }

#line 801
                    terminalEscapeSurface_0 = false;

#line 801
                    dt_6 = _S963;

#line 779
                }

#line 779
                bool terminalCaptureSurface_0;

#line 779
                bool terminalEscapeSurface_1;

#line 779
                float dt_8;

#line 810
                if(isWormhole_5)
                {

#line 810
                    thread Vec4Cart_0 _S976 = (&state_15)->x_5;

#line 810
                    thread Vec4Cart_0 _S977 = (&state_15)->u_1;

#line 810
                    thread Vec4Cart_0 _S978 = (&next_1)->x_5;

#line 810
                    thread Vec4Cart_0 _S979 = (&next_1)->u_1;
                    thread float captureFraction_0;
                    thread Vec4Cart_0 capturePosition_0;
                    thread Vec4Cart_0 captureTangent_0;

#line 813
                    bool _S980 = FindSphericalCaptureEvent_0(&_S976, &_S977, &_S978, &_S979, dt_6, horizonRadius_1, &captureFraction_0, &capturePosition_0, &captureTangent_0);
                    if(_S980)
                    {

                        (&next_1)->x_5 = capturePosition_0;
                        (&next_1)->u_1 = captureTangent_0;
                        float dt_9 = dt_6 * captureFraction_0;

#line 819
                        terminalCaptureSurface_0 = false;

#line 819
                        dt_8 = dt_9;

#line 814
                    }
                    else
                    {

#line 814
                        terminalCaptureSurface_0 = terminalEscapeSurface_0;

#line 814
                        dt_8 = dt_6;

#line 814
                    }

#line 765
                    bool _S981 = terminalCaptureSurface_0;

#line 765
                    terminalCaptureSurface_0 = _S980;

#line 765
                    terminalEscapeSurface_1 = _S981;

#line 810
                }
                else
                {

#line 810
                    terminalCaptureSurface_0 = false;

#line 810
                    terminalEscapeSurface_1 = terminalEscapeSurface_0;

#line 810
                    dt_8 = dt_6;

#line 810
                }

#line 810
                bool _S982;

#line 826
                if(_S875)
                {

#line 826
                    _S982 = !isWarp_5;

#line 826
                }
                else
                {

#line 826
                    _S982 = false;

#line 826
                }

#line 826
                bool _S983;

#line 826
                if(_S982)
                {

#line 826
                    _S983 = !isWormhole_5;

#line 826
                }
                else
                {

#line 826
                    _S983 = false;

#line 826
                }

#line 826
                float cr_0;

#line 826
                if(_S983)
                {

#line 827
                    if(M_20 > 9.99999997475242708e-07)
                    {

#line 827
                        cr_0 = a_29 / M_20;

#line 827
                    }
                    else
                    {

#line 827
                        cr_0 = a_29;

#line 827
                    }

#line 827
                    thread Vec4Cart_0 _S984 = (&state_15)->x_5;

#line 827
                    thread Vec4Cart_0 _S985 = (&next_1)->x_5;

#line 827
                    thread Vec4Cart_0 _S986 = (&state_15)->u_1;

#line 827
                    thread Vec4Cart_0 _S987 = (&next_1)->u_1;

#line 827
                    AccumulateVolumeSegment_0(&_S984, &_S985, &_S986, &_S987, dt_8, _S879, _S875, _S880, _S874, _S869, _S870, _S871, _S873, M_20, a_29, cr_0, _S876, _S877, _S878, 1.0, &volumeOpticalDepth_0, &volumeEmission_0);

#line 826
                }

#line 826
                bool _S988;

#line 841
                if(_S868)
                {

#line 841
                    _S988 = !_S875;

#line 841
                }
                else
                {

#line 841
                    _S988 = false;

#line 841
                }

#line 841
                bool _S989;

#line 841
                if(_S988)
                {

#line 841
                    _S989 = !isWarp_5;

#line 841
                }
                else
                {

#line 841
                    _S989 = false;

#line 841
                }

#line 841
                bool _S990;

#line 841
                if(_S989)
                {

#line 841
                    _S990 = !isWormhole_5;

#line 841
                }
                else
                {

#line 841
                    _S990 = false;

#line 841
                }

#line 841
                bool _S991;

#line 841
                if(_S990)
                {

#line 841
                    _S991 = dt_8 > 0.0;

#line 841
                }
                else
                {

#line 841
                    _S991 = false;

#line 841
                }

#line 841
                bool _S992;
                if(_S991)
                {

#line 842
                    _S992 = pupilRight_1 > 0.00999999977648258;

#line 842
                }
                else
                {

#line 842
                    _S992 = false;

#line 842
                }

#line 842
                bool _S993;

#line 842
                bool _S994;

#line 842
                bool coplanar_0;

#line 842
                bool terminalDiskHit_0;

#line 842
                float dt_10;

#line 842
                int pair_0;

#line 841
                if(_S992)
                {
                    float dzStart_0 = dt_8 * (&(&state_15)->u_1)->z_0;
                    float dzEnd_0 = dt_8 * (&(&next_1)->u_1)->z_0;

                    float coefficientA_1 = 2.0 * (&(&state_15)->x_5)->z_0 - 2.0 * (&(&next_1)->x_5)->z_0 + dzStart_0 + dzEnd_0;

                    float coefficientB_1 = -3.0 * (&(&state_15)->x_5)->z_0 + 3.0 * (&(&next_1)->x_5)->z_0 - 2.0 * dzStart_0 - dzEnd_0;

                    float coefficientD_1 = (&(&state_15)->x_5)->z_0;

                    float _S995 = abs(coefficientA_1);

#line 852
                    float _S996 = abs(coefficientB_1);
                    float _S997 = abs(dzStart_0);

#line 853
                    float _S998 = abs((&(&state_15)->x_5)->z_0);

#line 858
                    float rootTolerance_0 = 0.0000079999999798 * max(1.0, max(max(_S995, _S996), max(_S997, _S998)));


                    thread array<float, int(4)> rootFractions_0;
                    rootFractions_0[int(0)] = 0.0;
                    rootFractions_0[int(1)] = 0.0;
                    rootFractions_0[int(2)] = 0.0;
                    rootFractions_0[int(3)] = 0.0;

                    if(_S995 <= rootTolerance_0)
                    {

#line 867
                        _S993 = _S996 <= rootTolerance_0;

#line 867
                    }
                    else
                    {

#line 867
                        _S993 = false;

#line 867
                    }
                    if(_S993)
                    {

#line 868
                        _S994 = _S997 <= rootTolerance_0;

#line 868
                    }
                    else
                    {

#line 868
                        _S994 = false;

#line 868
                    }
                    if(_S994)
                    {

#line 869
                        coplanar_0 = _S998 <= rootTolerance_0;

#line 869
                    }
                    else
                    {

#line 869
                        coplanar_0 = false;

#line 869
                    }

#line 869
                    bool foundDiskCrossing_0;

#line 869
                    bool _S999;

#line 869
                    bool midpointInside_0;

#line 869
                    float frac_0;

#line 869
                    int radialIteration_0;

#line 869
                    int rootCount_3;

                    if(coplanar_0)
                    {

                        float startRadius_0 = ComputeKerrRadiusCart_0((&(&state_15)->x_5)->x_2, (&(&state_15)->x_5)->y_0, 0.0, a_29);

                        if(startRadius_0 >= _S869)
                        {

#line 876
                            foundDiskCrossing_0 = startRadius_0 <= _S870;

#line 876
                        }
                        else
                        {

#line 876
                            foundDiskCrossing_0 = false;

#line 876
                        }

#line 876
                        if(foundDiskCrossing_0)
                        {

#line 877
                            rootFractions_0[int(0)] = 0.0;

#line 877
                            rootCount_3 = int(1);

#line 876
                        }
                        else
                        {

#line 876
                            cr_0 = 0.0;

#line 876
                            pair_0 = int(1);

#line 876
                            rootCount_3 = int(0);



                            for(;;)
                            {

#line 880
                                if(pair_0 <= int(64))
                                {

#line 880
                                    _S999 = rootCount_3 == int(0);

#line 880
                                }
                                else
                                {

#line 880
                                    _S999 = false;

#line 880
                                }

#line 880
                                if(_S999)
                                {
                                }
                                else
                                {

#line 880
                                    break;
                                }
                                float radialUpper_0 = float(pair_0) / 64.0;

#line 882
                                thread Vec4Cart_0 _S1000 = (&state_15)->x_5;

#line 882
                                thread Vec4Cart_0 _S1001 = (&state_15)->u_1;

#line 882
                                thread Vec4Cart_0 _S1002 = (&next_1)->x_5;

#line 882
                                thread Vec4Cart_0 _S1003 = (&next_1)->u_1;

#line 882
                                Vec4Cart_0 _S1004 = AcceptedSegmentPosition_0(&_S1000, &_S1001, &_S1002, &_S1003, dt_8, radialUpper_0);


                                float upperRadius_0 = ComputeKerrRadiusCart_0(_S1004.x_2, _S1004.y_0, 0.0, a_29);

                                if(upperRadius_0 >= _S869)
                                {

#line 887
                                    terminalDiskHit_0 = upperRadius_0 <= _S870;

#line 887
                                }
                                else
                                {

#line 887
                                    terminalDiskHit_0 = false;

#line 887
                                }
                                if(terminalDiskHit_0)
                                {

#line 888
                                    frac_0 = cr_0;

#line 888
                                    transmittance_1 = radialUpper_0;

#line 888
                                    radialIteration_0 = int(0);
                                    for(;;)
                                    {

#line 889
                                        if(radialIteration_0 < int(40))
                                        {
                                        }
                                        else
                                        {

#line 889
                                            break;
                                        }

                                        float radialMidpoint_0 = (frac_0 + transmittance_1) * 0.5;

#line 892
                                        thread Vec4Cart_0 _S1005 = (&state_15)->x_5;

#line 892
                                        thread Vec4Cart_0 _S1006 = (&state_15)->u_1;

#line 892
                                        thread Vec4Cart_0 _S1007 = (&next_1)->x_5;

#line 892
                                        thread Vec4Cart_0 _S1008 = (&next_1)->u_1;

#line 892
                                        Vec4Cart_0 _S1009 = AcceptedSegmentPosition_0(&_S1005, &_S1006, &_S1007, &_S1008, dt_8, radialMidpoint_0);


                                        float midpointRadius_0 = ComputeKerrRadiusCart_0(_S1009.x_2, _S1009.y_0, 0.0, a_29);

                                        if(midpointRadius_0 >= _S869)
                                        {

#line 897
                                            midpointInside_0 = midpointRadius_0 <= _S870;

#line 897
                                        }
                                        else
                                        {

#line 897
                                            midpointInside_0 = false;

#line 897
                                        }

                                        if(midpointInside_0)
                                        {

#line 899
                                            transmittance_1 = radialMidpoint_0;

#line 899
                                        }
                                        else
                                        {

#line 899
                                            frac_0 = radialMidpoint_0;

#line 899
                                        }

#line 889
                                        radialIteration_0 = radialIteration_0 + int(1);

#line 889
                                    }

#line 905
                                    int _S1010 = rootCount_3 + int(1);

#line 905
                                    rootFractions_0[rootCount_3] = transmittance_1;

#line 905
                                    rootCount_3 = _S1010;

#line 888
                                }
                                else
                                {

#line 888
                                    transmittance_1 = radialUpper_0;

#line 888
                                }

#line 881
                                int radialInterval_0 = pair_0 + int(1);

#line 881
                                cr_0 = transmittance_1;

#line 881
                                pair_0 = radialInterval_0;

#line 880
                            }

#line 876
                        }

#line 871
                    }
                    else
                    {

#line 915
                        thread array<float, int(4)> breakpoints_0;
                        breakpoints_0[int(0)] = 0.0;
                        breakpoints_0[int(1)] = 1.0;


                        float _S1011 = coefficientB_1 * coefficientB_1;
                        float _S1012 = 3.0 * coefficientA_1;

#line 921
                        float _S1013 = _S1012 * dzStart_0;

#line 920
                        float derivativeDiscriminant_0 = _S1011 - _S1013;

#line 928
                        float discriminantTolerance_0 = 0.0000079999999798 * max(1.0, max(_S1011, abs(_S1013)));

                        if(_S995 > rootTolerance_0)
                        {

#line 930
                            foundDiskCrossing_0 = derivativeDiscriminant_0 >= (- discriminantTolerance_0);

#line 930
                        }
                        else
                        {

#line 930
                            foundDiskCrossing_0 = false;

#line 930
                        }

#line 930
                        if(foundDiskCrossing_0)
                        {
                            float derivativeRoot_0 = sqrt(max(derivativeDiscriminant_0, 0.0));
                            float _S1014 = - coefficientB_1;

#line 933
                            float firstCritical_0 = (_S1014 - derivativeRoot_0) / _S1012;

                            float secondCritical_0 = (_S1014 + derivativeRoot_0) / _S1012;

                            if(firstCritical_0 > 0.0)
                            {

#line 937
                                _S999 = firstCritical_0 < 1.0;

#line 937
                            }
                            else
                            {

#line 937
                                _S999 = false;

#line 937
                            }

#line 937
                            if(_S999)
                            {

#line 938
                                breakpoints_0[int(2)] = firstCritical_0;

#line 938
                                pair_0 = int(3);

#line 937
                            }
                            else
                            {

#line 937
                                pair_0 = int(2);

#line 937
                            }

                            if(secondCritical_0 > 0.0)
                            {

#line 939
                                terminalDiskHit_0 = secondCritical_0 < 1.0;

#line 939
                            }
                            else
                            {

#line 939
                                terminalDiskHit_0 = false;

#line 939
                            }

#line 939
                            if(terminalDiskHit_0)
                            {

#line 940
                                int _S1015 = pair_0 + int(1);

#line 940
                                breakpoints_0[pair_0] = secondCritical_0;

#line 940
                                radialIteration_0 = _S1015;

#line 939
                            }
                            else
                            {

#line 939
                                radialIteration_0 = pair_0;

#line 939
                            }

#line 939
                            pair_0 = radialIteration_0;

#line 930
                        }
                        else
                        {

#line 941
                            if(_S996 > rootTolerance_0)
                            {

#line 942
                                float critical_5 = - dzStart_0 / (2.0 * coefficientB_1);
                                if(critical_5 > 0.0)
                                {

#line 943
                                    _S999 = critical_5 < 1.0;

#line 943
                                }
                                else
                                {

#line 943
                                    _S999 = false;

#line 943
                                }

#line 943
                                if(_S999)
                                {

#line 944
                                    breakpoints_0[int(2)] = critical_5;

#line 944
                                    pair_0 = int(3);

#line 943
                                }
                                else
                                {

#line 943
                                    pair_0 = int(2);

#line 943
                                }

#line 941
                            }
                            else
                            {

#line 941
                                pair_0 = int(2);

#line 941
                            }

#line 930
                        }

#line 930
                        radialIteration_0 = int(0);

#line 946
                        for(;;)
                        {

#line 946
                            if(radialIteration_0 < pair_0)
                            {
                            }
                            else
                            {

#line 946
                                break;
                            }

#line 947
                            int _S1016 = radialIteration_0 + int(1);

#line 947
                            right_2 = _S1016;

#line 947
                            for(;;)
                            {

#line 947
                                if(right_2 < pair_0)
                                {
                                }
                                else
                                {

#line 947
                                    break;
                                }

#line 948
                                if((breakpoints_0[right_2]) < breakpoints_0[radialIteration_0])
                                {

#line 949
                                    float swapValue_0 = breakpoints_0[radialIteration_0];
                                    breakpoints_0[radialIteration_0] = breakpoints_0[right_2];
                                    breakpoints_0[right_2] = swapValue_0;

#line 948
                                }

#line 947
                                right_2 = right_2 + int(1);

#line 947
                            }

#line 946
                            radialIteration_0 = _S1016;

#line 946
                        }

#line 946
                        right_2 = int(0);

#line 946
                        rootCount_3 = int(0);

#line 955
                        for(;;)
                        {

#line 955
                            int _S1017 = right_2 + int(1);

#line 955
                            if(_S1017 < pair_0)
                            {
                            }
                            else
                            {

#line 955
                                break;
                            }
                            float lower_1 = breakpoints_0[right_2];
                            float upper_1 = breakpoints_0[_S1017];
                            float lowerZ_0 = EvaluateCubic_0(coefficientA_1, coefficientB_1, dzStart_0, coefficientD_1, breakpoints_0[right_2]);

                            float upperZ_0 = EvaluateCubic_0(coefficientA_1, coefficientB_1, dzStart_0, coefficientD_1, breakpoints_0[_S1017]);

                            if((abs(lowerZ_0)) <= rootTolerance_0)
                            {

#line 963
                                _S999 = rootCount_3 < int(4);

#line 963
                            }
                            else
                            {

#line 963
                                _S999 = false;

#line 963
                            }

#line 963
                            if(_S999)
                            {

#line 964
                                if(rootCount_3 == int(0))
                                {

#line 964
                                    terminalDiskHit_0 = true;

#line 964
                                }
                                else
                                {

#line 964
                                    terminalDiskHit_0 = (abs(lower_1 - rootFractions_0[rootCount_3 - int(1)])) > 0.00001999999949476;

#line 964
                                }

#line 963
                            }
                            else
                            {

#line 963
                                terminalDiskHit_0 = false;

#line 963
                            }

#line 963
                            if(terminalDiskHit_0)
                            {


                                int _S1018 = rootCount_3 + int(1);

#line 967
                                rootFractions_0[rootCount_3] = lower_1;

#line 967
                                rootCount_2 = _S1018;

#line 963
                            }
                            else
                            {

#line 963
                                rootCount_2 = rootCount_3;

#line 963
                            }

#line 969
                            if((lowerZ_0 < 0.0) != (upperZ_0 < 0.0))
                            {

#line 969
                                cr_0 = lower_1;

#line 969
                                frac_0 = upper_1;

#line 969
                                transmittance_1 = lowerZ_0;

#line 969
                                int rootIteration_0 = int(0);



                                for(;;)
                                {

#line 973
                                    if(rootIteration_0 < int(40))
                                    {
                                    }
                                    else
                                    {

#line 973
                                        break;
                                    }

#line 974
                                    float midpoint_1 = (cr_0 + frac_0) * 0.5;
                                    float midpointZ_0 = EvaluateCubic_0(coefficientA_1, coefficientB_1, dzStart_0, coefficientD_1, midpoint_1);


                                    if((transmittance_1 < 0.0) == (midpointZ_0 < 0.0))
                                    {

#line 978
                                        cr_0 = midpoint_1;

#line 978
                                        transmittance_1 = midpointZ_0;

#line 978
                                    }
                                    else
                                    {

#line 978
                                        frac_0 = midpoint_1;

#line 978
                                    }

#line 973
                                    rootIteration_0 = rootIteration_0 + int(1);

#line 973
                                }

#line 986
                                float rootFraction_0 = (cr_0 + frac_0) * 0.5;
                                if(rootCount_2 < int(4))
                                {

#line 988
                                    if(rootCount_2 == int(0))
                                    {

#line 988
                                        midpointInside_0 = true;

#line 988
                                    }
                                    else
                                    {

#line 988
                                        midpointInside_0 = (abs(rootFraction_0 - rootFractions_0[rootCount_2 - int(1)])) > 0.00001999999949476;

#line 988
                                    }

#line 987
                                }
                                else
                                {

#line 987
                                    midpointInside_0 = false;

#line 987
                                }

#line 987
                                int rootCount_4;

#line 987
                                if(midpointInside_0)
                                {


                                    int _S1019 = rootCount_2 + int(1);

#line 991
                                    rootFractions_0[rootCount_2] = rootFraction_0;

#line 991
                                    rootCount_4 = _S1019;

#line 987
                                }
                                else
                                {

#line 987
                                    rootCount_4 = rootCount_2;

#line 987
                                }

#line 987
                                rootCount_3 = rootCount_4;

#line 969
                            }
                            else
                            {

#line 969
                                rootCount_3 = rootCount_2;

#line 969
                            }

#line 955
                            right_2 = _S1017;

#line 955
                        }

#line 997
                        if((abs(EvaluateCubic_0(coefficientA_1, coefficientB_1, dzStart_0, coefficientD_1, 1.0))) <= rootTolerance_0)
                        {

#line 997
                            _S999 = rootCount_3 < int(4);

#line 997
                        }
                        else
                        {

#line 997
                            _S999 = false;

#line 997
                        }

#line 997
                        if(_S999)
                        {

#line 998
                            if(rootCount_3 == int(0))
                            {

#line 998
                                terminalDiskHit_0 = true;

#line 998
                            }
                            else
                            {

#line 998
                                terminalDiskHit_0 = (abs(1.0 - rootFractions_0[rootCount_3 - int(1)])) > 0.00001999999949476;

#line 998
                            }

#line 997
                        }
                        else
                        {

#line 997
                            terminalDiskHit_0 = false;

#line 997
                        }

#line 997
                        if(terminalDiskHit_0)
                        {


                            int _S1020 = rootCount_3 + int(1);

#line 1001
                            rootFractions_0[rootCount_3] = 1.0;

#line 1001
                            rootCount_2 = _S1020;

#line 997
                        }
                        else
                        {

#line 997
                            rootCount_2 = rootCount_3;

#line 997
                        }

#line 997
                        rootCount_3 = rootCount_2;

#line 871
                    }

#line 871
                    foundDiskCrossing_0 = false;

#line 871
                    cr_0 = 0.0;

#line 871
                    frac_0 = 0.0;

#line 871
                    pair_0 = int(0);

#line 1010
                    for(;;)
                    {

#line 1010
                        if(pair_0 < rootCount_3)
                        {

#line 1010
                            _S999 = !foundDiskCrossing_0;

#line 1010
                        }
                        else
                        {

#line 1010
                            _S999 = false;

#line 1010
                        }

#line 1010
                        if(_S999)
                        {
                        }
                        else
                        {

#line 1010
                            break;
                        }
                        float candidateFraction_2 = rootFractions_0[pair_0];

#line 1012
                        thread Vec4Cart_0 _S1021 = (&state_15)->x_5;

#line 1012
                        thread Vec4Cart_0 _S1022 = (&state_15)->u_1;

#line 1012
                        thread Vec4Cart_0 _S1023 = (&next_1)->x_5;

#line 1012
                        thread Vec4Cart_0 _S1024 = (&next_1)->u_1;

#line 1012
                        Vec4Cart_0 _S1025 = AcceptedSegmentPosition_0(&_S1021, &_S1022, &_S1023, &_S1024, dt_8, rootFractions_0[pair_0]);


                        float candidateRadius_0 = ComputeKerrRadiusCart_0(_S1025.x_2, _S1025.y_0, 0.0, a_29);

                        if(candidateRadius_0 >= _S869)
                        {

#line 1017
                            terminalDiskHit_0 = candidateRadius_0 <= _S870;

#line 1017
                        }
                        else
                        {

#line 1017
                            terminalDiskHit_0 = false;

#line 1017
                        }

#line 1017
                        if(terminalDiskHit_0)
                        {

#line 1017
                            thread Vec4Cart_0 _S1026 = (&state_15)->x_5;

#line 1017
                            thread Vec4Cart_0 _S1027 = (&state_15)->u_1;

#line 1017
                            thread Vec4Cart_0 _S1028 = (&next_1)->x_5;

#line 1017
                            thread Vec4Cart_0 _S1029 = (&next_1)->u_1;

#line 1017
                            Vec4Cart_0 _S1030 = AcceptedSegmentTangent_0(&_S1026, &_S1027, &_S1028, &_S1029, dt_8, candidateFraction_2);

#line 1017
                            foundDiskCrossing_0 = true;

#line 1017
                            pos0_0 = _S1025;

#line 1017
                            cr_0 = candidateRadius_0;

#line 1017
                            compU_2 = _S1030;

#line 1017
                            frac_0 = candidateFraction_2;

#line 1017
                        }

#line 1010
                        pair_0 = pair_0 + int(1);

#line 1010
                    }

#line 1026
                    if(foundDiskCrossing_0)
                    {

#line 1007
                        Vec4Cart_0 _S1031 = pos0_0;

#line 1007
                        Vec4Cart_0 _S1032 = pos0_0;

#line 1031
                        if(M_20 > 9.99999997475242708e-07)
                        {

#line 1031
                            transmittance_1 = a_29 / M_20;

#line 1031
                        }
                        else
                        {

#line 1031
                            transmittance_1 = a_29;

#line 1031
                        }

                        if(_S873 == int(1))
                        {

#line 1033
                            dt_10 = ShakuraSunyaevTemperature_0(_S871, cr_0, _S869);

#line 1033
                        }
                        else
                        {

#line 1033
                            dt_10 = FullPageThorneTemperature_0(_S871, cr_0, _S869, M_20, transmittance_1);

#line 1033
                        }

#line 1039
                        float cphi_0 = atan2(_S1032.y_0, _S1031.x_2);

#line 1039
                        thread Vec4Cart_0 _S1033 = compU_2;

#line 1039
                        float2 _S1034 = ComputeKerrDiskTransferCart_0(cr_0, cphi_0, &_S1033, 1.0, M_20, a_29);

#line 1039
                        float gFactor_1;


                        if(_S874)
                        {

#line 1042
                            gFactor_1 = _S1034.x;

#line 1042
                        }
                        else
                        {

#line 1042
                            gFactor_1 = _S1034.y;

#line 1042
                        }



                        float3 emit_0 = BlackbodyColor_0(dt_10 * gFactor_1) * float3(ObservedDiskBolometricSource_0(pow(dt_10 / _S871, 4.0), gFactor_1)) ;

#line 1051
                        (&next_1)->x_5 = pos0_0;
                        (&next_1)->u_1 = compU_2;
                        float dt_11 = dt_8 * frac_0;

#line 1053
                        terminalDiskHit_0 = true;

#line 1053
                        transmittance_1 = 0.0;

#line 1053
                        color_4 = emit_0;

#line 1053
                        dt_10 = dt_11;

#line 1026
                    }
                    else
                    {

#line 1026
                        terminalDiskHit_0 = false;

#line 1026
                        transmittance_1 = pupilRight_1;

#line 1026
                        color_4 = screenDirection_1;

#line 1026
                        dt_10 = dt_8;

#line 1026
                    }

#line 841
                }
                else
                {

#line 841
                    terminalDiskHit_0 = false;

#line 841
                    transmittance_1 = pupilRight_1;

#line 841
                    color_4 = screenDirection_1;

#line 841
                    dt_10 = dt_8;

#line 841
                }

#line 1057
                if(beamsEnabled_0)
                {

#line 1057
                    _S993 = !isWarp_5;

#line 1057
                }
                else
                {

#line 1057
                    _S993 = false;

#line 1057
                }

#line 1057
                if(_S993)
                {

#line 1057
                    _S994 = !isWormhole_5;

#line 1057
                }
                else
                {

#line 1057
                    _S994 = false;

#line 1057
                }

#line 1057
                if(_S994)
                {

#line 1057
                    coplanar_0 = dt_10 > 0.0;

#line 1057
                }
                else
                {

#line 1057
                    coplanar_0 = false;

#line 1057
                }

#line 1057
                if(coplanar_0)
                {

#line 1057
                    pair_0 = int(0);
                    for(;;)
                    {

#line 1058
                        if(pair_0 < int(2))
                        {
                        }
                        else
                        {

#line 1058
                            break;
                        }

#line 1059
                        int firstColumn_0 = int(2) * pair_0;



                        int _S1035 = firstColumn_0 + int(1);

#line 1063
                        thread Vec4Cart_0 _S1036 = (&state_15)->x_5;

#line 1063
                        thread Vec4Cart_0 _S1037 = (&state_15)->u_1;

#line 1063
                        thread Vec4Cart_0 _S1038 = (&next_1)->x_5;

#line 1063
                        thread Vec4Cart_0 _S1039 = (&next_1)->u_1;

#line 1063
                        AdvanceJacobiPairAcceptedSegment_0(&_S1036, &_S1037, &_S1038, &_S1039, dt_10, M_20, a_29, Q_13, Lambda_13, &X_1[firstColumn_0], &V_2[firstColumn_0], &X_1[_S1035], &V_2[_S1035]);

#line 1058
                        pair_0 = pair_0 + int(1);

#line 1058
                    }

#line 1057
                }

#line 1067
                float affineLength_10 = affineLength_9 + dt_10;

#line 1067
                uint accepted_5 = accepted_4 + 1U;
                state_15 = next_1;
                if(terminalDiskHit_0)
                {

#line 1069
                    reason_1 = 2U;

#line 1069
                    pupilUp_1 = _S969;

#line 1069
                    lastStep_3 = dt_10;

#line 1069
                    attempts_1 = attempts_3;

#line 1069
                    affineLength_9 = affineLength_10;

#line 1069
                    accepted_4 = accepted_5;

#line 1069
                    break;
                }

#line 1070
                if(terminalCaptureSurface_0)
                {

#line 1070
                    reason_1 = 1U;

#line 1070
                    pupilUp_1 = _S969;

#line 1070
                    lastStep_3 = dt_10;

#line 1070
                    attempts_1 = attempts_3;

#line 1070
                    affineLength_9 = affineLength_10;

#line 1070
                    accepted_4 = accepted_5;

                    break;
                }

#line 1080
                float nextRadius_0 = sqrt((&(&state_15)->x_5)->x_2 * (&(&state_15)->x_5)->x_2 + (&(&state_15)->x_5)->y_0 * (&(&state_15)->x_5)->y_0 + (&(&state_15)->x_5)->z_0 * (&(&state_15)->x_5)->z_0);

                if(isWormhole_5)
                {

#line 1083
                    if(nextRadius_0 <= horizonRadius_1)
                    {

#line 1083
                        reason_1 = 1U;

#line 1083
                        pupilUp_1 = _S969;

#line 1083
                        lastStep_3 = dt_10;

#line 1083
                        attempts_1 = attempts_3;

#line 1083
                        affineLength_9 = affineLength_10;

#line 1083
                        accepted_4 = accepted_5;

                        break;
                    }

#line 1082
                }
                else
                {


                    if(!isWarp_5)
                    {

                        if((ComputeKerrRadiusCart_0((&(&state_15)->x_5)->x_2, (&(&state_15)->x_5)->y_0, (&(&state_15)->x_5)->z_0, a_29)) <= horizonRadius_1)
                        {

#line 1090
                            reason_1 = 1U;

#line 1090
                            pupilUp_1 = _S969;

#line 1090
                            lastStep_3 = dt_10;

#line 1090
                            attempts_1 = attempts_3;

#line 1090
                            affineLength_9 = affineLength_10;

#line 1090
                            accepted_4 = accepted_5;

                            break;
                        }

#line 1087
                    }

#line 1082
                }

#line 1095
                if(terminalEscapeSurface_1)
                {

#line 1095
                    reason_1 = 3U;

#line 1095
                    pupilUp_1 = _S969;

#line 1095
                    lastStep_3 = dt_10;

#line 1095
                    attempts_1 = attempts_3;

#line 1095
                    affineLength_9 = affineLength_10;

#line 1095
                    accepted_4 = accepted_5;

                    break;
                }

#line 1097
                pupilUp_1 = _S969;

#line 1097
                lastStep_3 = dt_10;

#line 1097
                affineLength_9 = affineLength_10;

#line 1097
                accepted_4 = accepted_5;

#line 719
                int _S968 = step_8 + int(1);

#line 719
                pupilRight_1 = transmittance_1;

#line 719
                screenDirection_1 = color_4;

#line 719
                step_8 = _S968;

#line 719
                attempts_1 = attempts_3;

#line 719
            }

#line 719
            pos0_0 = compX_2;

#line 719
            compU_2 = compU_1;

#line 698
        }

#line 1102
        if(reason_1 == 0U)
        {

#line 1102
            escaped_0 = attempts_1 >= uint(maxSteps_0);

#line 1102
        }
        else
        {

#line 1102
            escaped_0 = false;

#line 1102
        }

#line 1102
        if(escaped_0)
        {

#line 1102
            reason_1 = 4U;

#line 1102
        }

#line 1102
        thread Vec4Cart_0 _S1040 = (&state_15)->x_5;

#line 1102
        float4 _S1041 = TracePack_0(&_S1040);

        (&record_1)->x_1 = _S1041;

#line 1104
        thread Vec4Cart_0 _S1042 = (&state_15)->u_1;

#line 1104
        float4 _S1043 = TracePack_0(&_S1042);

#line 1104
        (&record_1)->u_0 = _S1043;

#line 1104
        column_0 = int(0);
        for(;;)
        {

#line 1105
            if(column_0 < int(4))
            {
            }
            else
            {

#line 1105
                break;
            }

#line 1106
            (&record_1)->X_0[column_0] = X_1[column_0];

#line 1106
            (&record_1)->V_0[column_0] = V_2[column_0];

#line 1105
            column_0 = column_0 + int(1);

#line 1105
        }

#line 1105
        thread Vec4Cart_0 _S1044 = pos0_0;

#line 1105
        float4 _S1045 = TracePack_0(&_S1044);


        (&record_1)->compX_0 = _S1045;

#line 1108
        thread Vec4Cart_0 _S1046 = compU_2;

#line 1108
        float4 _S1047 = TracePack_0(&_S1046);

#line 1108
        (&record_1)->compU_0 = _S1047;
        (&record_1)->integration_0 = float4(pupilUp_1, affineLength_9, lastStep_3, 0.00100000004749745);
        (&record_1)->source_0 = float4(color_4, transmittance_1);
        (&record_1)->volume_0 = float4(volumeEmission_0, volumeOpticalDepth_0);

#line 1111
        uint _S1048;
        if(reason_1 == 0U)
        {

#line 1112
            _S1048 = 1U;

#line 1112
        }
        else
        {

#line 1112
            _S1048 = 2U;

#line 1112
        }

#line 1112
        (&record_1)->control_0 = uint4(_S1048, reason_1, attempts_1, accepted_4);

        (&record_1)->identity_0 = identity_1;

#line 1114
        thread TraceContinuation_0 _S1049 = record_1;

#line 1114
        bool _S1050 = TraceStateFinite_0(&_S1049);
        if(!_S1050)
        {

#line 1116
            (&record_1)->control_0.x = 2U;
            (&record_1)->control_0.y = 9U;

#line 1115
        }

#line 1115
        TraceContinuation_natural_0 device* _S1051 = (&kernelContext_6)->continuation_0+slot_0;

#line 1115
        _S1051->x_1 = packed_float4(record_1.x_1) ;

#line 1115
        _S1051->u_0 = packed_float4(record_1.u_0) ;

#line 1115
        (&_S1051->X_0)->data_0[int(0)] = packed_float4(record_1.X_0[int(0)]) ;

#line 1115
        (&_S1051->X_0)->data_0[int(1)] = packed_float4(record_1.X_0[int(1)]) ;

#line 1115
        (&_S1051->X_0)->data_0[int(2)] = packed_float4(record_1.X_0[int(2)]) ;

#line 1115
        (&_S1051->X_0)->data_0[int(3)] = packed_float4(record_1.X_0[int(3)]) ;

#line 1115
        (&_S1051->V_0)->data_0[int(0)] = packed_float4(record_1.V_0[int(0)]) ;

#line 1115
        (&_S1051->V_0)->data_0[int(1)] = packed_float4(record_1.V_0[int(1)]) ;

#line 1115
        (&_S1051->V_0)->data_0[int(2)] = packed_float4(record_1.V_0[int(2)]) ;

#line 1115
        (&_S1051->V_0)->data_0[int(3)] = packed_float4(record_1.V_0[int(3)]) ;

#line 1115
        _S1051->compX_0 = packed_float4(record_1.compX_0) ;

#line 1115
        _S1051->compU_0 = packed_float4(record_1.compU_0) ;

#line 1115
        _S1051->integration_0 = packed_float4(record_1.integration_0) ;

#line 1115
        _S1051->source_0 = packed_float4(record_1.source_0) ;

#line 1115
        _S1051->volume_0 = packed_float4(record_1.volume_0) ;

#line 1115
        _S1051->control_0 = packed_uint4(record_1.control_0) ;

#line 1115
        _S1051->identity_0 = packed_uint4(record_1.identity_0) ;

#line 1120
        return;
    }
    if(reason_2 == 3U)
    {

#line 1122
        escaped_0 = true;

#line 1122
    }
    else
    {

#line 1122
        escaped_0 = reason_2 == 10U;

#line 1122
    }
    bool oppositeSheet_0 = reason_2 == 10U;
    bool captured_0 = reason_2 == 1U;

    float3 _S1052 = float3(0.0, 0.0, 1.0);
    thread float beamOrientation_0 = 0.0;

#line 1134
    if(escaped_0)
    {

#line 1134
        terminalEscapeSurface_0 = true;

#line 1134
    }
    else
    {

#line 1134
        if(beamsEnabled_0)
        {

#line 1134
            terminalEscapeSurface_0 = !captured_0;

#line 1134
        }
        else
        {

#line 1134
            terminalEscapeSurface_0 = false;

#line 1134
        }

#line 1134
    }

#line 1134
    if(terminalEscapeSurface_0)
    {

#line 1135
        thread array<array<float, int(4)>, int(4)> terminalMetric_0;
        thread array<array<float, int(4)>, int(4)> terminalInverse_0;
        if(isWarp_5)
        {

#line 1137
            thread Vec4Cart_0 _S1053 = (&state_15)->x_5;

#line 1137
            GetWarpDriveMetric_0(&_S1053, warpVs_5, warpSigma_5, warpR_5, 0.0, 0.0, 0.0, &terminalMetric_0, &terminalInverse_0);

#line 1137
        }
        else
        {
            if(isWormhole_5)
            {

#line 1140
                thread Vec4Cart_0 _S1054 = (&state_15)->x_5;

#line 1140
                GetMorrisThorneCartesianMetric_0(&_S1054, throatB0_5, 0.0, wormShape_8, &terminalMetric_0, &terminalInverse_0);

#line 1140
            }
            else
            {

#line 1140
                thread Vec4Cart_0 _S1055 = (&state_15)->x_5;

#line 1140
                GetKerrSchildFamilyMetric_0(&_S1055, M_20, a_29, Q_13, Lambda_13, &terminalMetric_0, &terminalInverse_0);

#line 1140
            }

#line 1137
        }

#line 1147
        float3 _S1056 = float3(1.0, 0.0, 0.0);
        float3 _S1057 = float3(0.0, -1.0, 0.0);

#line 1148
        float3 _S1058 = float3(0.0) ;

#line 1148
        thread array<array<float, int(4)>, int(4)> _S1059 = terminalMetric_0;

#line 1148
        thread array<array<float, int(4)>, int(4)> _S1060 = terminalInverse_0;

#line 1148
        CameraObserverFrame_0 _S1061 = BuildCameraObserverFrame_0(&_S1059, &_S1060, _S1056, _S1057, _S1052, _S1058);

#line 1148
        thread array<array<float, int(4)>, int(4)> _S1062 = terminalMetric_0;

#line 1148
        thread Vec4Cart_0 _S1063 = (&state_15)->u_1;

#line 1148
        thread Vec4Cart_0 _S1064 = _S1061.time_0;

#line 1148
        float _S1065 = MetricDotCart_0(&_S1062, &_S1063, &_S1064);

        if(isfinite(_S1065))
        {

#line 1150
            terminalEscapeSurface_0 = _S1065 > 0.0;

#line 1150
        }
        else
        {

#line 1150
            terminalEscapeSurface_0 = false;

#line 1150
        }

#line 1150
        if(!terminalEscapeSurface_0)
        {

#line 1151
            (&record_1)->control_0.x = 2U;

#line 1151
            (&record_1)->control_0.y = 9U;

#line 1151
            TraceContinuation_natural_0 device* _S1066 = (&kernelContext_6)->continuation_0+slot_0;

#line 1151
            _S1066->x_1 = packed_float4(record_1.x_1) ;

#line 1151
            _S1066->u_0 = packed_float4(record_1.u_0) ;

#line 1151
            (&_S1066->X_0)->data_0[int(0)] = packed_float4(record_1.X_0[int(0)]) ;

#line 1151
            (&_S1066->X_0)->data_0[int(1)] = packed_float4(record_1.X_0[int(1)]) ;

#line 1151
            (&_S1066->X_0)->data_0[int(2)] = packed_float4(record_1.X_0[int(2)]) ;

#line 1151
            (&_S1066->X_0)->data_0[int(3)] = packed_float4(record_1.X_0[int(3)]) ;

#line 1151
            (&_S1066->V_0)->data_0[int(0)] = packed_float4(record_1.V_0[int(0)]) ;

#line 1151
            (&_S1066->V_0)->data_0[int(1)] = packed_float4(record_1.V_0[int(1)]) ;

#line 1151
            (&_S1066->V_0)->data_0[int(2)] = packed_float4(record_1.V_0[int(2)]) ;

#line 1151
            (&_S1066->V_0)->data_0[int(3)] = packed_float4(record_1.V_0[int(3)]) ;

#line 1151
            _S1066->compX_0 = packed_float4(record_1.compX_0) ;

#line 1151
            _S1066->compU_0 = packed_float4(record_1.compU_0) ;

#line 1151
            _S1066->integration_0 = packed_float4(record_1.integration_0) ;

#line 1151
            _S1066->source_0 = packed_float4(record_1.source_0) ;

#line 1151
            _S1066->volume_0 = packed_float4(record_1.volume_0) ;

#line 1151
            _S1066->control_0 = packed_uint4(record_1.control_0) ;

#line 1151
            _S1066->identity_0 = packed_uint4(record_1.identity_0) ;
            return;
        }

#line 1152
        thread array<array<float, int(4)>, int(4)> _S1067 = terminalMetric_0;

#line 1152
        thread Vec4Cart_0 _S1068 = (&state_15)->u_1;

#line 1152
        thread Vec4Cart_0 _S1069 = _S1061.forward_0;

#line 1152
        float _S1070 = MetricDotCart_0(&_S1067, &_S1068, &_S1069);


        float _S1071 = _S1070 / _S1065;

#line 1155
        thread array<array<float, int(4)>, int(4)> _S1072 = terminalMetric_0;

#line 1155
        thread Vec4Cart_0 _S1073 = (&state_15)->u_1;

#line 1155
        thread Vec4Cart_0 _S1074 = _S1061.up_0;

#line 1155
        float _S1075 = MetricDotCart_0(&_S1072, &_S1073, &_S1074);

        float _S1076 = _S1075 / _S1065;

#line 1157
        thread array<array<float, int(4)>, int(4)> _S1077 = terminalMetric_0;

#line 1157
        thread Vec4Cart_0 _S1078 = (&state_15)->u_1;

#line 1157
        thread Vec4Cart_0 _S1079 = _S1061.right_0;

#line 1157
        float _S1080 = MetricDotCart_0(&_S1077, &_S1078, &_S1079);

#line 1154
        float3 terminalDirection_0 = normalize(float3(_S1071, _S1076, _S1080 / _S1065));

#line 1160
        if(!all(isfinite(terminalDirection_0)))
        {

#line 1161
            (&record_1)->control_0.x = 2U;

#line 1161
            (&record_1)->control_0.y = 9U;

#line 1161
            TraceContinuation_natural_0 device* _S1081 = (&kernelContext_6)->continuation_0+slot_0;

#line 1161
            _S1081->x_1 = packed_float4(record_1.x_1) ;

#line 1161
            _S1081->u_0 = packed_float4(record_1.u_0) ;

#line 1161
            (&_S1081->X_0)->data_0[int(0)] = packed_float4(record_1.X_0[int(0)]) ;

#line 1161
            (&_S1081->X_0)->data_0[int(1)] = packed_float4(record_1.X_0[int(1)]) ;

#line 1161
            (&_S1081->X_0)->data_0[int(2)] = packed_float4(record_1.X_0[int(2)]) ;

#line 1161
            (&_S1081->X_0)->data_0[int(3)] = packed_float4(record_1.X_0[int(3)]) ;

#line 1161
            (&_S1081->V_0)->data_0[int(0)] = packed_float4(record_1.V_0[int(0)]) ;

#line 1161
            (&_S1081->V_0)->data_0[int(1)] = packed_float4(record_1.V_0[int(1)]) ;

#line 1161
            (&_S1081->V_0)->data_0[int(2)] = packed_float4(record_1.V_0[int(2)]) ;

#line 1161
            (&_S1081->V_0)->data_0[int(3)] = packed_float4(record_1.V_0[int(3)]) ;

#line 1161
            _S1081->compX_0 = packed_float4(record_1.compX_0) ;

#line 1161
            _S1081->compU_0 = packed_float4(record_1.compU_0) ;

#line 1161
            _S1081->integration_0 = packed_float4(record_1.integration_0) ;

#line 1161
            _S1081->source_0 = packed_float4(record_1.source_0) ;

#line 1161
            _S1081->volume_0 = packed_float4(record_1.volume_0) ;

#line 1161
            _S1081->control_0 = packed_uint4(record_1.control_0) ;

#line 1161
            _S1081->identity_0 = packed_uint4(record_1.identity_0) ;
            return;
        }
        if(oppositeSheet_0)
        {

#line 1164
            float3 _S1082;

#line 1164
            thread Vec4Cart_0 _S1083 = (&state_15)->x_5;
            thread float3 mappedDirection_1;

#line 1165
            bool _S1084 = MapEllisSecondSheetSkyDirection_0(&_S1083, terminalDirection_0, &mappedDirection_1);
            if(_S1084)
            {

#line 1166
                _S1082 = normalize(mappedDirection_1);

#line 1166
            }
            else
            {

                (&record_1)->control_0.x = 2U;
                (&record_1)->control_0.y = 6U;

#line 1171
                TraceContinuation_natural_0 device* _S1085 = (&kernelContext_6)->continuation_0+slot_0;

#line 1171
                _S1085->x_1 = packed_float4(record_1.x_1) ;

#line 1171
                _S1085->u_0 = packed_float4(record_1.u_0) ;

#line 1171
                (&_S1085->X_0)->data_0[int(0)] = packed_float4(record_1.X_0[int(0)]) ;

#line 1171
                (&_S1085->X_0)->data_0[int(1)] = packed_float4(record_1.X_0[int(1)]) ;

#line 1171
                (&_S1085->X_0)->data_0[int(2)] = packed_float4(record_1.X_0[int(2)]) ;

#line 1171
                (&_S1085->X_0)->data_0[int(3)] = packed_float4(record_1.X_0[int(3)]) ;

#line 1171
                (&_S1085->V_0)->data_0[int(0)] = packed_float4(record_1.V_0[int(0)]) ;

#line 1171
                (&_S1085->V_0)->data_0[int(1)] = packed_float4(record_1.V_0[int(1)]) ;

#line 1171
                (&_S1085->V_0)->data_0[int(2)] = packed_float4(record_1.V_0[int(2)]) ;

#line 1171
                (&_S1085->V_0)->data_0[int(3)] = packed_float4(record_1.V_0[int(3)]) ;

#line 1171
                _S1085->compX_0 = packed_float4(record_1.compX_0) ;

#line 1171
                _S1085->compU_0 = packed_float4(record_1.compU_0) ;

#line 1171
                _S1085->integration_0 = packed_float4(record_1.integration_0) ;

#line 1171
                _S1085->source_0 = packed_float4(record_1.source_0) ;

#line 1171
                _S1085->volume_0 = packed_float4(record_1.volume_0) ;

#line 1171
                _S1085->control_0 = packed_uint4(record_1.control_0) ;

#line 1171
                _S1085->identity_0 = packed_uint4(record_1.identity_0) ;

                return;
            }

#line 1173
            screenDirection_1 = _S1082;

#line 1164
        }
        else
        {

#line 1164
            screenDirection_1 = terminalDirection_0;

#line 1164
        }

#line 1164
        thread CameraObserverFrame_0 _S1086 = _S1061;

#line 1176
        thread Vec4Cart_0 terminalScreenA_0;
        thread Vec4Cart_0 terminalScreenB_0;

#line 1177
        SachsScreenBasis_0(&_S1086, screenDirection_1, &terminalScreenA_0, &terminalScreenB_0);

        if(pointStarfieldEnabled_0)
        {

#line 1179
            step_8 = int(0);

#line 1179
        }
        else
        {

#line 1179
            step_8 = int(2);

#line 1179
        }
        Vec4Cart_0 xiACart_0 = TraceUnpack_0(X_1[step_8]);
        Vec4Cart_0 xiBCart_0 = TraceUnpack_0(X_1[step_8 + int(1)]);

#line 1181
        thread array<array<float, int(4)>, int(4)> _S1087 = terminalMetric_0;

#line 1181
        thread Vec4Cart_0 _S1088 = xiACart_0;

#line 1181
        thread Vec4Cart_0 _S1089 = terminalScreenA_0;

#line 1181
        float _S1090 = MetricDotCart_0(&_S1087, &_S1088, &_S1089);

#line 1181
        thread array<array<float, int(4)>, int(4)> _S1091 = terminalMetric_0;

#line 1181
        thread Vec4Cart_0 _S1092 = xiBCart_0;

#line 1181
        thread Vec4Cart_0 _S1093 = terminalScreenA_0;

#line 1181
        float _S1094 = MetricDotCart_0(&_S1091, &_S1092, &_S1093);

#line 1181
        thread array<array<float, int(4)>, int(4)> _S1095 = terminalMetric_0;

#line 1181
        thread Vec4Cart_0 _S1096 = xiACart_0;

#line 1181
        thread Vec4Cart_0 _S1097 = terminalScreenB_0;

#line 1181
        float _S1098 = MetricDotCart_0(&_S1095, &_S1096, &_S1097);

#line 1181
        thread array<array<float, int(4)>, int(4)> _S1099 = terminalMetric_0;

#line 1181
        thread Vec4Cart_0 _S1100 = xiBCart_0;

#line 1181
        thread Vec4Cart_0 _S1101 = terminalScreenB_0;

#line 1181
        float _S1102 = MetricDotCart_0(&_S1099, &_S1100, &_S1101);

#line 1186
        thread float beamSemiMajor_0;
        thread float beamSemiMinor_0;
        thread float beamDeterminant_0;
        BeamEllipseFromTransverseMap_0(_S1090, _S1094, _S1098, _S1102, &beamSemiMajor_0, &beamSemiMinor_0, &beamOrientation_0, &beamDeterminant_0);



        float terminalRadius_0 = length(float3((&(&state_15)->x_5)->x_2, (&(&state_15)->x_5)->y_0, (&(&state_15)->x_5)->z_0));

        if(pointStarfieldEnabled_0)
        {

#line 1195
            pupilUp_1 = max(terminalRadius_0 * 0.00100000004749745, 9.999999960041972e-13);

#line 1195
        }
        else
        {

#line 1195
            pupilUp_1 = 0.00100000004749745;

#line 1195
        }

        float beamMajorExpansion_0 = beamSemiMajor_0 / pupilUp_1;
        float beamMinorExpansion_0 = beamSemiMinor_0 / pupilUp_1;
        float _S1103 = sqrt(max(beamMajorExpansion_0 * beamMinorExpansion_0, 0.0));

#line 1199
        pupilUp_1 = beamMajorExpansion_0;

#line 1199
        pupilRight_1 = beamMinorExpansion_0;

#line 1199
        lastStep_3 = _S1103;

#line 1134
    }
    else
    {

#line 1134
        screenDirection_1 = _S1052;

#line 1134
        pupilUp_1 = 1.0;

#line 1134
        pupilRight_1 = 1.0;

#line 1134
        lastStep_3 = 1.0;

#line 1134
    }

#line 1205
    if(escaped_0)
    {

        if(beamsEnabled_0)
        {

#line 1208
            affineLength_9 = max(pixelAngularSize_0 * pupilUp_1, pixelAngularSize_0);

#line 1208
        }
        else
        {

#line 1208
            affineLength_9 = 0.30000001192092896 * pixelAngularSize_0;

#line 1208
        }


        if(beamsEnabled_0)
        {

#line 1211
            exoticScale_0 = max(pixelAngularSize_0 * pupilRight_1, pixelAngularSize_0);

#line 1211
        }
        else
        {

#line 1211
            exoticScale_0 = 0.30000001192092896 * pixelAngularSize_0;

#line 1211
        }


        if(pointStarfieldEnabled_0)
        {

#line 1214
            float3 _S1104 = SamplePointStarfield_0(screenDirection_1, affineLength_9, exoticScale_0, beamOrientation_0, pointBrightnessScale_0, &kernelContext_6);

#line 1214
            color_4 = _S1104;

#line 1214
        }
        else
        {
            if(starfieldEnabled_0)
            {

#line 1217
                float3 _S1105 = SampleStarfield_0(screenDirection_1, starfieldWidth_0, starfieldHeight_0, &kernelContext_6);

#line 1217
                screenDirection_1 = _S1105;

#line 1217
            }
            else
            {

#line 1217
                screenDirection_1 = SampleBackgroundGradient_0(screenDirection_1);

#line 1217
            }

#line 1217
            color_4 = screenDirection_1;

#line 1214
        }

#line 1214
        color_4 = color_3 + float3(transmittance_0)  * color_4;

#line 1205
    }
    else
    {

#line 1205
        color_4 = color_3;

#line 1205
    }

#line 1221
    if(volumeOpticalDepth_0 > 0.0)
    {

#line 1221
        color_4 = color_4 * float3(exp(- volumeOpticalDepth_0))  + volumeEmission_0;

#line 1221
    }

#line 1229
    if(beamsEnabled_0)
    {

#line 1229
        escaped_0 = !captured_0;

#line 1229
    }
    else
    {

#line 1229
        escaped_0 = false;

#line 1229
    }

#line 1229
    if(escaped_0)
    {

#line 1229
        affineLength_9 = lastStep_3;

#line 1229
    }
    else
    {

#line 1229
        affineLength_9 = 1.0;

#line 1229
    }

#line 1234
    if(sampleIndex_0 > 0.5)
    {
        float inverseCount_0 = 1.0 / (sampleIndex_0 + 1.0);
        int _S1106 = _S881 * int(4);



        float3 _S1107 = (float3(*((&kernelContext_6)->radiance_0+_S1106), *((&kernelContext_6)->radiance_0+(_S1106 + int(1))), *((&kernelContext_6)->radiance_0+(_S1106 + int(2)))) * float3(sampleIndex_0)  + color_4) * float3(inverseCount_0) ;

#line 1241
        affineLength_9 = (*((&kernelContext_6)->radiance_0+(_S1106 + int(3))) * sampleIndex_0 + affineLength_9) * inverseCount_0;

#line 1241
        color_4 = _S1107;

#line 1234
    }

#line 1244
    if(!all(isfinite(color_4)))
    {

#line 1244
        escaped_0 = true;

#line 1244
    }
    else
    {

#line 1244
        escaped_0 = !isfinite(affineLength_9);

#line 1244
    }

#line 1244
    if(escaped_0)
    {

#line 1245
        (&record_1)->control_0.x = 2U;

#line 1245
        (&record_1)->control_0.y = 9U;

#line 1245
        TraceContinuation_natural_0 device* _S1108 = (&kernelContext_6)->continuation_0+slot_0;

#line 1245
        _S1108->x_1 = packed_float4(record_1.x_1) ;

#line 1245
        _S1108->u_0 = packed_float4(record_1.u_0) ;

#line 1245
        (&_S1108->X_0)->data_0[int(0)] = packed_float4(record_1.X_0[int(0)]) ;

#line 1245
        (&_S1108->X_0)->data_0[int(1)] = packed_float4(record_1.X_0[int(1)]) ;

#line 1245
        (&_S1108->X_0)->data_0[int(2)] = packed_float4(record_1.X_0[int(2)]) ;

#line 1245
        (&_S1108->X_0)->data_0[int(3)] = packed_float4(record_1.X_0[int(3)]) ;

#line 1245
        (&_S1108->V_0)->data_0[int(0)] = packed_float4(record_1.V_0[int(0)]) ;

#line 1245
        (&_S1108->V_0)->data_0[int(1)] = packed_float4(record_1.V_0[int(1)]) ;

#line 1245
        (&_S1108->V_0)->data_0[int(2)] = packed_float4(record_1.V_0[int(2)]) ;

#line 1245
        (&_S1108->V_0)->data_0[int(3)] = packed_float4(record_1.V_0[int(3)]) ;

#line 1245
        _S1108->compX_0 = packed_float4(record_1.compX_0) ;

#line 1245
        _S1108->compU_0 = packed_float4(record_1.compU_0) ;

#line 1245
        _S1108->integration_0 = packed_float4(record_1.integration_0) ;

#line 1245
        _S1108->source_0 = packed_float4(record_1.source_0) ;

#line 1245
        _S1108->volume_0 = packed_float4(record_1.volume_0) ;

#line 1245
        _S1108->control_0 = packed_uint4(record_1.control_0) ;

#line 1245
        _S1108->identity_0 = packed_uint4(record_1.identity_0) ;
        return;
    }
    int _S1109 = _S881 * int(4);

#line 1248
    *((&kernelContext_6)->radiance_0+_S1109) = color_4.x;
    *((&kernelContext_6)->radiance_0+(_S1109 + int(1))) = color_4.y;
    *((&kernelContext_6)->radiance_0+(_S1109 + int(2))) = color_4.z;
    *((&kernelContext_6)->radiance_0+(_S1109 + int(3))) = affineLength_9;
    (&record_1)->control_0.x = 3U;

#line 1252
    TraceContinuation_natural_0 device* _S1110 = (&kernelContext_6)->continuation_0+slot_0;

#line 1252
    _S1110->x_1 = packed_float4(record_1.x_1) ;

#line 1252
    _S1110->u_0 = packed_float4(record_1.u_0) ;

#line 1252
    (&_S1110->X_0)->data_0[int(0)] = packed_float4(record_1.X_0[int(0)]) ;

#line 1252
    (&_S1110->X_0)->data_0[int(1)] = packed_float4(record_1.X_0[int(1)]) ;

#line 1252
    (&_S1110->X_0)->data_0[int(2)] = packed_float4(record_1.X_0[int(2)]) ;

#line 1252
    (&_S1110->X_0)->data_0[int(3)] = packed_float4(record_1.X_0[int(3)]) ;

#line 1252
    (&_S1110->V_0)->data_0[int(0)] = packed_float4(record_1.V_0[int(0)]) ;

#line 1252
    (&_S1110->V_0)->data_0[int(1)] = packed_float4(record_1.V_0[int(1)]) ;

#line 1252
    (&_S1110->V_0)->data_0[int(2)] = packed_float4(record_1.V_0[int(2)]) ;

#line 1252
    (&_S1110->V_0)->data_0[int(3)] = packed_float4(record_1.V_0[int(3)]) ;

#line 1252
    _S1110->compX_0 = packed_float4(record_1.compX_0) ;

#line 1252
    _S1110->compU_0 = packed_float4(record_1.compU_0) ;

#line 1252
    _S1110->integration_0 = packed_float4(record_1.integration_0) ;

#line 1252
    _S1110->source_0 = packed_float4(record_1.source_0) ;

#line 1252
    _S1110->volume_0 = packed_float4(record_1.volume_0) ;

#line 1252
    _S1110->control_0 = packed_uint4(record_1.control_0) ;

#line 1252
    _S1110->identity_0 = packed_uint4(record_1.identity_0) ;

    return;
}

