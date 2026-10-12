// Dedicated fp64 diagnostic for the existing nine-packet coupled regression.
// The original float buffer ABI and opcode 25 contract are retained.
import gr_types;
import gr_metrics;
import gr_integrator;
import gr_coupled;

StructuredBuffer<float> samples;     // binding 0
RWStructuredBuffer<float> results;   // binding 1
StructuredBuffer<float> params;      // binding 2

static const uint OP_COUPLED_STEP = 25;

[shader("compute")]
[numthreads(64, 1, 1)]
void ComputeMain(uint3 tid: SV_DispatchThreadID) {
    uint opcode = uint(params[0]);
    uint count = uint(params[1]);
    uint sStride = uint(params[2]);
    uint rStride = uint(params[3]);
    uint i = tid.x;
    if (i >= count) return;

    uint sb = i * sStride;
    uint rb = i * rStride;

    uint metricId = uint(samples[sb + 0]);
    float p1 = samples[sb + 1];
    float p2 = samples[sb + 2];
    float p3 = samples[sb + 3];
    float p4 = samples[sb + 4];
    float c0 = samples[sb + 5];
    float c1 = samples[sb + 6];
    float c2 = samples[sb + 7];

    if (opcode == OP_COUPLED_STEP) {
        // Actual same-stage four-column production candidate versus independently
        // retraced scalar neighbours. q is held fixed through each whole trial;
        // no finite difference of the candidate's own variation is used.
        CoupledState raw=(CoupledState)0;
        raw.x=Real4(Real(0),Real(c0),Real(c1),Real(c2));
        raw.k=Real4(Real(samples[sb+8]),Real(samples[sb+9]),
                    Real(samples[sb+10]),Real(samples[sb+11]));
        raw.K[0].z=Real(1);raw.K[1].w=Real(1);
        raw.X[2].z=Real(1);raw.X[3].w=Real(1);
        Real h=Real(samples[sb+12]);Real epsilon=Real(samples[sb+13]);
        for(uint j=0;j<64;++j) results[rb+j]=0.0f;
        CoupledState initial;Real4 initialV[4];int selected;
        bool represented=TryProjectCoupled(raw,p1,p2,p3,p4,initial,initialV,selected);
        results[rb+7]=represented?1.0f:0.0f;
        results[rb+54]=float(selected);
        if(!represented) { results[rb]=float(kCoupledDerivativeDomain); return; }
        CoupledControl control;
        control.lengthScale=p1>0.0f?Real(p1):length(raw.x.yzw);
        control.frequencyScale=Real(1);
        control.variationTolerance=Real(1.0e-4)/(Real(4)*Real(100000));
        control.centralTolerance=AdaptiveIntegrationTolerance();
        control.nullTolerance=AdaptiveNullTolerance();
        CoupledTrial trial=AttemptCoupledGeodesicStep(initial,h,p1,p2,p3,p4,control);
        results[rb]=float(trial.status);results[rb+1]=float(trial.stages);
        results[rb+2]=float(trial.errorRatio);results[rb+3]=float(trial.unprojectedError);
        results[rb+4]=float(trial.unprojectedNull);results[rb+5]=float(trial.projectedNull);
        if(trial.status!=kCoupledNone&&trial.status!=kCoupledAboveBudget) return;
        AdaptiveGeodesicStep scalar=AttemptAdaptiveGeodesicStep(CoupledCentral(initial),h,
            false,false,p1,p2,p3,p4,0,1,1,1,0,control.centralTolerance,control.nullTolerance);
        Real4 scalarX=CoupledPack(scalar.candidate.x),scalarK=CoupledPack(scalar.candidate.u);
        Real scalarDifference=Real(0);
        for(int mu=0;mu<4;++mu) {
            scalarDifference=max(scalarDifference,abs(scalarX[mu]-trial.fine.x[mu]));
            scalarDifference=max(scalarDifference,abs(scalarK[mu]-trial.fine.k[mu]));
        }
        results[rb+6]=float(scalarDifference);
        if(!(epsilon>Real(0))||!isfinite(epsilon)) return;
        for(int column=0;column<4;++column) {
            GeodesicStateCart neighbours[2];
            bool good=true;
            for(int side=0;side<2;++side) {
                Real sign=side==0?Real(-1):Real(1);
                GeodesicStateCart launch;
                launch.x=CoupledVector(raw.x+sign*epsilon*raw.X[column]);
                launch.u=CoupledVector(raw.k+sign*epsilon*raw.K[column]);
                Real g[4][4],ginv[4][4];
                GetKerrSchildFamilyMetric(launch.x,p1,p2,p3,p4,g,ginv);
                Vec4Cart projected;
                if(!ProjectNullTangentPreservingBranch(g,launch.u,projected)) { good=false;break; }
                launch.u=projected;
                AdaptiveGeodesicStep neighbour=AttemptAdaptiveGeodesicStep(launch,h,
                    false,false,p1,p2,p3,p4,0,1,1,1,0,control.centralTolerance,control.nullTolerance);
                // Inspect the actual fixed-h projected map, independently of
                // acceptance; above-budget refusal is reported separately.
                neighbours[side]=neighbour.candidate;
            }
            if(!good) return;
            Real4 dX=(CoupledPack(neighbours[1].x)-CoupledPack(neighbours[0].x))/(Real(2)*epsilon);
            Real4 dK=(CoupledPack(neighbours[1].u)-CoupledPack(neighbours[0].u))/(Real(2)*epsilon);
            for(int mu=0;mu<4;++mu) {
                results[rb+8+column*8+mu]=float(dX[mu]-trial.fine.X[column][mu]);
                results[rb+12+column*8+mu]=float(dK[mu]-trial.fine.K[column][mu]);
            }
            results[rb+40+column*2]=float(length(trial.fine.X[column]));
            results[rb+41+column*2]=float(length(trial.fine.K[column]));
        }
        results[rb+48]=1.0f; // All four independent neighbours represented.
        // Contract the actual differentiated final null constraint before any
        // float output narrowing. This is diagnostic, not the only oracle.
        CoupledJet4 metric[4][4],connection[4][4][4];
        if(TryKerrSchildConnectionJet(CoupledVector(trial.fine.x),trial.fine.X,
                                     p1,p2,p3,p4,metric,connection)) {
            for(int column=0;column<4;++column) {
                Real residual=Real(0),scale=Real(0);
                for(int mu=0;mu<4;++mu) for(int nu=0;nu<4;++nu) {
                    Real aTerm=metric[mu][nu].derivative[column]*trial.fine.k[mu]*trial.fine.k[nu];
                    Real bTerm=metric[mu][nu].value*(trial.fine.K[column][mu]*trial.fine.k[nu]+
                        trial.fine.k[mu]*trial.fine.K[column][nu]);
                    residual+=aTerm+bTerm;scale+=abs(aTerm)+abs(bTerm);
                }
                results[rb+49+column]=float(abs(residual)/max(Real(1),scale));
            }
            results[rb+53]=1.0f;
        }
    }
}
