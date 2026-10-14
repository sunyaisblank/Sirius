// Joint RK4: same-stage coordinate variations of the actual central
// acceleration. Dense/event acceptance and production publication are separate
// required integration work; a successful local trial is not a source map.
module gr_coupled;
import gr_types;
import gr_metrics;
import gr_integrator;

public static const uint kCoupledNone=0;
public static const uint kCoupledInvalid=1;
public static const uint kCoupledDerivativeDomain=2;
public static const uint kCoupledProjection=3;
public static const uint kCoupledAboveBudget=4;

public struct CoupledState {
    public Real4 x;
    public Real4 k;
    public Real4 X[4];
    public Real4 K[4]; // Coordinate derivative, NOT covariant D X/dlambda.
}
public struct CoupledIncrement {
    public Real4 dx;
    public Real4 dX[4];
}
public struct CoupledControl {
    public Real lengthScale;
    public Real frequencyScale;
    public Real variationTolerance;
    public Real centralTolerance;
    public Real nullTolerance;
}
public struct CoupledTrial {
    public CoupledState full;
    public CoupledState midpoint;
    public CoupledState fine;
    public CoupledIncrement fullIncrement;
    public CoupledIncrement midpointIncrement;
    public CoupledIncrement fineIncrement;
    public uint status;
    public uint stages;
    public Real errorRatio;
    public Real unprojectedError;
    public Real unprojectedNull;
    public Real projectedNull;
    public Real stepFactor;
}
public Vec4Cart CoupledVector(Real4 value) {
    return MakeVec4Cart(value.x,value.y,value.z,value.w);
}
public Real4 CoupledPack(Vec4Cart value) {
    return Real4(value.t,value.x,value.y,value.z);
}
public GeodesicStateCart CoupledCentral(CoupledState state) {
    GeodesicStateCart result;result.x=CoupledVector(state.x);result.u=CoupledVector(state.k);
    return result;
}
public bool CoupledFinite(CoupledState state) {
    if(!all(isfinite(state.x))||!all(isfinite(state.k))) return false;
    for(int i=0;i<4;++i) if(!all(isfinite(state.X[i]))||!all(isfinite(state.K[i]))) return false;
    return true;
}
public bool CoupledControlRepresented(CoupledControl control) {
    return isfinite(control.lengthScale)&&control.lengthScale>Real(0)&&
        isfinite(control.frequencyScale)&&control.frequencyScale>Real(0)&&
        isfinite(control.variationTolerance)&&control.variationTolerance>Real(0)&&
        isfinite(control.centralTolerance)&&control.centralTolerance>Real(0)&&
        isfinite(control.nullTolerance)&&control.nullTolerance>Real(0);
}

public bool TryCoupledDerivative(CoupledState state,float M,float a,float Q,float Lambda,
                                  out CoupledState slope) {
    slope=(CoupledState)0;
    if(!CoupledFinite(state)) return false;
    CoupledJet4 metric[4][4];
    CoupledJet4 connection[4][4][4];
    if(!TryKerrSchildConnectionJet(CoupledVector(state.x),state.X,M,a,Q,Lambda,
                                  metric,connection)) return false;
    slope.x=state.k;
    for(int i=0;i<4;++i) slope.X[i]=state.K[i];
    for(int mu=0;mu<4;++mu) for(int nu=0;nu<4;++nu) for(int rho=0;rho<4;++rho) {
        CoupledJet4 gamma=connection[mu][nu][rho];
        slope.k[mu]-=gamma.value*state.k[nu]*state.k[rho];
        for(int i=0;i<4;++i) {
            slope.K[i][mu]-=gamma.derivative[i]*state.k[nu]*state.k[rho]+
                gamma.value*(state.K[i][nu]*state.k[rho]+state.k[nu]*state.K[i][rho]);
        }
    }
    return CoupledFinite(slope);
}
public CoupledState CoupledAdd(CoupledState state,CoupledState slope,Real step,Real weight) {
    CoupledState result;
    result.x=state.x+(slope.x*step)*weight;result.k=state.k+(slope.k*step)*weight;
    for(int i=0;i<4;++i) {
        result.X[i]=state.X[i]+(slope.X[i]*step)*weight;
        result.K[i]=state.K[i]+(slope.K[i]*step)*weight;
    }
    return result;
}

// Protect both the retained final product and its endpoint addition.
// The captured arithmetic control requires both decorations on this route.
// This does not declare a global rounding mode or change stage arithmetic.
public Real4 CoupledAddRetainedIncrement(Real4 start,Real4 increment) {
    return spirv_asm {
        OpFAdd $$Real4 %sum $start $increment;
        OpDecorate %sum NoContraction;
        OpCopyObject $$Real4 result %sum;
    };
}

public Real4 CoupledRetainedProduct(Real4 weighted,Real sixth) {
    Real4 scale=Real4(sixth);
    return spirv_asm {
        OpFMul $$Real4 %product $weighted $scale;
        OpDecorate %product NoContraction;
        OpCopyObject $$Real4 result %product;
    };
}

// Four derivatives use the same central states/tableau. Retain increments
// before adding coordinates, for the shared dense event authority.
public bool TryCoupledRK4(CoupledState start,Real h,float M,float a,float Q,float Lambda,
                          out CoupledState finish,out CoupledIncrement increment,
                          inout uint stageCount) {
    finish=start;increment=(CoupledIncrement)0;
    if(!isfinite(h)||!(h>Real(0))) return false;
    CoupledState s1,s2,s3,s4;
    ++stageCount;if(!TryCoupledDerivative(start,M,a,Q,Lambda,s1)) return false;
    CoupledState second=CoupledAdd(start,s1,h,Real(0.5));
    ++stageCount;if(!TryCoupledDerivative(second,M,a,Q,Lambda,s2)) return false;
    CoupledState third=CoupledAdd(start,s2,h,Real(0.5));
    ++stageCount;if(!TryCoupledDerivative(third,M,a,Q,Lambda,s3)) return false;
    CoupledState fourth=CoupledAdd(start,s3,h,Real(1));
    ++stageCount;if(!TryCoupledDerivative(fourth,M,a,Q,Lambda,s4)) return false;
    // Parentheses follow the central Delta authority, not an alternative tableau.
    Real sixth=Real(1)/Real(6);
    increment.dx=CoupledRetainedProduct((s1.x*h+Real(2)*(s2.x*h))+(Real(2)*(s3.x*h)+s4.x*h),sixth);
    Real4 dk=((s1.k*h+Real(2)*(s2.k*h))+(Real(2)*(s3.k*h)+s4.k*h))*sixth;
    finish.x=CoupledAddRetainedIncrement(start.x,increment.dx);finish.k=start.k+dk;
    for(int i=0;i<4;++i) {
        increment.dX[i]=CoupledRetainedProduct((s1.X[i]*h+Real(2)*(s2.X[i]*h))+
                         (Real(2)*(s3.X[i]*h)+s4.X[i]*h),sixth);
        Real4 dK=((s1.K[i]*h+Real(2)*(s2.K[i]*h))+
                  (Real(2)*(s3.K[i]*h)+s4.K[i]*h))*sixth;
        finish.X[i]=CoupledAddRetainedIncrement(start.X[i],increment.dX[i]);finish.K[i]=start.K[i]+dK;
    }
    return CoupledFinite(finish)&&IsKerrSchildFamilyEventRepresented(
        CoupledVector(finish.x),M,a,Q,Lambda);
}

// Differentiate the exact selected component/root. Other K components remain
// the incoming coordinate variations; no independent projection of V is made.
public bool TryProjectCoupled(CoupledState input,float M,float a,float Q,float Lambda,
    out CoupledState projected,out Real4 covariantV[4],out int selectedComponent) {
    projected=input;selectedComponent=-1;
    for(int i=0;i<4;++i) covariantV[i]=Real4(Real(0));
    if(!CoupledFinite(input)) return false;
    CoupledJet4 metric[4][4];CoupledJet4 connection[4][4][4];
    if(!TryKerrSchildConnectionJet(CoupledVector(input.x),input.X,M,a,Q,Lambda,
                                  metric,connection)) return false;
    Real g[4][4];
    for(int mu=0;mu<4;++mu) for(int nu=0;nu<4;++nu) g[mu][nu]=metric[mu][nu].value;
    Vec4Cart tangent;
    if(!ProjectNullTangentWithSelection(g,CoupledVector(input.k),tangent,selectedComponent))
        return false;
    projected.k=CoupledPack(tangent);
    Real denominator=Real(0),scale=Real(0);
    for(int mu=0;mu<4;++mu) {
        Real term=(g[selectedComponent][mu]+g[mu][selectedComponent])*projected.k[mu];
        denominator+=term;scale+=abs(term);
    }
#ifdef SIRIUS_FP64
    Real epsilon=Real(2.2204460492503131e-16);
#else
    Real epsilon=Real(1.1920928955078125e-7);
#endif
    if(!isfinite(denominator)||!isfinite(scale)||!(scale>Real(0))||
       abs(denominator)<=Real(256)*epsilon*scale) return false;
    for(int i=0;i<4;++i) {
        Real numerator=Real(0);
        for(int mu=0;mu<4;++mu) for(int nu=0;nu<4;++nu) {
            numerator+=metric[mu][nu].derivative[i]*projected.k[mu]*projected.k[nu];
            if(nu!=selectedComponent)
                numerator+=(g[mu][nu]+g[nu][mu])*projected.k[mu]*input.K[i][nu];
        }
        if(!isfinite(numerator)) return false;
        projected.K[i][selectedComponent]=-numerator/denominator;
        covariantV[i]=projected.K[i];
        for(int mu=0;mu<4;++mu) for(int nu=0;nu<4;++nu) for(int rho=0;rho<4;++rho)
            covariantV[i][mu]+=connection[mu][nu][rho].value*projected.k[nu]*projected.X[i][rho];
        if(!all(isfinite(covariantV[i]))) return false;
    }
    return CoupledFinite(projected);
}

public bool TryCoupledEndpointError(CoupledState first,Real4 firstV[4],
    CoupledState second,Real4 secondV[4],CoupledControl control,out Real ratio) {
    ratio=Real(0);
    if(!CoupledControlRepresented(control)||!CoupledFinite(first)||!CoupledFinite(second))
        return false;
    for(int i=0;i<4;++i) for(int mu=0;mu<4;++mu) {
        Real sx=control.variationTolerance*(control.lengthScale+
            max(abs(first.X[i][mu]),abs(second.X[i][mu])));
        Real sv=control.variationTolerance*(control.frequencyScale+
            max(abs(firstV[i][mu]),abs(secondV[i][mu])));
        Real ex=abs(first.X[i][mu]-second.X[i][mu])/sx;
        Real ev=abs(firstV[i][mu]-secondV[i][mu])/sv;
        if(!isfinite(sx)||!isfinite(sv)||!(sx>Real(0))||!(sv>Real(0))||
           !isfinite(ex)||!isfinite(ev)) return false;
        ratio=max(ratio,max(ex,ev));
    }
    // Keep the existing central convention in addition to the unprojected
    // ratio. Projection amplification is compared without a 1/15 assumption.
    for(int mu=0;mu<4;++mu) {
        Real sx=control.centralTolerance*max(Real(1),max(abs(first.x[mu]),abs(second.x[mu])));
        Real sk=control.centralTolerance*max(Real(1),max(abs(first.k[mu]),abs(second.k[mu])));
        Real ex=abs(first.x[mu]-second.x[mu])/sx;
        Real ek=abs(first.k[mu]-second.k[mu])/sk;
        if(!isfinite(sx)||!isfinite(sk)||!(sx>Real(0))||!(sk>Real(0))||
           !isfinite(ex)||!isfinite(ek)) return false;
        ratio=max(ratio,max(ex,ek));
    }
    return isfinite(ratio);
}

public CoupledTrial AttemptCoupledGeodesicStep(CoupledState start,Real step,
    float M,float a,float Q,float Lambda,CoupledControl control) {
    CoupledTrial result=(CoupledTrial)0;
    result.full=start;result.midpoint=start;result.fine=start;
    result.status=kCoupledInvalid;result.stepFactor=Real(0.5);
    if(!CoupledControlRepresented(control)||!CoupledFinite(start)||
       !isfinite(step)||!(step>Real(0))) return result;
    if(!TryCoupledRK4(start,step,M,a,Q,Lambda,result.full,result.fullIncrement,result.stages)||
       !TryCoupledRK4(start,step*Real(0.5),M,a,Q,Lambda,result.midpoint,
                      result.midpointIncrement,result.stages)||
       !TryCoupledRK4(result.midpoint,step*Real(0.5),M,a,Q,Lambda,result.fine,
                      result.fineIncrement,result.stages)) {
        result.status=kCoupledDerivativeDomain;return result;
    }
    result.unprojectedError=RelativeStepDifference(CoupledCentral(start),
        CoupledCentral(result.full),CoupledCentral(result.fine))/control.centralTolerance;
    Real g[4][4],inverse[4][4];
    GetKerrSchildFamilyMetric(CoupledVector(result.fine.x),M,a,Q,Lambda,g,inverse);
    result.unprojectedNull=RelativeNullConstraint(g,CoupledVector(result.fine.k));
    CoupledState projectedFull,projectedFine;
    Real4 fullV[4],fineV[4];int fullComponent,fineComponent;
    if(!TryProjectCoupled(result.full,M,a,Q,Lambda,projectedFull,fullV,fullComponent)||
       !TryProjectCoupled(result.fine,M,a,Q,Lambda,projectedFine,fineV,fineComponent)||
       fullComponent!=fineComponent) {
        result.status=kCoupledProjection;return result;
    }
    result.full=projectedFull;result.fine=projectedFine;
    result.projectedNull=RelativeNullConstraint(g,CoupledVector(result.fine.k));
    if(!TryCoupledEndpointError(result.full,fullV,result.fine,fineV,control,result.errorRatio)||
       !isfinite(result.unprojectedError)||!isfinite(result.unprojectedNull)||
       !isfinite(result.projectedNull)) return result;
    result.errorRatio=max(result.errorRatio,result.unprojectedError);
    float error=max(float(result.errorRatio),1.0e-12f);
    result.stepFactor=Real(clamp(0.9f*pow(1.0f/error,0.2f),0.2f,2.0f));
    if(result.errorRatio>Real(1)||result.unprojectedNull>control.nullTolerance||
       result.projectedNull>control.nullTolerance) {
        result.status=kCoupledAboveBudget;result.stepFactor=min(result.stepFactor,Real(0.5));
        return result;
    }
    result.status=kCoupledNone;return result;
}
