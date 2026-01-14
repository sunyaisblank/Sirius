// PHGD002A.cpp - Numerical Geodesic Integrator Implementation
#include "PHGD002A.h"
#include <PHMT000A.h>
#include <MTTN001A.h>
#include <cmath>
#include <limits>

namespace Sirius {
namespace Physics {

NumericalGeodesicIntegrator::NumericalGeodesicIntegrator() = default;
NumericalGeodesicIntegrator::~NumericalGeodesicIntegrator() = default;

void NumericalGeodesicIntegrator::setMetric(IMetric* metric) {
    m_Metric = metric;
}

double NumericalGeodesicIntegrator::getRadius(const GeodesicState& state) const {
    // For spherical coordinates: x = (t, r, theta, phi)
    return state.x(1);
}

bool NumericalGeodesicIntegrator::isValidState(const GeodesicState& state) const {
    for (int i = 0; i < 4; ++i) {
        if (std::isnan(state.x(i)) || std::isinf(state.x(i))) return false;
        if (std::isnan(state.u(i)) || std::isinf(state.u(i))) return false;
    }
    return true;
}

GeodesicState NumericalGeodesicIntegrator::computeDerivative(const GeodesicState& state) {
    GeodesicState deriv;
    
    // dx^μ/dλ = u^μ
    deriv.x = state.u;
    
    // Get metric and its derivatives at current position
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    
    // Convert state position to Tensor format
    Tensor<double, 4> pos;
    for (int i = 0; i < 4; ++i) {
        pos(i) = state.x(i);
    }
    
    m_Metric->evaluate(pos, g, dg);
    
    // Use TensorOps to compute Christoffel symbols and geodesic acceleration
    ChristoffelSymbols gamma = TensorOps::christoffel(g, dg);
    Vec4 acceleration = TensorOps::geodesicAcceleration(state.u, gamma);
    
    deriv.u = acceleration;
    
    return deriv;
}

GeodesicState NumericalGeodesicIntegrator::rk4Step(const GeodesicState& state, double h) {
    // k1 = f(y_n)
    GeodesicState k1 = computeDerivative(state);
    
    // k2 = f(y_n + h/2 * k1)
    GeodesicState y2;
    for (int i = 0; i < 4; ++i) {
        y2.x(i) = state.x(i) + 0.5 * h * k1.x(i);
        y2.u(i) = state.u(i) + 0.5 * h * k1.u(i);
    }
    GeodesicState k2 = computeDerivative(y2);
    
    // k3 = f(y_n + h/2 * k2)
    GeodesicState y3;
    for (int i = 0; i < 4; ++i) {
        y3.x(i) = state.x(i) + 0.5 * h * k2.x(i);
        y3.u(i) = state.u(i) + 0.5 * h * k2.u(i);
    }
    GeodesicState k3 = computeDerivative(y3);
    
    // k4 = f(y_n + h * k3)
    GeodesicState y4;
    for (int i = 0; i < 4; ++i) {
        y4.x(i) = state.x(i) + h * k3.x(i);
        y4.u(i) = state.u(i) + h * k3.u(i);
    }
    GeodesicState k4 = computeDerivative(y4);
    
    // y_{n+1} = y_n + h/6 * (k1 + 2*k2 + 2*k3 + k4)
    GeodesicState next;
    for (int i = 0; i < 4; ++i) {
        next.x(i) = state.x(i) + (h / 6.0) * (k1.x(i) + 2.0*k2.x(i) + 2.0*k3.x(i) + k4.x(i));
        next.u(i) = state.u(i) + (h / 6.0) * (k1.u(i) + 2.0*k2.u(i) + 2.0*k3.u(i) + k4.u(i));
    }
    
    return next;
}

IntegrationResult NumericalGeodesicIntegrator::integrate(const GeodesicState& initial) {
    IntegrationResult result;
    result.finalState = initial;
    result.affineDist = 0.0;
    result.steps = 0;
    result.hitHorizon = false;
    result.escaped = false;
    result.error = false;
    
    if (!m_Metric) {
        result.error = true;
        return result;
    }
    
    GeodesicState current = initial;
    double h = m_StepSize;
    
    for (int step = 0; step < m_MaxSteps; ++step) {
        result.steps = step + 1;
        
        // Check termination conditions
        double r = getRadius(current);
        
        // Horizon crossing
        if (r < m_HorizonRadius) {
            result.hitHorizon = true;
            result.finalState = current;
            return result;
        }
        
        // Escape to infinity
        if (r > m_MaxDistance) {
            result.escaped = true;
            result.finalState = current;
            return result;
        }
        
        // Adaptive step size based on distance to horizon
        double distToHorizon = r - m_HorizonRadius;
        if (distToHorizon > 0.0 && distToHorizon < 1.0) {
            h = m_StepSize * distToHorizon;
            if (h < m_StepSize * 0.01) h = m_StepSize * 0.01;
        } else {
            h = m_StepSize;
        }
        
        // RK4 integration step
        current = rk4Step(current, h);
        result.affineDist += h;
        
        // Numerical stability check
        if (!isValidState(current)) {
            result.error = true;
            return result;
        }
        
        // Keep theta in valid range for spherical coords
        if (current.x(2) < 0.01) current.x(2) = 0.01;
        if (current.x(2) > M_PI - 0.01) current.x(2) = M_PI - 0.01;
        
        // Wrap phi
        while (current.x(3) < 0.0) current.x(3) += 2.0 * M_PI;
        while (current.x(3) >= 2.0 * M_PI) current.x(3) -= 2.0 * M_PI;
    }
    
    result.finalState = current;
    return result;
}

} // namespace Physics
} // namespace Sirius
