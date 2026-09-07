#pragma once

#include "sirius/core/metrics/kerr_schild_family.h"

#include <array>
#include <optional>

namespace sirius::core {

// Past-directed rays cross the past horizon in the outgoing chart. The
// ingoing chart used by the public metric covers the future horizon instead.
// See Bozzola, Chan & Paschalidis, Phys. Rev. D 108, 084004 (2023).
// A time/azimuth reversal pulls back the existing analytic metric and its
// derivatives; the extra azimuth reversal retains the physical spin sign.
class OutgoingKerrSchild final : public IMetric {
  public:
    struct ChartMap {
        Vec4 position;
        std::array<std::array<double, 4>, 4> jacobian{};

        [[nodiscard]] Vec4 Apply(const Vec4& vector) const {
            Vec4 result;
            for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu) result(mu) += jacobian[mu][nu] * vector(nu);
            return result;
        }
    };

    explicit OutgoingKerrSchild(KerrSchildFamily& source) : source_(source) {}

    void Evaluate(const Vec4& position, Metric4d& metric,
                  Tensor<Dual<double>, 4, 4, 4>& derivative) override {
        source_.Evaluate(Reflect(position), metric, derivative);
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                metric(mu, nu) *= kReflection[mu] * kReflection[nu];
                for (int rho = 0; rho < 4; ++rho)
                    derivative(mu, nu, rho) *= kReflection[mu] * kReflection[nu] * kReflection[rho];
            }
        }
    }

    bool InverseMetric(const Vec4& position, Metric4d& inverse) const override {
        if (!source_.InverseMetric(Reflect(position), inverse)) return false;
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) inverse(mu, nu) *= kReflection[mu] * kReflection[nu];
        return true;
    }

    bool IsValidEvent(const Vec4& position) const override {
        return source_.IsValidEvent(Reflect(position));
    }
    bool InsideCaptureSurface(const Vec4& position, double margin) const override {
        return source_.InsideCaptureSurface(Reflect(position), margin);
    }
    const Config& GetParameters() const override { return source_.GetParameters(); }
    void SetParameter(const std::string& key, double value) override {
        source_.SetParameter(key, value);
    }
    const char* GetName() const override { return source_.GetName(); }
    KerrSchildFamily& Source() const { return source_; }

    [[nodiscard]] std::optional<ChartMap> FromIngoing(const Vec4& position) const {
        return Map(position, 1.0);
    }
    [[nodiscard]] std::optional<ChartMap> ToIngoing(const Vec4& position) const {
        return Map(position, -1.0);
    }

    // Only azimuth-dependent matter needs this scalar map. Circular emitter
    // frequencies depend on the unchanged stationary/axial Killing quantities.
    [[nodiscard]] double IngoingAzimuth(double radius, double azimuth) const {
        const auto shift = RadialShift(radius);
        return shift ? azimuth - (*shift)[1] : std::numeric_limits<double>::quiet_NaN();
    }

  private:
    KerrSchildFamily& source_;
    static constexpr std::array<double, 4> kReflection{-1.0, 1.0, -1.0, 1.0};

    static Vec4 Reflect(Vec4 position) {
        for (int mu = 0; mu < 4; ++mu) position(mu) *= kReflection[mu];
        return position;
    }

    // Integral of D=(1-f)/f. Additive constants in time and azimuth are fixed
    // at an exterior reference event and have no physical significance.
    [[nodiscard]] std::optional<std::array<double, 2>> RadialShift(double radius) const {
        if (!source_.HasHorizon()) return std::nullopt;
        const auto p = source_.GetParams();
        const double horizon = source_.OuterHorizonRadius();
        if (!(radius > horizon)) return std::nullopt;
        double reference = 2.0 * horizon;
        double time_integral = 0.0;
        double angle_integral = 0.0;
        if (p.Lambda > 0.0) {
            const double cosmological = source_.CosmologicalHorizonRadius();
            if (!(radius < cosmological) || !(horizon < cosmological)) return std::nullopt;
            reference = 0.5 * (horizon + cosmological);
            time_integral = -(radius - reference);
            for (const double root : {horizon, cosmological, -horizon - cosmological}) {
                const double slope = 2.0 * p.M / (root * root) - 2.0 * p.Lambda * root / 3.0;
                time_integral += std::log(std::abs((radius - root) / (reference - root))) / slope;
            }
        } else {
            const double inner = source_.InnerHorizonRadius();
            const double separation = horizon - inner;
            double inverse_delta_integral;
            if (separation == 0.0) {
                const double inverse_difference =
                    1.0 / (radius - horizon) - 1.0 / (reference - horizon);
                time_integral = 2.0 * p.M * std::log((radius - horizon) / (reference - horizon)) -
                                (p.M * p.M + p.a * p.a) * inverse_difference;
                inverse_delta_integral = -inverse_difference;
            } else {
                const double outer_log = std::log((radius - horizon) / (reference - horizon));
                const double inner_log = std::log((radius - inner) / (reference - inner));
                time_integral = ((2.0 * p.M * horizon - p.Q * p.Q) * outer_log -
                                 (2.0 * p.M * inner - p.Q * p.Q) * inner_log) /
                                separation;
                inverse_delta_integral = (outer_log - inner_log) / separation;
            }
            if (p.a != 0.0) {
                angle_integral = -2.0 * p.a * inverse_delta_integral +
                                 2.0 * (std::atan2(radius, p.a) - std::atan2(reference, p.a));
            }
        }
        if (!std::isfinite(time_integral) || !std::isfinite(angle_integral)) return std::nullopt;
        return std::array<double, 2>{-2.0 * time_integral, angle_integral};
    }

    [[nodiscard]] std::optional<ChartMap> Map(const Vec4& position, double orientation) const {
        const auto p = source_.GetParams();
        const coordinates::Vec4Cart cart{position(0), position(1), position(2), position(3)};
        const auto geometry = coordinates::TryKerrSchildRadiusDifferential(cart, p.a);
        if (!geometry) return std::nullopt;
        const double r = geometry->radius;
        const auto shift = RadialShift(r);
        if (!shift) return std::nullopt;
        const double delta =
            r * r - 2.0 * p.M * r + p.a * p.a + p.Q * p.Q - p.Lambda * r * r * r * r / 3.0;
        const double d = (r * r + p.a * p.a) / delta - 1.0;
        const double time_derivative = -2.0 * orientation * d;
        const double angle_derivative = -2.0 * orientation * p.a * d / (r * r + p.a * p.a);
        const double angle = orientation * (*shift)[1];
        const double cosine = std::cos(angle);
        const double sine = std::sin(angle);
        ChartMap result;
        result.position = position;
        result.position(0) += orientation * (*shift)[0];
        result.position(1) = cosine * position(1) - sine * position(2);
        result.position(2) = sine * position(1) + cosine * position(2);
        const std::array<double, 4> radial{0.0, geometry->dx, geometry->dy, geometry->dz};
        result.jacobian[0][0] = 1.0;
        result.jacobian[1][1] = cosine;
        result.jacobian[1][2] = -sine;
        result.jacobian[2][1] = sine;
        result.jacobian[2][2] = cosine;
        result.jacobian[3][3] = 1.0;
        for (int nu = 1; nu < 4; ++nu) {
            result.jacobian[0][nu] += time_derivative * radial[nu];
            result.jacobian[1][nu] -= result.position(2) * angle_derivative * radial[nu];
            result.jacobian[2][nu] += result.position(1) * angle_derivative * radial[nu];
        }
        for (int mu = 0; mu < 4; ++mu) {
            if (!std::isfinite(result.position(mu))) return std::nullopt;
            for (int nu = 0; nu < 4; ++nu)
                if (!std::isfinite(result.jacobian[mu][nu])) return std::nullopt;
        }
        return result;
    }
};

}  // namespace sirius::core
