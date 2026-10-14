#pragma once

// Independent test geometry: double chord lengths determine the central angle.
// Production uses a float cross/dot angle. Both consume the same represented
// input vectors, including their permitted departure from unit length.
#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <vector>

namespace sirius::test::point_star_oracle {
using Vector = std::array<double, 3>;
inline double Dot(Vector a, Vector b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline Vector Unit(Vector v) {
    const double length = std::hypot(v[0], v[1], v[2]);
    for (double& x : v) x /= length;
    return v;
}
inline std::array<Vector, 2> Basis(Vector v) {
    const auto at = std::min_element(v.begin(), v.end(),
                                     [](double x, double y) { return std::abs(x) < std::abs(y); });
    const auto axis = static_cast<std::size_t>(at - v.begin());
    Vector a{};
    for (std::size_t k = 0; k < 3; ++k) a[k] = (k == axis ? 1.0 : 0.0) - v[axis] * v[k];
    a = Unit(a);
    return {
        a, Vector{v[1] * a[2] - v[2] * a[1], v[2] * a[0] - v[0] * a[2], v[0] * a[1] - v[1] * a[0]}};
}
struct Case {
    std::array<float, 3> direction;
    std::array<float, 3> star;
    float major;
    float minor;
    float orientation;
};
inline double Weight(const Case& c, bool circular = false) {
    const Vector u = Unit({c.direction[0], c.direction[1], c.direction[2]});
    const Vector v = Unit({c.star[0], c.star[1], c.star[2]});
    Vector difference{}, sum{};
    for (std::size_t k = 0; k < 3; ++k) {
        difference[k] = v[k] - u[k];
        sum[k] = v[k] + u[k];
    }
    const double angle = 2 * std::atan2(std::hypot(difference[0], difference[1], difference[2]),
                                        std::hypot(sum[0], sum[1], sum[2]));
    if (angle > std::min(4.0 * c.major, std::numbers::pi)) return 0;
    if (circular) return std::exp(-0.5 * std::pow(angle / c.major, 2));
    if (angle == 0) return 1;
    const auto basis = Basis(u);
    const double cosine = Dot(u, v);
    Vector tangent{};
    for (std::size_t k = 0; k < 3; ++k) tangent[k] = v[k] - cosine * u[k];
    if (Dot(tangent, tangent) == 0) {
        // The exact antipodal anisotropic log map is outside the represented
        // footprint; equal axes retain the direction-independent circle.
        return c.major == c.minor ? std::exp(-0.5 * std::pow(angle / c.major, 2)) : 0;
    }
    tangent = Unit(tangent);
    const double x = angle * Dot(tangent, basis[0]);
    const double y = angle * Dot(tangent, basis[1]);
    const double a = (std::cos(c.orientation) * x + std::sin(c.orientation) * y) / c.major;
    const double b = (-std::sin(c.orientation) * x + std::cos(c.orientation) * y) / c.minor;
    return std::exp(-0.5 * (a * a + b * b));
}
inline Case RecordedImaxFailure() {
    return {{-0.7856376767158508f, -0.6129410862922668f, 0.0841231718659401f},
            {-0.785654604434967f, -0.612954318523407f, 0.08386841416358948f},
            0.0002556634717620909f,
            0.0002556634717620909f,
            0.0f};
}
inline std::vector<Case> Cases() {
    std::vector<Case> cases{RecordedImaxFailure(),
                            {{1, 0, 0}, {-1, 0, 0}, 1.0f, 1.0f, 0.0f},
                            {{1, 0, 0}, {-1, 0, 0}, 1.0f, 0.5f, 0.0f},
                            {{1e10f, 1e10f, 0}, {1, 0, 0}, 0.1f, 0.1f, 0.0f},
                            {{1e-10f, 1e-10f, 0}, {1, 0, 0}, 0.1f, 0.1f, 0.0f}};
    // Axis, off-axis, longitude seam, pole and latitude-cell boundary.
    const double boundary = 31 * std::numbers::pi / 256;
    const Vector centres[] = {{1, 0, 0},     {0.4, 0.5, 0.7},
                              {-1, 1e-7, 0}, {1e-5, 0, 1},
                              {1e-5, 0, -1}, {std::sin(boundary), 0, std::cos(boundary)}};
    for (Vector centre : centres) {
        centre = Unit(centre);
        const auto basis = Basis(centre);
        for (int height : {1080, 4096}) {
            for (float pixels : {0.3f, 1.0f}) {
                const float sigma = static_cast<float>(std::numbers::pi / 3 / height) * pixels;
                for (float ratio : {1.0f, 4.0f}) {
                    for (float orientation : {0.0f, 0.7f}) {
                        for (double position_angle :
                             {0.0, std::numbers::pi / 4, std::numbers::pi / 2}) {
                            for (double offset : {0.0, 0.25, 0.5, 1.0, 2.0, 3.99, 4.01}) {
                                Case c{};
                                c.major = sigma * ratio;
                                c.minor = sigma;
                                c.orientation = orientation;
                                const double angle = offset * c.major;
                                for (std::size_t k = 0; k < 3; ++k) {
                                    c.direction[k] = static_cast<float>(centre[k]);
                                    c.star[k] = static_cast<float>(
                                        std::cos(angle) * centre[k] +
                                        std::sin(angle) *
                                            (std::cos(orientation + position_angle) * basis[0][k] +
                                             std::sin(orientation + position_angle) * basis[1][k]));
                                }
                                cases.push_back(c);
                            }
                        }
                    }
                }
            }
        }
    }
    return cases;
}
}  // namespace sirius::test::point_star_oracle
