#pragma once

#include <cmath>
#include <optional>
#include <iostream>
#include <iomanip>

#include "region.hpp"
#include "ring.hpp"
#include "mymath.hpp"

namespace gridsynth
{

    /**
     * C++ implementation of EpsilonRegion from Python gridsynth project.
     * This is a concrete implementation of ConvexSet used for testing TDGP.
     */
    class EpsilonRegion : public ConvexSet
    {
    private:
        Float _theta;
        Float _epsilon;
        Float _d;
        Float _z_x;
        Float _z_y;
        Ellipse _ellipse;

        static Ellipse create_proper_ellipse(Float theta, Float epsilon)
        {
            Float d = Float(1.0) - epsilon * epsilon / Float(2.0);
            Float z_x = cos(-theta / Float(2.0));
            Float z_y = sin(-theta / Float(2.0));

            // Matrix multiplication D_1 * D_2 * D_3 from Python version
            // D_1 = [[z_x, -z_y], [z_y, z_x]]
            // D_2 = [[4 * (1/epsilon)^4, 0], [0, (1/epsilon)^2]]
            // D_3 = [[z_x, z_y], [-z_y, z_x]]

            Float inv_eps2 = Float(1.0) / (epsilon * epsilon);
            Float inv_eps4 = inv_eps2 * inv_eps2;

            // D_2 values
            Float d2_00 = Float(4.0) * inv_eps4;
            Float d2_11 = inv_eps2;

            // Final matrix: D_1 * D_2 * D_3
            // This is a complex 2x2 matrix multiplication
            Float final_a = d2_00 * z_x * z_x + d2_11 * z_y * z_y;
            Float final_b = (d2_00 - d2_11) * z_x * z_y;
            Float final_d = d2_00 * z_y * z_y + d2_11 * z_x * z_x;

            // Center point p
            Float center_x = d * z_x;
            Float center_y = d * z_y;

            return Ellipse(final_a, final_b, final_d, center_x, center_y);
        }

    public:
        EpsilonRegion(Float theta, Float epsilon)
            : _theta(theta), _epsilon(epsilon), _d(Float(1.0) - epsilon * epsilon / Float(2.0)),
              _z_x(cos(-theta / Float(2.0))), _z_y(sin(-theta / Float(2.0))),
              _ellipse(create_proper_ellipse(theta, epsilon))
        {
        }

        // Properties
        Float theta() const { return _theta; }
        Float epsilon() const { return _epsilon; }
        const Ellipse &ellipse() const { return _ellipse; }

        bool inside(const std::array<Float, 2> &u) const override
        {
            Float real_part = u[0];
            Float imag_part = u[1];
            Float norm_squared = real_part * real_part + imag_part * imag_part;
            Float tol = Float(1e-30);
            if (norm_squared > Float(1.0) + tol)
                return false;
            Float cos_similarity = _z_x * real_part + _z_y * imag_part;
            return cos_similarity + tol >= _d; // allow tiny undershoot due to rounding
        }

        std::optional<std::pair<Float, Float>> intersect(
            const std::array<Float, 2> &u0,
            const std::array<Float, 2> &v) const override
        {
            // Line parameterized as u0 + t*v
            // We need to intersect with |u|^2 <= 1 and cos_similarity >= _d

            // For |u0 + t*v|^2 <= 1:
            // |u0|^2 + 2*t*Re(u0*conj(v)) + t^2*|v|^2 <= 1
            Float a = v[0] * v[0] + v[1] * v[1];
            Float b = Float(2.0) * (u0[0] * v[0] + u0[1] * v[1]);
            Float c = u0[0] * u0[0] + u0[1] * u0[1] - Float(1.0);

            auto quad_solution = solve_quadratic_optional(a, b, c);
            if (!quad_solution.has_value())
                return std::nullopt;

            Float t0 = quad_solution->first;
            Float t1 = quad_solution->second;

            // For cos similarity constraint: _z_x * (u0[0] + t*v[0]) + _z_y * (u0[1] + t*v[1]) >= _d
            // _z_x * u0[0] + _z_y * u0[1] + t * (_z_x * v[0] + _z_y * v[1]) >= _d
            Float vz = _z_x * v[0] + _z_y * v[1];
            Float rhs = _d - _z_x * u0[0] - _z_y * u0[1];
            Float tolerance = Float(1e-30);

            if (vz > tolerance)
            {
                Float t2 = rhs / vz;
                Float t_start = max(t0, t2);
                Float t_end = t1;
                if (t_start > t_end)
                    return std::nullopt;
                return std::make_pair(t_start, t_end);
            }
            else if (vz < -tolerance)
            {
                Float t2 = rhs / vz;
                Float t_start = t0;
                Float t_end = min(t1, t2);
                if (t_start > t_end)
                    return std::nullopt;
                return std::make_pair(t_start, t_end);
            }
            else
            {
                if (rhs <= tolerance)
                    return std::make_pair(t0, t1);
                return std::nullopt;
            }
        }

        bool intersects_line_segment(const std::array<Float, 2> &p, const std::array<Float, 2> &q) const override
        {
            // Simple implementation: check endpoints and midpoint
            if (inside(p) || inside(q))
                return true;

            std::array<Float, 2> midpoint = {(p[0] + q[0]) / Float(2.0), (p[1] + q[1]) / Float(2.0)};
            return inside(midpoint);
        }

        std::unique_ptr<ConvexSet> clone() const override
        {
            return std::make_unique<EpsilonRegion>(*this);
        }
    };

} // namespace gridsynth
