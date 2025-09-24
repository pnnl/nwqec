#pragma once

#include <array>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <memory>
#include <utility>
#include <optional>

// Include GridOp after class declarations to avoid circular dependency
#include "grid_op.hpp"

namespace gridsynth
{

    // Forward declaration for GridOp to avoid circular dependency
    class GridOp;

    // Abstract base class for convex sets
    class ConvexSet
    {
    public:
        virtual ~ConvexSet() = default;
        virtual bool inside(const std::array<Float, 2> &v) const = 0;
        virtual bool intersects_line_segment(const std::array<Float, 2> &p, const std::array<Float, 2> &q) const = 0;
        virtual std::unique_ptr<ConvexSet> clone() const = 0;

        // Method for line intersection - returns parameter values where line intersects set boundary
        virtual std::optional<std::pair<Float, Float>> intersect(const std::array<Float, 2> &u0, const std::array<Float, 2> &v) const = 0;
    };

    class Interval
    {
    private:
        Float _a, _b;

    public:
        Interval(Float a, Float b) : _a(a), _b(b)
        {
            if (a > b)
            {
                throw std::invalid_argument("Interval: a must be <= b");
            }
        }

        inline const Float &a() const noexcept { return _a; }
        inline const Float &b() const noexcept { return _b; }

        // Aliases for different naming conventions
        inline const Float &l() const noexcept { return _a; }          // left endpoint
        inline const Float &r() const noexcept { return _b; }          // right endpoint
        inline Float width() const noexcept { return _b - _a; } // alias for length

        inline Float length() const noexcept { return _b - _a; }
        inline Float center() const noexcept { return (_a + _b) / Float(2.0); }

        inline bool inside(const Float &x) const noexcept { return x >= _a && x <= _b; }
        inline bool within(const Float &x) const noexcept { return inside(x); } // alias for inside

        // Create a fattened (expanded) interval
        inline Interval fatten(const Float &amount) const
        {
            return Interval(_a - amount, _b + amount);
        }

        // Arithmetic operators
        inline Interval operator*(const Float &scale) const
        {
            if (scale >= 0)
            {
                return Interval(_a * scale, _b * scale);
            }
            else
            {
                return Interval(_b * scale, _a * scale); // reverse order for negative scale
            }
        }

        inline Interval operator+(const Float &shift) const
        {
            return Interval(_a + shift, _b + shift);
        }

        inline Interval operator-(const Float &shift) const
        {
            return Interval(_a - shift, _b - shift);
        }

        std::vector<Float> tolist() const
        {
            return {_a, _b};
        }

        // For printing/debugging
        std::string to_string() const
        {
            return "[" + _a.to_string() + ", " + _b.to_string() + "]";
        }
    };

    class Rectangle : public ConvexSet
    {
    private:
        Interval _x_range;
        Interval _y_range;

    public:
        Rectangle(Float x1, Float x2, Float y1, Float y2)
            : _x_range(x1, x2), _y_range(y1, y2) {}

        Rectangle(const Interval &x_range, const Interval &y_range)
            : _x_range(x_range), _y_range(y_range) {}

        const Interval &x_range() const { return _x_range; }
        const Interval &y_range() const { return _y_range; }

        // Aliases used in TDGP
        const Interval &I_x() const { return _x_range; }
        const Interval &I_y() const { return _y_range; }

        Float width() const { return _x_range.length(); }
        Float height() const { return _y_range.length(); }
        Float area() const { return width() * height(); }

        bool inside(const std::array<Float, 2> &v) const override
        {
            return _x_range.inside(v[0]) && _y_range.inside(v[1]);
        }

        bool intersects_line_segment(const std::array<Float, 2> &p, const std::array<Float, 2> &q) const override
        {
            // Use Liang-Barsky line clipping algorithm
            Float dx = q[0] - p[0];
            Float dy = q[1] - p[1];

            // Use precision-adaptive tolerance instead of hardcoded 1e-12
            Float tolerance = Float(1e-15); // More conservative for high precision
            if (abs(dx) < tolerance && abs(dy) < tolerance)
            {
                return inside(p);
            }

            Float t_min = 0.0, t_max = 1.0;

            // Check against each boundary
            std::array<Float, 4> p_vals = {-dx, dx, -dy, dy};
            std::array<Float, 4> q_vals = {p[0] - _x_range.a(), _x_range.b() - p[0],
                                           p[1] - _y_range.a(), _y_range.b() - p[1]};

            for (Integer i = 0; i < 4; ++i)
            {
                Float tolerance = Float(1e-15); // Use high precision tolerance
                if (abs(p_vals[i]) < tolerance)
                {
                    if (q_vals[i] < 0)
                        return false;
                }
                else
                {
                    Float t = q_vals[i] / p_vals[i];
                    if (p_vals[i] < 0)
                    {
                        t_max = min(t_max, t);
                    }
                    else
                    {
                        t_min = max(t_min, t);
                    }
                    if (t_min > t_max)
                        return false;
                }
            }

            return true;
        }

        std::optional<std::pair<Float, Float>> intersect(const std::array<Float, 2> &u0, const std::array<Float, 2> &v) const override
        {
            // Return full (possibly unbounded) parameter interval along infinite line:
            // { t | u0 + t v inside rectangle }. This matches Python TDGP expectation.
            // Compute intersection against each slab independently.
            Float tol = Float(1e-30); // tighter tolerance for high precision scenarios
            Float t_low = -Float(std::numeric_limits<double>::infinity());
            Float t_high = Float(std::numeric_limits<double>::infinity());

            auto update_bounds = [&](Float p0, Float dp, Float a, Float b) -> bool
            {
                if (abs(dp) < tol)
                {
                    // Line parallel to slab; reject if outside
                    if (p0 < a - tol || p0 > b + tol)
                        return false;
                    return true; // No constraint added
                }
                Float t1 = (a - p0) / dp;
                Float t2 = (b - p0) / dp;
                if (t1 > t2)
                    std::swap(t1, t2);
                if (t1 > t_low)
                    t_low = t1;
                if (t2 < t_high)
                    t_high = t2;
                return t_low <= t_high + tol;
            };

            if (!update_bounds(u0[0], v[0], _x_range.a(), _x_range.b()))
                return std::nullopt;
            if (!update_bounds(u0[1], v[1], _y_range.a(), _y_range.b()))
                return std::nullopt;

            if (t_low > t_high)
                return std::nullopt;
            return std::make_pair(t_low, t_high);
        }

        std::unique_ptr<ConvexSet> clone() const override
        {
            return std::make_unique<Rectangle>(*this);
        }

        std::string to_string() const
        {
            return "Rectangle(" + _x_range.to_string() + ", " + _y_range.to_string() + ")";
        }
    };

    class Ellipse : public ConvexSet
    {
    private:
        std::array<std::array<Float, 2>, 2> _D; // 2x2 matrix
        std::array<Float, 2> _p;                // center point [px, py]

    public:
        Ellipse(const std::array<std::array<Float, 2>, 2> &D, const std::array<Float, 2> &p)
            : _D(D), _p(p)
        {
            // Validate that the matrix represents a valid ellipse (positive definite)
            Float det = _D[0][0] * _D[1][1] - _D[0][1] * _D[1][0];
            if (det <= 0 || _D[0][0] <= 0 || _D[1][1] <= 0)
            {
                throw std::domain_error("Ellipse matrix must be positive definite");
            }
        }

        Ellipse(Float a, Float b, Float d, Float px, Float py)
            : _D{{{{a, b}}, {{b, d}}}}, _p{{px, py}}
        {
            // Validate that the matrix represents a valid ellipse (positive definite)
            Float det = a * d - b * b;
            if (det <= 0 || a <= 0 || d <= 0)
            {
                throw std::domain_error("Ellipse parameters must form a positive definite matrix");
            }
        }

        // Properties
        const std::array<std::array<Float, 2>, 2> &D() const { return _D; }
        const std::array<Float, 2> &p() const { return _p; }

        Float px() const { return _p[0]; }
        void set_px(Float px) { _p[0] = px; }

        Float py() const { return _p[1]; }
        void set_py(Float py) { _p[1] = py; }

        // Set both position components at once (faster / fewer calls)
        void set_p(Float px, Float py)
        {
            _p[0] = px;
            _p[1] = py;
        }

        Float a() const { return _D[0][0]; }
        void set_a(Float a) { _D[0][0] = a; }
        void scale_a(Float factor)
        {
            _D[0][0] *= factor;
        }

        Float b() const { return _D[0][1]; }
        void set_b(Float b)
        {
            _D[0][1] = b;
            _D[1][0] = b;
        }

        void flip_b()
        {
            _D[0][1] = -_D[0][1];
            _D[1][0] = -_D[1][0];
        }

        Float d() const { return _D[1][1]; }
        void set_d(Float d) { _D[1][1] = d; }
        void scale_d(Float factor)
        {
            _D[1][1] *= factor;
        }

        // Set D matrix components (a,b,d) in one call
        void set_D(Float a, Float b, Float d)
        {
            _D[0][0] = a;
            _D[0][1] = b;
            _D[1][0] = b;
            _D[1][1] = d;
        }

        // Transform this ellipse in-place by a GridOp (used by EllipsePair)
        // use_conj: true when g_local is the conjugated GridOp for the second ellipse
        // use_fallback: if true, prefer calling g_local.inv() to obtain inverse
        inline void transform_by_gridop(const GridOp &g_local,
                                        bool use_conj,
                                        bool use_fallback,
                                        Float preinv00, Float preinv01, Float preinv10, Float preinv11,
                                        Float tol)
        {
            // Extract forward matrix (real parts)
            auto mat = g_local.to_mat();
            Float F00 = mat[0][0].real();
            Float F01 = mat[0][1].real();
            Float F10 = mat[1][0].real();
            Float F11 = mat[1][1].real();

            // Compute inverse entries I00..I11
            Float I00, I01, I10, I11;
            if (use_fallback)
            {
                GridOp g_inv = g_local.inv();
                auto invm = g_inv.to_mat();
                I00 = invm[0][0].real();
                I01 = invm[0][1].real();
                I10 = invm[1][0].real();
                I11 = invm[1][1].real();
            }
            else
            {
                if (!use_conj)
                {
                    I00 = preinv00;
                    I01 = preinv01;
                    I10 = preinv10;
                    I11 = preinv11;
                }
                else
                {
                    Float d = F00 * F11 - F01 * F10;
                    if (abs(d) < tol)
                    {
                        GridOp g_inv = g_local.inv();
                        auto invm = g_inv.to_mat();
                        I00 = invm[0][0].real();
                        I01 = invm[0][1].real();
                        I10 = invm[1][0].real();
                        I11 = invm[1][1].real();
                    }
                    else
                    {
                        I00 = F11 / d;
                        I01 = -F01 / d;
                        I10 = -F10 / d;
                        I11 = F00 / d;
                    }
                }
            }

            // Compute D' = I^T * D * I where D = [A B; B C]
            Float A = a();
            Float B = b();
            Float C = d();

            Float S00 = I00 * A + I10 * B;
            Float S01 = I00 * B + I10 * C;
            Float S10 = I01 * A + I11 * B;
            Float S11 = I01 * B + I11 * C;

            Float na = S00 * I00 + S01 * I10;
            Float nb = S00 * I01 + S01 * I11;
            Float nd = S10 * I01 + S11 * I11;

            // Preserve old position before updating D
            Float old_px = px();
            Float old_py = py();

            // Update D and position
            set_D(na, nb, nd);
            Float px_new = F00 * old_px + F01 * old_py;
            Float py_new = F10 * old_px + F11 * old_py;
            set_p(px_new, py_new);
        }

        // Variant that takes the forward matrix entries directly to avoid recomputing
        // the underlying GridOp->to_mat() for conjugated GridOps.
        inline void transform_by_gridop_mat(Float F00, Float F01, Float F10, Float F11,
                                            bool use_conj,
                                            bool use_fallback,
                                            Float preinv00, Float preinv01, Float preinv10, Float preinv11,
                                            Float tol,
                                            const GridOp *fallback_g = nullptr)
        {
            // Compute inverse entries I00..I11
            Float I00, I01, I10, I11;
            if (use_fallback)
            {
                if (!fallback_g)
                    throw std::invalid_argument("transform_by_gridop_mat: fallback requested but no GridOp provided");
                GridOp g_inv = fallback_g->inv();
                auto invm = g_inv.to_mat();
                I00 = invm[0][0].real();
                I01 = invm[0][1].real();
                I10 = invm[1][0].real();
                I11 = invm[1][1].real();
            }
            else
            {
                if (!use_conj)
                {
                    I00 = preinv00;
                    I01 = preinv01;
                    I10 = preinv10;
                    I11 = preinv11;
                }
                else
                {
                    Float d = F00 * F11 - F01 * F10;
                    if (abs(d) < tol)
                    {
                        if (!fallback_g)
                            throw std::domain_error("transform_by_gridop_mat: singular matrix and no fallback GridOp");
                        GridOp g_inv = fallback_g->inv();
                        auto invm = g_inv.to_mat();
                        I00 = invm[0][0].real();
                        I01 = invm[0][1].real();
                        I10 = invm[1][0].real();
                        I11 = invm[1][1].real();
                    }
                    else
                    {
                        I00 = F11 / d;
                        I01 = -F01 / d;
                        I10 = -F10 / d;
                        I11 = F00 / d;
                    }
                }
            }

            // Compute D' = I^T * D * I where D = [A B; B C]
            Float A = a();
            Float B = b();
            Float C = d();

            Float S00 = I00 * A + I10 * B;
            Float S01 = I00 * B + I10 * C;
            Float S10 = I01 * A + I11 * B;
            Float S11 = I01 * B + I11 * C;

            Float na = S00 * I00 + S01 * I10;
            Float nb = S00 * I01 + S01 * I11;
            Float nd = S10 * I01 + S11 * I11;

            // Preserve old position before updating D
            Float old_px = px();
            Float old_py = py();

            // Update D and position
            set_D(na, nb, nd);
            Float px_new = F00 * old_px + F01 * old_py;
            Float py_new = F10 * old_px + F11 * old_py;
            set_p(px_new, py_new);
        }

        Float sqrt_det() const
        {
            Float det = d() * a() - b() * b();
            if (det <= 0)
            {
                throw std::domain_error("Ellipse determinant must be positive");
            }
            return sqrt(det);
        }

        Float area() const
        {
            return PI / sqrt_det();
        }

        Float skew() const
        {
            return b() * b();
        }

        Float bias() const
        {
            return d() / a();
        }

        // Methods
        bool inside(const std::array<Float, 2> &v) const override
        {
            Float x = v[0] - px();
            Float y = v[1] - py();
            Float tmp = a() * x * x + 2 * b() * x * y + d() * y * y;
            return tmp <= 1.0;
        }

        Rectangle bbox() const
        {
            Float sqrt_det_val = sqrt_det();
            Float w = sqrt(d()) / sqrt_det_val;
            Float h = sqrt(a()) / sqrt_det_val;
            return Rectangle(px() - w, px() + w, py() - h, py() + h);
        }

        Ellipse normalize() const
        {
            Float sqrt_det_val = sqrt_det();
            std::array<std::array<Float, 2>, 2> new_D;
            new_D[0][0] = _D[0][0] / sqrt_det_val;
            new_D[0][1] = _D[0][1] / sqrt_det_val;
            new_D[1][0] = _D[1][0] / sqrt_det_val;
            new_D[1][1] = _D[1][1] / sqrt_det_val;
            return Ellipse(new_D, _p);
        }

        bool intersects_line_segment(const std::array<Float, 2> &p_start, const std::array<Float, 2> &p_end) const override
        {
            // Translate so ellipse is centered at origin
            Float x1 = p_start[0] - px();
            Float y1 = p_start[1] - py();
            Float x2 = p_end[0] - px();
            Float y2 = p_end[1] - py();

            // Check if either endpoint is inside the ellipse
            if (a() * x1 * x1 + 2 * b() * x1 * y1 + d() * y1 * y1 <= 1.0)
                return true;
            if (a() * x2 * x2 + 2 * b() * x2 * y2 + d() * y2 * y2 <= 1.0)
                return true;

            // Check if line intersects ellipse boundary
            Float dx = x2 - x1;
            Float dy = y2 - y1;

            // Quadratic equation coefficients for line-ellipse intersection
            Float A = a() * dx * dx + Float(2.0) * b() * dx * dy + d() * dy * dy;
            Float B = Float(2.0) * (a() * x1 * dx + b() * (x1 * dy + y1 * dx) + d() * y1 * dy);
            Float C = a() * x1 * x1 + Float(2.0) * b() * x1 * y1 + d() * y1 * y1 - Float(1.0);

            Float discriminant = B * B - Float(4.0) * A * C;
            if (discriminant < Float(0.0))
                return false;

            // Check if intersection points are within the line segment
            Float sqrt_disc = sqrt(discriminant);
            Float t1 = (-B - sqrt_disc) / (Float(2.0) * A);
            Float t2 = (-B + sqrt_disc) / (Float(2.0) * A);

            return (t1 >= 0 && t1 <= 1) || (t2 >= 0 && t2 <= 1);
        }

        std::optional<std::pair<Float, Float>> intersect(const std::array<Float, 2> &u0, const std::array<Float, 2> &v) const override
        {
            // Line intersection with ellipse
            // Translate so ellipse is centered at origin
            Float x0 = u0[0] - px();
            Float y0 = u0[1] - py();
            Float dx = v[0];
            Float dy = v[1];

            // Quadratic equation coefficients for line-ellipse intersection
            Float A = a() * dx * dx + Float(2.0) * b() * dx * dy + d() * dy * dy;
            Float B = Float(2.0) * (a() * x0 * dx + b() * (x0 * dy + y0 * dx) + d() * y0 * dy);
            Float C = a() * x0 * x0 + Float(2.0) * b() * x0 * y0 + d() * y0 * y0 - Float(1.0);

            Float discriminant = B * B - Float(4.0) * A * C;
            if (discriminant < Float(0.0))
                return std::nullopt;

            Float tolerance = Float(1e-15); // Use high precision tolerance
            if (abs(A) < tolerance)
            {
                // Linear case
                if (abs(B) < tolerance)
                    return std::nullopt;
                Float t = -C / B;
                return std::make_pair(t, t);
            }

            // Two solutions
            Float sqrt_disc = sqrt(discriminant);
            Float t1 = (-B - sqrt_disc) / (Float(2.0) * A);
            Float t2 = (-B + sqrt_disc) / (Float(2.0) * A);

            if (t1 > t2)
            {
                Float temp = t1;
                t1 = t2;
                t2 = temp;
            }
            return std::make_pair(t1, t2);
        }

        std::unique_ptr<ConvexSet> clone() const override
        {
            return std::make_unique<Ellipse>(*this);
        }

        std::string to_string() const
        {
            return "Ellipse(a=" + std::to_string(a()) + ", b=" + std::to_string(b()) +
                   ", d=" + std::to_string(d()) + ", px=" + std::to_string(px()) +
                   ", py=" + std::to_string(py()) + ")";
        }
    };

    /**
     * UnitDisk: Represents the unit disk |z| <= 1
     * Simple convex set implementation for the unit disk constraint
     */
    class UnitDisk : public ConvexSet
    {
    private:
        Ellipse _ellipse;

    public:
        UnitDisk() : _ellipse(1.0, 0.0, 1.0, 0.0, 0.0) {}

        // Ellipse property for to_upright_set_pair compatibility
        const Ellipse &ellipse() const { return _ellipse; }

        bool inside(const std::array<Float, 2> &u) const override
        {
            // throw std::runtime_error("inside");

            // |u|^2 <= 1 with adaptive tolerance (important at very high precision)
            Float norm_squared = u[0] * u[0] + u[1] * u[1];
            Float tol = Float(1e-30);
            return norm_squared <= Float(1.0) + tol;
        }

        std::optional<std::pair<Float, Float>> intersect(const std::array<Float, 2> &u0, const std::array<Float, 2> &v) const override
        {
            // throw std::runtime_error("intersect");

            // Faithful replication of Python UnitDisk.intersect via solve_quadratic
            // a = |v|^2, b = 2 * Re(v̄ u0), c = |u0|^2 - 1
            Float a = v[0] * v[0] + v[1] * v[1];
            Float b = Float(2.0) * (u0[0] * v[0] + u0[1] * v[1]);
            Float c = u0[0] * u0[0] + u0[1] * u0[1] - Float(1.0);

            // Use numerically stable quadratic solver identical in branching to Python
            auto roots_opt = solve_quadratic_optional(a, b, c);
            if (!roots_opt.has_value())
            {
                return std::nullopt;
            }
            auto [t1, t2] = roots_opt.value();
            if (t1 > t2)
                std::swap(t1, t2); // Ensure ascending order
            return std::make_pair(t1, t2);
        }

        bool intersects_line_segment(const std::array<Float, 2> &p, const std::array<Float, 2> &q) const override
        {

            // Simple implementation: check endpoints
            if (inside(p) || inside(q))
                return true;

            // For now, just a simple approximation
            // A more complete implementation would check actual line-circle intersection
            std::array<Float, 2> midpoint = {(p[0] + q[0]) / Float(2.0), (p[1] + q[1]) / Float(2.0)};
            return inside(midpoint);
        }

        std::unique_ptr<ConvexSet> clone() const override
        {
            return std::make_unique<UnitDisk>(*this);
        }
    };

    // EllipsePair class (moved from grid_op.hpp to avoid circular dependency)
    class EllipsePair
    {
    private:
        Ellipse _e1;
        Ellipse _e2;

    public:
        EllipsePair(const Ellipse &e1, const Ellipse &e2) : _e1(e1), _e2(e2) {}

        const Ellipse &e1() const { return _e1; }
        const Ellipse &e2() const { return _e2; }

        void set_e1(const Ellipse &e1) { _e1 = e1; }
        void set_e2(const Ellipse &e2) { _e2 = e2; }

        Float area() const
        {
            return _e1.area() + _e2.area();
        }

        Float skew() const
        {
            return _e1.skew() + _e2.skew();
        }

        Float bias() const
        {
            return _e2.bias() / _e1.bias();
        }

        // In-place variant that updates internal ellipses without extra allocations
        void apply_grid_op(const GridOp &g)
        {
            // Forward matrix (real parts)
            auto g_mat = g.to_mat();
            Float M00 = g_mat[0][0].real();
            Float M01 = g_mat[0][1].real();
            Float M10 = g_mat[1][0].real();
            Float M11 = g_mat[1][1].real();

            // Precompute inverse analytically for 2x2
            Float det = M00 * M11 - M01 * M10;
            Float tol = Float(1e-30);

            Float inv00, inv01, inv10, inv11;
            bool use_fallback = false;
            if (abs(det) < tol)
            {
                use_fallback = true;
            }
            else
            {
                inv00 = M11 / det;
                inv01 = -M01 / det;
                inv10 = -M10 / det;
                inv11 = M00 / det;
            }

            // Update first ellipse with g using matrix-variant
            _e1.transform_by_gridop_mat(M00, M01, M10, M11, false, use_fallback, inv00, inv01, inv10, inv11, tol, &g);

            // For the second ellipse we need the conjugated forward matrix (conj under sqrt(2)).
            // Compute it exactly from the underlying ZOmega components to avoid approximation errors.
            if (use_fallback)
            {
                // Fallback path: construct conjugated GridOp and let transform_by_gridop handle inverse computation
                GridOp g_conj = g.conj_sq2();
                _e2.transform_by_gridop(g_conj, true, use_fallback, inv00, inv01, inv10, inv11, tol);
            }
            else
            {
                // Use the underlying ZOmega elements (we have const accessors) and compute their conj_sq2 complex values
                const ZOmega &u0 = g.u0();
                const ZOmega &u1 = g.u1();
                ZOmega u0c_zo = u0.conj_sq2();
                ZOmega u1c_zo = u1.conj_sq2();
                Float M00c, M10c, M01c, M11c;
                u0c_zo.to_real_imag(M00c, M10c);
                u1c_zo.to_real_imag(M01c, M11c);

                _e2.transform_by_gridop_mat(M00c, M01c, M10c, M11c, true, false, inv00, inv01, inv10, inv11, tol, &g);
            }
        }

        std::string to_string() const
        {
            return "EllipsePair(e1=" + _e1.to_string() + ", e2=" + _e2.to_string() + ")";
        }
    };

} // namespace gridsynth
