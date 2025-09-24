#pragma once

#include <array>
#include <complex>
#include <optional>
#include <iostream>
#include <stdexcept>

#include "mymath.hpp"
#include "ring.hpp"
#include "types.hpp"

namespace gridsynth
{

    // Forward declarations
    class EllipsePair;

    /**
     * GridOp: Represents a grid operation with two ZOmega elements
     * Constraints: d0 + b0 + d1 + b1 : even
     *             a0 + c0 + a1 + c1 : even
     */
    class GridOp
    {
    private:
        ZOmega _u0;
        ZOmega _u1;

    public:
        GridOp(const ZOmega &u0, const ZOmega &u1) : _u0(u0), _u1(u1) {}

        // Properties
        // Return const references to avoid unnecessary copies in hot paths
        const ZOmega &u0() const { return _u0; }
        const ZOmega &u1() const { return _u1; }

        // Convenience accessors for components
        inline const Integer &a0() const noexcept { return _u0.a(); }
        inline const Integer &b0() const noexcept { return _u0.b(); }
        inline const Integer &c0() const noexcept { return _u0.c(); }
        inline const Integer &d0() const noexcept { return _u0.d(); }

        inline const Integer &a1() const noexcept { return _u1.a(); }
        inline const Integer &b1() const noexcept { return _u1.b(); }
        inline const Integer &c1() const noexcept { return _u1.c(); }
        inline const Integer &d1() const noexcept { return _u1.d(); }

        // String representation
        friend std::ostream &operator<<(std::ostream &os, const GridOp &g);

        // Properties
        ZOmega det_vec() const
        {
            return _u0.conj() * _u1;
        }

        bool is_special() const
        {
            ZOmega v = det_vec();
            return (v.a() + v.c() == 0) && (v.b() == 1 || v.b() == -1);
        }

        std::array<std::array<std::complex<Float>, 2>, 2> to_mat() const
        {
            // Compute real and imag directly to avoid constructing std::complex twice
            Float u0_r, u0_i, u1_r, u1_i;
            _u0.to_real_imag(u0_r, u0_i);
            _u1.to_real_imag(u1_r, u1_i);
            return {{{{u0_r, u1_r}},
                     {{u0_i, u1_i}}}};
        }

        // Multiplication operations
        GridOp operator*(const GridOp &other) const
        {
            // Use const refs to avoid copies from accessors
            const ZOmega &o0 = other.u0();
            const ZOmega &o1 = other.u1();
            return GridOp((*this) * o0, (*this) * o1);
        }

        ZOmega operator*(const ZOmega &other) const
        {
            // Cache local components to avoid repeated accessor calls
            const Integer a0_ = a0();
            const Integer b0_ = b0();
            const Integer c0_ = c0();
            const Integer d0_ = d0();

            const Integer a1_ = a1();
            const Integer b1_ = b1();
            const Integer c1_ = c1();
            const Integer d1_ = d1();

            const Integer oa = other.a();
            const Integer ob = other.b();
            const Integer oc = other.c();
            const Integer od = other.d();

            const Integer t1 = floordiv(c1_ - a1_ + c0_ - a0_, 2);
            const Integer t2 = floordiv(c1_ - a1_ - c0_ + a0_, 2);
            const Integer t3 = floordiv(b1_ + d1_ + b0_ + d0_, 2);
            const Integer t4 = floordiv(b1_ + d1_ - b0_ - d0_, 2);
            const Integer t5 = floordiv(c1_ + a1_ + c0_ + a0_, 2);
            const Integer t6 = floordiv(c1_ + a1_ - c0_ - a0_, 2);
            const Integer t7 = floordiv(b1_ - d1_ + b0_ - d0_, 2);
            const Integer t8 = floordiv(b1_ - d1_ - b0_ + d0_, 2);

            Integer new_d = (d0_ * od + d1_ * ob + t1 * oc + t2 * oa);
            Integer new_c = (c0_ * od + c1_ * ob + t3 * oc + t4 * oa);
            Integer new_b = (b0_ * od + b1_ * ob + t5 * oc + t6 * oa);
            Integer new_a = (a0_ * od + a1_ * ob + t7 * oc + t8 * oa);

            return ZOmega(new_a, new_b, new_c, new_d);
        }

        DOmega operator*(const DOmega &other) const
        {
            // other.u() is a ZOmega (by value); avoid unnecessary temporaries by using it directly
            return DOmega((*this) * other.u(), other.k());
        }

        // Inverse (for special matrices only)
        GridOp inv() const
        {
            if (!is_special())
            {
                throw std::domain_error("GridOp::inv: not a special matrix");
            }

            // Cache components
            const Integer a0_ = a0();
            const Integer b0_ = b0();
            const Integer c0_ = c0();
            const Integer d0_ = d0();

            const Integer a1_ = a1();
            const Integer b1_ = b1();
            const Integer c1_ = c1();
            const Integer d1_ = d1();

            Integer new_c0 = floordiv(c1_ + a1_ - c0_ - a0_, 2);
            Integer new_a0 = floordiv(-c1_ - a1_ - c0_ - a0_, 2);
            ZOmega new_u0(new_a0, -b0_, new_c0, b1_);

            Integer new_c1 = floordiv(-c1_ + a1_ + c0_ - a0_, 2);
            Integer new_a1 = floordiv(c1_ - a1_ + c0_ - a0_, 2);
            ZOmega new_u1(new_a1, d0_, new_c1, -d1_);

            if (det_vec().b() == -1)
            {
                new_u0 = ZOmega(0, 0, 0, 0) - new_u0;
                new_u1 = ZOmega(0, 0, 0, 0) - new_u1;
            }

            return GridOp(new_u0, new_u1);
        }

        // Power operation
        GridOp pow(Integer exp) const
        {
            if (exp < 0)
            {
                return inv().pow(-exp);
            }

            GridOp result(ZOmega(0, 0, 0, 1), ZOmega(0, 1, 0, 0)); // Identity
            GridOp base = *this;

            while (exp > 0)
            {
                if (exp & 1)
                {
                    result = result * base;
                }
                base = base * base;
                exp >>= 1;
            }
            return result;
        }

        // Adjoint
        GridOp adj() const
        {
            Integer new_c0 = floordiv(c1() - a1() + c0() - a0(), 2);
            Integer new_a0 = floordiv(c1() - a1() - c0() + a0(), 2);
            ZOmega new_u0(new_a0, d1(), new_c0, d0());

            Integer new_c1 = floordiv(c1() + a1() + c0() + a0(), 2);
            Integer new_a1 = floordiv(c1() + a1() - c0() - a0(), 2);
            ZOmega new_u1(new_a1, b1(), new_c1, b0());

            return GridOp(new_u0, new_u1);
        }

        // Conjugate under √2
        GridOp conj_sq2() const
        {
            return GridOp(_u0.conj_sq2(), _u1.conj_sq2());
        }
    };

    // String representation for GridOp
    inline std::ostream &operator<<(std::ostream &os, const GridOp &g)
    {
        os << "[[" << g.d0() << (g.c0() - g.a0() >= 0 ? "+" : "") << (g.c0() - g.a0()) << "/√2, "
           << g.d1() << (g.c1() - g.a1() >= 0 ? "+" : "") << (g.c1() - g.a1()) << "/√2],\n"
           << " [" << g.b0() << (g.c0() + g.a0() >= 0 ? "+" : "") << (g.c0() + g.a0()) << "/√2, "
           << g.b1() << (g.c1() + g.a1() >= 0 ? "+" : "") << (g.c1() + g.a1()) << "/√2]]";
        return os;
    }

} // namespace gridsynth
