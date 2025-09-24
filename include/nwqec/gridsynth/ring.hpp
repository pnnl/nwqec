#pragma once

#include <iostream>
#include <array>
#include <stdexcept>
#include <complex>
#include <optional>
#include <sstream>
#include <tuple>
#include <cassert>

#include "mymath.hpp"
#include "types.hpp"

namespace gridsynth
{

    // Forward declarations
    class ZRootTwo;
    class DRootTwo;
    class ZOmega;
    class DOmega;

    /**
     * ZRootTwo: Numbers of the form a + b√2 where a, b are integers
     */
    class ZRootTwo
    {
    private:
        Integer _a, _b;

    public:
        ZRootTwo(Integer a = 0, Integer b = 0) noexcept : _a(a), _b(b) {}

        // Properties
        const Integer &a() const noexcept { return _a; }
        const Integer &b() const noexcept { return _b; }

        // Factory methods
        static ZRootTwo from_int(Integer x)
        {
            return ZRootTwo(x, 0);
        }

        static ZRootTwo from_zomega(const ZOmega &x);

        // String representation
        friend std::ostream &operator<<(std::ostream &os, const ZRootTwo &z)
        {
            os << z._a;
            if (z._b >= 0)
            {
                os << "+" << z._b << "√2";
            }
            else
            {
                os << z._b << "√2";
            }
            return os;
        }

        // Equality operators
        bool operator==(const ZRootTwo &other) const
        {
            return _a == other._a && _b == other._b;
        }

        bool operator==(Integer other) const
        {
            return *this == ZRootTwo::from_int(other);
        }

        bool operator!=(const ZRootTwo &other) const
        {
            return !(*this == other);
        }

        // Comparison operators - using safer arithmetic to avoid overflow
        bool operator<(const ZRootTwo &other) const
        {
            Integer da = _a - other._a;
            Integer db = _b - other._b;

            if (db == 0)
            {
                return da < 0;
            }

            if (db > 0)
            {
                // need da < -db*sqrt2 -> requires da < 0 and da*da > 2*db*db
                if (da >= 0)
                    return false;
                Integer lhs = da * da;
                Integer rhs = Integer(2) * db * db;
                return lhs > rhs;
            }
            else // db < 0
            {
                Integer dbp = -db;
                if (da < 0)
                    return true;
                Integer lhs = da * da;
                Integer rhs = Integer(2) * dbp * dbp;
                return lhs < rhs;
            }
        }

        bool operator<=(const ZRootTwo &other) const
        {
            return *this < other || *this == other;
        }

        bool operator>(const ZRootTwo &other) const
        {
            return other < *this;
        }

        bool operator>=(const ZRootTwo &other) const
        {
            return !(*this < other);
        }

        // Arithmetic operations
        inline ZRootTwo operator+(const ZRootTwo &other) const noexcept
        {
            return ZRootTwo(_a + other._a, _b + other._b);
        }

        inline ZRootTwo operator-(const ZRootTwo &other) const noexcept
        {
            return ZRootTwo(_a - other._a, _b - other._b);
        }

        inline ZRootTwo operator-(Integer other) const noexcept
        {
            return *this - ZRootTwo::from_int(other);
        }

        inline ZRootTwo operator-() const noexcept
        {
            return ZRootTwo(-_a, -_b);
        }

        inline ZRootTwo operator*(const ZRootTwo &other) const noexcept
        {
            Integer new_a = _a * other._a + 2 * _b * other._b;
            Integer new_b = _a * other._b + _b * other._a;
            return ZRootTwo(new_a, new_b);
        }

        // Scalar multiplication: scale both integer coefficients directly.
        inline ZRootTwo operator*(Integer other) const noexcept
        {
            return ZRootTwo(_a * other, _b * other);
        }

        // Properties computed on demand
        Integer parity() const noexcept
        {
            return _a & 1;
        }

        Integer norm() const noexcept
        {
            return _a * _a - 2 * _b * _b;
        }

        Float to_real() const noexcept
        {
            return Float(_a) + SQRT2 * Float(_b);
        }

        ZRootTwo conj_sq2() const
        {
            return ZRootTwo(_a, -_b);
        }

        // Inverse (only for units)
        ZRootTwo inv() const
        {
            Integer n = norm();
            if (n == 1)
            {
                return conj_sq2();
            }
            else if (n == -1)
            {
                return -conj_sq2();
            }
            else
            {
                throw std::domain_error("ZRootTwo::inv: not a unit");
            }
        }

        // Power operation
        ZRootTwo pow(Integer exp) const
        {
            if (exp < 0)
            {
                return inv().pow(-exp);
            }
            ZRootTwo result = ZRootTwo::from_int(1);
            ZRootTwo base = *this;
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

        // Square root (returns optional value to avoid memory management issues)
        std::optional<ZRootTwo> sqrt() const
        {
            Integer n = norm();
            if (n < 0 || _a < 0)
            {
                return std::nullopt;
            }

            Integer r = floorsqrt(n);
            Integer a1 = floorsqrt(floordiv(_a + r, 2));
            Integer b1 = floorsqrt(floordiv(_a - r, 4));
            Integer a2 = floorsqrt(floordiv(_a - r, 2));
            Integer b2 = floorsqrt(floordiv(_a + r, 4));

            ZRootTwo w1, w2;
            if (sign(_a) * sign(_b) >= 0)
            {
                w1 = ZRootTwo(a1, b1);
                w2 = ZRootTwo(a2, b2);
            }
            else
            {
                w1 = ZRootTwo(a1, -b1);
                w2 = ZRootTwo(a2, -b2);
            }

            if (*this == w1 * w1)
            {
                return w1;
            }
            else if (*this == w2 * w2)
            {
                return w2;
            }
            else
            {
                return std::nullopt;
            }
        }

        // Division with remainder
        std::pair<ZRootTwo, ZRootTwo> divmod(const ZRootTwo &other) const
        {
            ZRootTwo p = (*this) * other.conj_sq2();
            Integer k = other.norm();
            ZRootTwo q = ZRootTwo(rounddiv(p.a(), k), rounddiv(p.b(), k));
            ZRootTwo r = *this - other * q;

            return {q, r};
        }

        std::pair<ZRootTwo, ZRootTwo> divmod(Integer other) const
        {
            return divmod(ZRootTwo::from_int(other));
        }

        ZRootTwo operator/(const ZRootTwo &other) const
        {
            return divmod(other).first;
        }

        ZRootTwo operator%(const ZRootTwo &other) const
        {
            return divmod(other).second;
        }

        // Static methods for GCD operations
        static bool sim(const ZRootTwo &a, const ZRootTwo &b)
        {
            return (a % b == ZRootTwo(0, 0)) && (b % a == ZRootTwo(0, 0));
        }

        static ZRootTwo gcd(ZRootTwo a, ZRootTwo b)
        {
            const ZRootTwo zero = ZRootTwo::from_int(0);
            while (!(b == zero))
            {
                auto r = a.divmod(b).second;
                a = b;
                b = r;
            }
            return a;
        }
    };

    /**
     * DRootTwo: Numbers of the form α / √2^k where α is a ZRootTwo
     */
    class DRootTwo
    {
    private:
        ZRootTwo _alpha;
        Integer _k;

    public:
        DRootTwo(const ZRootTwo &alpha = ZRootTwo(), Integer k = 0) noexcept : _alpha(alpha), _k(k) {}

        // Constructor from integer (for convenience)
        DRootTwo(int val) : _alpha(ZRootTwo::from_int(Integer(val))), _k(0) {}
        DRootTwo(Integer val) : _alpha(ZRootTwo::from_int(val)), _k(0) {}

        // Properties
        const ZRootTwo &alpha() const noexcept { return _alpha; }
        const Integer &k() const noexcept { return _k; }

        // Factory methods
        static DRootTwo from_int(Integer x)
        {
            return DRootTwo(ZRootTwo::from_int(x), 0);
        }

        static DRootTwo from_zroottwo(const ZRootTwo &x)
        {
            return DRootTwo(x, 0);
        }

        static DRootTwo from_zomega(const ZOmega &x);

        static DRootTwo fromDOmega(const DOmega &x);

        // String representation
        friend std::ostream &operator<<(std::ostream &os, const DRootTwo &d)
        {
            os << d._alpha << " / √2^" << d._k;
            return os;
        }

        // Equality operators
        bool operator==(const DRootTwo &other) const
        {
            // Normalize to the same denominator exponent before comparing.
            if (_k < other._k)
            {
                return renew_denomexp(other._k) == other;
            }
            else if (_k > other._k)
            {
                return *this == other.renew_denomexp(_k);
            }
            else
            {
                return _alpha == other._alpha; // _k equal already
            }
        }

        // Arithmetic operations
        inline DRootTwo operator+(const DRootTwo &other) const noexcept
        {
            if (_k < other._k)
            {
                return renew_denomexp(other._k) + other;
            }
            else if (_k > other._k)
            {
                return *this + other.renew_denomexp(_k);
            }
            else
            {
                return DRootTwo(_alpha + other._alpha, _k);
            }
        }

        inline DRootTwo operator-(const DRootTwo &other) const noexcept
        {
            return *this + (-other);
        }

        inline DRootTwo operator-() const noexcept
        {
            return DRootTwo(-_alpha, _k);
        }

        inline DRootTwo operator*(const DRootTwo &other) const noexcept
        {
            return DRootTwo(_alpha * other._alpha, _k + other._k);
        }

        // Scalar multiplication: scale the underlying ZRootTwo alpha.
        inline DRootTwo operator*(Integer other) const noexcept
        {
            return DRootTwo(_alpha * other, _k);
        }

        // Utility methods
        DRootTwo renew_denomexp(Integer new_k) const
        {
            ZRootTwo new_alpha = mul_by_sqrt2_power(new_k - _k).alpha();
            return DRootTwo(new_alpha, new_k);
        }

        DRootTwo mul_by_inv_sqrt2() const
        {
            if (!(_alpha.a() & 1))
            {
                ZRootTwo new_alpha(_alpha.b(), _alpha.a() >> 1);
                return DRootTwo(new_alpha, _k);
            }
            else
            {
                throw std::domain_error("mul_by_inv_sqrt2: invalid operation");
            }
        }

        DRootTwo mul_by_sqrt2_power(Integer d) const
        {
            if (d < 0)
            {
                if (d == -1)
                {
                    return mul_by_inv_sqrt2();
                }
                Integer abs_d = -d;
                Integer d_div_2 = abs_d >> 1;
                Integer d_mod_2 = abs_d & 1;
                if (d_mod_2 == 0)
                {
                    Integer bit = (1LL << d_div_2) - 1;
                    if ((_alpha.a() & bit) == 0 && (_alpha.b() & bit) == 0)
                    {
                        ZRootTwo new_alpha(_alpha.a() >> d_div_2, _alpha.b() >> d_div_2);
                        return DRootTwo(new_alpha, _k);
                    }
                    else
                    {
                        throw std::domain_error("mul_by_sqrt2_power: invalid operation");
                    }
                }
                else
                {
                    Integer bit = (1LL << d_div_2) - 1;
                    Integer bit2 = (1LL << (d_div_2 + 1)) - 1;
                    if ((_alpha.a() & bit2) == 0 && (_alpha.b() & bit) == 0)
                    {
                        ZRootTwo new_alpha(_alpha.b() >> d_div_2, _alpha.a() >> (d_div_2 + 1));
                        return DRootTwo(new_alpha, _k);
                    }
                    else
                    {
                        throw std::domain_error("mul_by_sqrt2_power: invalid operation");
                    }
                }
            }
            else
            {
                Integer d_div_2 = d >> 1;
                Integer d_mod_2 = d & 1;
                // Scale by 2^(d_div_2) via exact left shifts of integer coefficients
                ZRootTwo new_alpha(_alpha.a() << static_cast<int>(d_div_2),
                                   _alpha.b() << static_cast<int>(d_div_2));
                if (d_mod_2)
                {
                    new_alpha = new_alpha * ZRootTwo(0, 1);
                }
                return DRootTwo(new_alpha, _k);
            }
        }

        DRootTwo mul_by_sqrt2_power_renewing_denomexp(Integer d) const
        {
            if (d > _k)
            {
                throw std::domain_error("mul_by_sqrt2_power_renewing_denomexp: invalid operation");
            }
            return DRootTwo(_alpha, _k - d);
        }

        // Properties
        Integer parity() const
        {
            return _alpha.parity();
        }

        Float scale() const
        {
            return pow_sqrt2(_k);
        }

        Float to_real() const
        {
            return _alpha.to_real() / scale();
        }

        DRootTwo conj_sq2() const
        {
            if (_k & 1)
            {
                return DRootTwo(-_alpha.conj_sq2(), _k);
            }
            else
            {
                return DRootTwo(_alpha.conj_sq2(), _k);
            }
        }

        static DRootTwo power_of_inv_sqrt2(Integer k)
        {
            return DRootTwo(ZRootTwo::from_int(1), k);
        }
    };

    /**
     * ZOmega: Numbers in Z[ω] where ω = e^(iπ/4)
     */
    class ZOmega
    {
    private:
        Integer _a, _b, _c, _d; // aω³ + bω² + cω + d

    public:
        ZOmega(Integer a = 0, Integer b = 0, Integer c = 0, Integer d = 0)
            : _a(a), _b(b), _c(c), _d(d) {}

        // Properties
        const Integer &a() const noexcept { return _a; }
        const Integer &b() const noexcept { return _b; }
        const Integer &c() const noexcept { return _c; }
        const Integer &d() const noexcept { return _d; }

        inline std::array<Integer, 4> coef() const noexcept
        {
            return {_d, _c, _b, _a};
        }

        // Factory methods
        static ZOmega from_int(Integer x)
        {
            return ZOmega(0, 0, 0, x);
        }

        static ZOmega from_zroottwo(const ZRootTwo &x)
        {
            return ZOmega(-x.b(), 0, x.b(), x.a());
        }

        // String representation
        friend std::ostream &operator<<(std::ostream &os, const ZOmega &z)
        {
            os << z._a << "ω³";
            if (z._b >= 0)
                os << "+";
            os << z._b << "ω²";
            if (z._c >= 0)
                os << "+";
            os << z._c << "ω";
            if (z._d >= 0)
                os << "+";
            os << z._d;
            return os;
        }

        // Equality operators
        bool operator==(const ZOmega &other) const
        {
            return _a == other._a && _b == other._b && _c == other._c && _d == other._d;
        }

        bool operator==(const ZRootTwo &other) const
        {
            return *this == ZOmega::from_zroottwo(other);
        }

        // Arithmetic operations
        inline ZOmega operator+(const ZOmega &other) const noexcept
        {
            return ZOmega(_a + other._a, _b + other._b, _c + other._c, _d + other._d);
        }

        inline ZOmega operator-(const ZOmega &other) const noexcept
        {
            return ZOmega(_a - other._a, _b - other._b, _c - other._c, _d - other._d);
        }

        inline ZOmega operator-() const noexcept
        {
            return ZOmega(-_a, -_b, -_c, -_d);
        }

        inline ZOmega operator*(const ZOmega &other) const noexcept
        {
            // Multiply modulo x^4 + 1 with coefficient order [d, c, b, a]
            const Integer &a0 = _a, &b0 = _b, &c0 = _c, &d0 = _d;
            const Integer &a1 = other._a, &b1 = other._b, &c1 = other._c, &d1 = other._d;

            const Integer r0 = d0 * d1;
            const Integer r1 = d0 * c1 + c0 * d1;
            const Integer r2 = d0 * b1 + c0 * c1 + b0 * d1;
            const Integer r3 = d0 * a1 + c0 * b1 + b0 * c1 + a0 * d1;
            const Integer r4 = c0 * a1 + b0 * b1 + a0 * c1;
            const Integer r5 = b0 * a1 + a0 * b1;
            const Integer r6 = a0 * a1;

            const Integer nd = r0 - r4;
            const Integer nc = r1 - r5;
            const Integer nb = r2 - r6;
            const Integer na = r3;

            return ZOmega(na, nb, nc, nd);
        }

        inline ZOmega operator*(Integer other) const noexcept
        {
            // Multiplying by an integer scalar x is equivalent to scaling each
            // coefficient by x.
            return ZOmega(_a * other, _b * other, _c * other, _d * other);
        }

        // Multiplicative inverse (for units only)
        ZOmega inv() const
        {
            if (norm() == 1)
            {
                return conj_sq2() * conj() * conj().conj_sq2();
            }
            else
            {
                throw std::domain_error("ZOmega::inv: not a unit");
            }
        }

        // Power operation
        ZOmega pow(Integer exp) const
        {
            if (exp < 0)
            {
                throw std::domain_error("ZOmega::pow: negative exponent not supported");
            }

            ZOmega result = ZOmega::from_int(1);
            ZOmega base = *this;
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

        // Division with remainder
        std::pair<ZOmega, ZOmega> divmod(const ZOmega &other) const
        {
            ZOmega p = (*this) * other.conj() * other.conj().conj_sq2() * other.conj_sq2();
            Integer k = other.norm();
            ZOmega q(rounddiv(p.a(), k), rounddiv(p.b(), k),
                     rounddiv(p.c(), k), rounddiv(p.d(), k));
            ZOmega r = *this - other * q;

            return {q, r};
        }

        // Omega operations
        ZOmega mul_by_omega() const
        {
            return ZOmega(_b, _c, _d, -_a);
        }

        ZOmega mul_by_omega_inv() const
        {
            return ZOmega(-_d, _a, _b, _c);
        }

        ZOmega mul_by_omega_power(Integer n) const noexcept
        {
            // Reduce modulo 8 and compute explicit rotations
            n &= 0b111;
            switch (n)
            {
            case 0:
                return *this;
            case 1:
                return mul_by_omega(); // (b,c,d,-a)
            case 2:
                // apply mul_by_omega twice: (c,d,-a,-b)
                return ZOmega(_c, _d, -_a, -_b);
            case 3:
                // (d,-a,-b,-c)
                return ZOmega(_d, -_a, -_b, -_c);
            case 4:
                return ZOmega(-_a, -_b, -_c, -_d);
            case 5:
                // (-b,-c,-d,a)
                return ZOmega(-_b, -_c, -_d, _a);
            case 6:
                // (-c,-d,a,b)
                return ZOmega(-_c, -_d, _a, _b);
            case 7:
                // (-d,a,b,c)
                return ZOmega(-_d, _a, _b, _c);
            default:
                return *this;
            }
        }

        // Properties
        Integer residue() const
        {
            return ((_a & 1) << 3) | ((_b & 1) << 2) | ((_c & 1) << 1) | (_d & 1);
        }

        Integer norm() const
        {
            Integer sum_squares = _a * _a + _b * _b + _c * _c + _d * _d;
            Integer cross_term = _a * _b + _b * _c + _c * _d - _d * _a;
            return sum_squares * sum_squares - 2 * cross_term * cross_term;
        }

        std::complex<Float> to_complex() const noexcept
        {
            Float real_part = Float(_d) + SQRT2 * Float(_c - _a) / Float(2.0);
            Float imag_part = Float(_b) + SQRT2 * Float(_c + _a) / Float(2.0);
            return std::complex<Float>(real_part, imag_part);
        }

        inline void to_real_imag(Float &out_real, Float &out_imag) const noexcept
        {
            out_real = Float(_d) + SQRT2 * Float(_c - _a) / Float(2.0);
            out_imag = Float(_b) + SQRT2 * Float(_c + _a) / Float(2.0);
        }

        ZOmega conj() const noexcept
        {
            return ZOmega(-_c, -_b, -_a, _d);
        }

        ZOmega conj_sq2() const noexcept
        {
            return ZOmega(-_a, _b, -_c, _d);
        }

        // Static GCD methods
        static std::tuple<ZOmega, ZOmega, ZOmega> ext_gcd(ZOmega a, ZOmega b)
        {
            ZOmega x = ZOmega::from_int(1);
            ZOmega y = ZOmega::from_int(0);
            ZOmega z = ZOmega::from_int(0);
            ZOmega w = ZOmega::from_int(1);

            while (!(b == ZOmega(0, 0, 0, 0)))
            {
                auto [q, r] = a.divmod(b);

                auto new_x = y;
                auto new_y = x - y * q;
                auto new_z = w;
                auto new_w = z - w * q;

                x = new_x;
                y = new_y;
                z = new_z;
                w = new_w;
                a = b;
                b = r;
            }
            return {x, z, a};
        }

        static ZOmega gcd(ZOmega a, ZOmega b)
        {
            const ZOmega zero(0,0,0,0);
            while (!(b == zero)) {
                auto r = a.divmod(b).second;
                a = b;
                b = r;
            }
            return a;
        }

        static bool sim(const ZOmega &a, const ZOmega &b)
        {
            return (a.divmod(b).second == ZOmega(0, 0, 0, 0)) &&
                   (b.divmod(a).second == ZOmega(0, 0, 0, 0));
        }
    };

    /**
     * DOmega: Numbers of the form u / √2^k where u is a ZOmega
     */
    class DOmega
    {
    private:
        ZOmega _u;
        Integer _k;

    public:
        DOmega(const ZOmega &u = ZOmega(), Integer k = 0) noexcept : _u(u), _k(k) {}

        // Properties
        const ZOmega &u() const noexcept { return _u; }
        const Integer &k() const noexcept { return _k; }

        // Factory methods
        static DOmega from_int(Integer x)
        {
            return DOmega(ZOmega::from_int(x), 0);
        }

        static DOmega from_zroottwo(const ZRootTwo &x)
        {
            return DOmega(ZOmega::from_zroottwo(x), 0);
        }

        static DOmega from_droottwo(const DRootTwo &x)
        {
            return DOmega(ZOmega::from_zroottwo(x.alpha()), x.k());
        }

        static DOmega from_droottwo_vector(const DRootTwo &x, const DRootTwo &y, Integer k)
        {
            DOmega dx = DOmega::from_droottwo(x);
            DOmega dy = DOmega::from_droottwo(y) * ZOmega(0, 1, 0, 0);
            return (dx + dy).renew_denomexp(k);
        }

        static DOmega from_zomega(const ZOmega &x)
        {
            return DOmega(x, 0);
        }

        // String representation
        friend std::ostream &operator<<(std::ostream &os, const DOmega &d)
        {
            os << d._u << " / √2^" << d._k;
            return os;
        }

        // Python-style string representation
        std::string to_string() const
        {
            std::ostringstream oss;
            oss << "DOmega(ZOmega(" << _u.a() << ", " << _u.b() << ", " << _u.c() << ", " << _u.d() << "), " << _k << ")";
            return oss.str();
        }

        // Equality operators
        bool operator==(const DOmega &other) const
        {
            if (_k < other._k)
            {
                return renew_denomexp(other._k) == other;
            }
            else if (_k > other._k)
            {
                return *this == other.renew_denomexp(_k);
            }
            else
            {
                return _u == other._u && _k == other._k;
            }
        }

        // Arithmetic operations
        inline DOmega operator+(const DOmega &other) const noexcept
        {
            if (_k < other._k)
            {
                return renew_denomexp(other._k) + other;
            }
            else if (_k > other._k)
            {
                return *this + other.renew_denomexp(_k);
            }
            else
            {
                return DOmega(_u + other._u, _k);
            }
        }

        inline DOmega operator-(const DOmega &other) const noexcept
        {
            return *this + (-other);
        }

        inline DOmega operator-() const noexcept
        {
            return DOmega(-_u, _k);
        }

        inline DOmega operator*(const DOmega &other) const noexcept
        {
            return DOmega(_u * other._u, _k + other._k);
        }

        // Scalar multiplication: scale the ZOmega numerator directly.
        inline DOmega operator*(Integer other) const noexcept
        {
            return DOmega(_u * other, _k);
        }

        // Utility methods
        DOmega renew_denomexp(Integer new_k) const
        {
            ZOmega new_u = mul_by_sqrt2_power(new_k - _k).u();
            return DOmega(new_u, new_k);
        }

        DOmega reduce_denomexp() const
        {
            Integer k_a = (_u.a() == 0) ? _k : ntz(_u.a());
            Integer k_b = (_u.b() == 0) ? _k : ntz(_u.b());
            Integer k_c = (_u.c() == 0) ? _k : ntz(_u.c());
            Integer k_d = (_u.d() == 0) ? _k : ntz(_u.d());
            Integer reduce_k = min(min(k_a, k_b), min(k_c, k_d));
            Integer new_k = _k - reduce_k * 2;
            Integer bit = (1LL << (reduce_k + 1)) - 1;
            if (((_u.c() + _u.a()) & bit) == 0 && ((_u.b() + _u.d()) & bit) == 0)
            {
                new_k -= 1;
            }
            return renew_denomexp(max(Integer(0), new_k));
        }

        DOmega mul_by_inv_sqrt2() const
        {
            if (!((_u.b() + _u.d()) & 1) && !((_u.c() + _u.a()) & 1))
            {
                ZOmega new_u((_u.b() - _u.d()) >> 1,
                             (_u.c() + _u.a()) >> 1,
                             (_u.b() + _u.d()) >> 1,
                             (_u.c() - _u.a()) >> 1);
                return DOmega(new_u, _k);
            }
            else
            {
                throw std::domain_error("mul_by_inv_sqrt2: invalid operation");
            }
        }

        DOmega mul_by_sqrt2_power(Integer d) const
        {
            if (d < 0)
            {
                if (d == -1)
                {
                    return mul_by_inv_sqrt2();
                }
                Integer abs_d = -d;
                Integer d_div_2 = abs_d >> 1;
                Integer d_mod_2 = abs_d & 1;
                if (d_mod_2 == 0)
                {
                    Integer bit = (1LL << d_div_2) - 1;
                    if ((_u.a() & bit) == 0 && (_u.b() & bit) == 0 &&
                        (_u.c() & bit) == 0 && (_u.d() & bit) == 0)
                    {
                        ZOmega new_u(_u.a() >> d_div_2, _u.b() >> d_div_2,
                                     _u.c() >> d_div_2, _u.d() >> d_div_2);
                        return DOmega(new_u, _k);
                    }
                    else
                    {
                        throw std::domain_error("mul_by_sqrt2_power: invalid operation");
                    }
                }
                else
                {
                    Integer bit = (1LL << (d_div_2 + 1)) - 1;
                    if (((_u.b() - _u.d()) & bit) == 0 &&
                        ((_u.c() + _u.a()) & bit) == 0 &&
                        ((_u.b() + _u.d()) & bit) == 0 &&
                        ((_u.c() - _u.a()) & bit) == 0)
                    {
                        ZOmega new_u((_u.b() - _u.d()) >> (d_div_2 + 1),
                                     (_u.c() + _u.a()) >> (d_div_2 + 1),
                                     (_u.b() + _u.d()) >> (d_div_2 + 1),
                                     (_u.c() - _u.a()) >> (d_div_2 + 1));
                        return DOmega(new_u, _k);
                    }
                    else
                    {
                        throw std::domain_error("mul_by_sqrt2_power: invalid operation");
                    }
                }
            }
            else
            {
                Integer d_div_2 = d >> 1;
                Integer d_mod_2 = d & 1;
                // Scale by 2^(d_div_2) via exact left shifts per component
                ZOmega new_u(_u.a() << static_cast<int>(d_div_2),
                             _u.b() << static_cast<int>(d_div_2),
                             _u.c() << static_cast<int>(d_div_2),
                             _u.d() << static_cast<int>(d_div_2));
                if (d_mod_2)
                {
                    new_u = new_u * ZOmega(-1, 0, 1, 0); // multiply by √2 in Z[ω]
                }
                return DOmega(new_u, _k);
            }
        }

        // Omega operations
        DOmega mul_by_omega() const noexcept
        {
            return DOmega(_u.mul_by_omega(), _k);
        }

        DOmega mul_by_omega_inv() const noexcept
        {
            return DOmega(_u.mul_by_omega_inv(), _k);
        }

        DOmega mul_by_omega_power(Integer n) const noexcept
        {
            return DOmega(_u.mul_by_omega_power(n), _k);
        }

        // Properties
        Float scale() const noexcept
        {
            return pow_sqrt2(_k);
        }

        Integer squared_scale() const noexcept
        {
            return 1LL << _k;
        }

        Integer residue() const noexcept
        {
            return _u.residue();
        }

        std::complex<Float> to_complex() const
        {
            Float ur, ui;
            _u.to_real_imag(ur, ui);
            const Float inv_scale = Float(1.0) / scale();
            return std::complex<Float>(ur * inv_scale, ui * inv_scale);
        }

        // Direct coordinate accessors (avoid constructing std::complex)
        inline Float real_part() const
        {
            // Real(u) = d + (c - a)/√2 then scaled by 1/scale()
            const Float numerator = Float(_u.d()) + (Float(_u.c() - _u.a()) * SQRT2 / Float(2.0));
            const Float inv_scale = Float(1.0) / scale();
            return numerator * inv_scale;
        }

        inline Float imag_part() const
        {
            // Imag(u) = b + (c + a)/√2 then scaled by 1/scale()
            const Float numerator = Float(_u.b()) + (Float(_u.c() + _u.a()) * SQRT2 / Float(2.0));
            const Float inv_scale = Float(1.0) / scale();
            return numerator * inv_scale;
        }

        inline std::array<Float, 2> coords() const
        {
            const Float inv_scale = Float(1.0) / scale();
            const Float real_numer = Float(_u.d()) + (Float(_u.c() - _u.a()) * SQRT2 / Float(2.0));
            const Float imag_numer = Float(_u.b()) + (Float(_u.c() + _u.a()) * SQRT2 / Float(2.0));
            return {real_numer * inv_scale, imag_numer * inv_scale};
        }

        // Write coords into provided outputs to avoid constructing std::array
        inline void coords_into(Float &out_real, Float &out_imag) const noexcept
        {
            const Float inv_scale = Float(1.0) / scale();
            const Float real_numer = Float(_u.d()) + (Float(_u.c() - _u.a()) * SQRT2 / Float(2.0));
            const Float imag_numer = Float(_u.b()) + (Float(_u.c() + _u.a()) * SQRT2 / Float(2.0));
            out_real = real_numer * inv_scale;
            out_imag = imag_numer * inv_scale;
        }

        // Variants that reuse precomputed invariants
        inline std::array<Float, 2> coords_with(const Float &inv_scale, const Float &sqrt2_over_2) const
        {
            const Float real_numer = Float(_u.d()) + (Float(_u.c() - _u.a()) * sqrt2_over_2);
            const Float imag_numer = Float(_u.b()) + (Float(_u.c() + _u.a()) * sqrt2_over_2);
            return {real_numer * inv_scale, imag_numer * inv_scale};
        }

        inline void coords_into_with(const Float &inv_scale, const Float &sqrt2_over_2, Float &out_real, Float &out_imag) const noexcept
        {
            const Float real_numer = Float(_u.d()) + (Float(_u.c() - _u.a()) * sqrt2_over_2);
            const Float imag_numer = Float(_u.b()) + (Float(_u.c() + _u.a()) * sqrt2_over_2);
            out_real = real_numer * inv_scale;
            out_imag = imag_numer * inv_scale;
        }

        DOmega conj() const
        {
            return DOmega(_u.conj(), _k);
        }

        DOmega conj_sq2() const
        {
            if (_k & 1)
            {
                return DOmega(-_u.conj_sq2(), _k);
            }
            else
            {
                return DOmega(_u.conj_sq2(), _k);
            }
        }
    };

    // Forward declaration implementations
    inline ZRootTwo ZRootTwo::from_zomega(const ZOmega &x)
    {
        assert(x.b() == 0 && x.a() == -x.c());
        return ZRootTwo(x.d(), x.c());
    }

    inline DRootTwo DRootTwo::from_zomega(const ZOmega &x)
    {
        return DRootTwo(ZRootTwo::from_zomega(x), 0);
    }

    inline DRootTwo DRootTwo::fromDOmega(const DOmega &x)
    {
        return DRootTwo(ZRootTwo::from_zomega(x.u()), x.k());
    }

    // Constants
    inline const ZRootTwo LAMBDA(1, 1);
    inline const ZOmega OMEGA(0, 0, 1, 0);
    inline const std::array<ZOmega, 8> OMEGA_POWER = {
        ZOmega(0, 0, 0, 1),  // ω^0 = 1
        ZOmega(0, 0, 1, 0),  // ω^1 = ω
        ZOmega(0, 1, 0, 0),  // ω^2
        ZOmega(1, 0, 0, 0),  // ω^3
        ZOmega(0, 0, 0, -1), // ω^4 = -1
        ZOmega(0, 0, -1, 0), // ω^5 = -ω
        ZOmega(0, -1, 0, 0), // ω^6 = -ω^2
        ZOmega(-1, 0, 0, 0)  // ω^7 = -ω^3
    };

} // namespace gridsynth
