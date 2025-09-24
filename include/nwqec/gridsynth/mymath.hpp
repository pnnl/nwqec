#pragma once

#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <optional>

#include "types.hpp"
// NOTE: Added floor, ceil, floorlog to support grid problem port.

namespace gridsynth
{

    // High precision SQRT2 constant (100+ digits)
    const Float SQRT2 = Float("1.4142135623730950488016887242096980785696718753769480731766797379907324784621070388503875343276415727");
    const Float PI = Float("3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679");
    // Number of trailing zeros in binary representation
    inline Integer ntz(Integer n) noexcept
    {
        if (n == Integer(0))
            return Integer(0);

        Integer count = Integer(0);
        while ((n & Integer(1)) == Integer(0))
        {
            n >>= 1;
            count = count + Integer(1);
        }
        return count;
    }

    inline Integer sign(Integer x) noexcept
    {
        return (x > Integer(0)) ? Integer(1) : ((x < Integer(0)) ? Integer(-1) : Integer(0));
    }

    inline Integer floorsqrt(const Integer &x)
    {
        if (x < Integer(0))
            throw std::invalid_argument("floorsqrt: negative input");
        if (x == Integer(0))
            return Integer(0);

        return x.floorsqrt();
    }

    // Helper functions for Float to Integer conversion (efficient versions)
    inline Integer round_to_integer(Float x)
    {
        return round_to_gmpinteger(x);
    }

    inline Integer ceil_to_integer(Float x)
    {
        return ceil_to_gmpinteger(x);
    }

    inline Integer floor_to_integer(Float x)
    {
        return floor_to_gmpinteger(x);
    }

    inline Integer floorsqrt(Float x)
    {
        return floor_to_integer(sqrt(x));
    }

    // Floor division helper - matches Python's // behavior
    inline Integer floordiv(Integer x, Integer y)
    {
        Integer result = x / y;
        // If there's a remainder and the signs differ, subtract 1
        if ((x % y != Integer(0)) && ((x < Integer(0)) != (y < Integer(0))))
        {
            result -= Integer(1);
        }
        return result;
    }

    inline Integer rounddiv(Integer x, Integer y)
    {
        if (y > Integer(0))
        {
            return floordiv(x + floordiv(y, Integer(2)), y);
        }
        else
        {
            return floordiv(x - floordiv(-y, Integer(2)), y);
        }
    }

    inline Integer gcd(Integer a, Integer b)
    {
        // Implement GCD algorithm since std::gcd may not work with our custom Integer
        if (a < Integer(0))
            a = -a;
        if (b < Integer(0))
            b = -b;

        while (b != Integer(0))
        {
            Integer temp = b;
            b = a % b;
            a = temp;
        }
        return a;
    }

    template <typename T>
    inline T max(T a, T b)
    {
        return (a > b) ? a : b;
    }

    template <typename T>
    inline T min(T a, T b)
    {
        return (a < b) ? a : b;
    }

    inline Float pow_sqrt2(Integer k)
    {
        // Original implementation: compute (sqrt(2))^k using arbitrary-precision
        // arithmetic without any static caching. Handles negative k by inversion.
        if (k < Integer(0))
        {
            return Float("1.0") / pow_sqrt2(-k);
        }

        Integer half = k >> 1;        // floor(k/2)
        Integer odd = k & Integer(1); // k % 2

        // Compute 2^half via exponentiation by squaring using arbitrary-precision Float.
        Float result("1.0");
        if (half != Integer(0))
        {
            Float base("2.0");
            Integer e = half;
            while (e > Integer(0))
            {
                if (e & Integer(1))
                {
                    result = result * base;
                }
                base = base * base;
                e >>= 1;
            }
        }
        if (odd != Integer(0))
        {
            result = result * SQRT2; // multiply by √2 for odd k
        }
        return result;
    }

    // Replicates Python floorlog(x, y) algorithm returning (n, r) where algorithmically
    inline std::pair<Integer, Float> floorlog(Float x, Float y)
    {
        if (x <= Float("0.0"))
        {
            throw std::runtime_error("floating point error");
        }
        Float tmp = y;
        Integer m = 0;
        while (x >= tmp || x * tmp < Float("1.0"))
        {
            // std::cout << x << " " << tmp << "\n";
            tmp *= tmp;
            ++m;
        }
        // Build powers list (reversed accumulation mimic)
        std::vector<Float> pow_y; // length m
        if (m > 0)
        {
            // emulate accumulate with squaring chain
            Float cur = y;
            pow_y.push_back(cur);
            for (Integer i = 1; i < m; ++i)
            {
                cur = cur * cur;
                pow_y.push_back(cur);
            }
            // reverse order
            std::reverse(pow_y.begin(), pow_y.end());
        }
        Integer n;
        Float r;
        if (x >= Float("1.0"))
        {
            n = Integer(0);
            r = x;
        }
        else
        {
            n = Integer(-1);
            r = x * tmp; // matches Python branch
        }
        for (Float p : pow_y)
        {
            n <<= 1;
            if (r > p)
            {
                r /= p;
                n += Integer(1);
            }
        }
        return {n, r};
    }

    inline std::pair<Float, Float> solve_quadratic(Float a, Float b, Float c)
    {
        if (a < 0)
        {
            a = -a;
            b = -b;
            c = -c;
        }
        Float discriminant = b * b - 4 * a * c;
        if (discriminant < 0)
        {
            throw std::invalid_argument("solve_quadratic: negative discriminant");
        }

        Float sqrt_disc = sqrt(discriminant);
        Float s1 = -b - sqrt_disc;
        Float s2 = -b + sqrt_disc;

        if (b >= 0)
        {
            return {s1 / (2 * a), s2 / (2 * a)};
        }
        else
        {
            return {(2 * c) / s2, (2 * c) / s1};
        }
    }

    // Version that returns optional instead of throwing exception
    inline std::optional<std::pair<Float, Float>> solve_quadratic_optional(Float a, Float b, Float c)
    {
        if (a < 0)
        {
            a = -a;
            b = -b;
            c = -c;
        }
        Float discriminant = b * b - 4 * a * c;
        if (discriminant < 0)
        {
            return std::nullopt;
        }

        Float sqrt_disc = sqrt(discriminant);
        Float s1 = -b - sqrt_disc;
        Float s2 = -b + sqrt_disc;

        if (b >= 0)
        {
            return std::make_pair(s1 / (2 * a), s2 / (2 * a));
        }
        else
        {
            return std::make_pair((2 * c) / s2, (2 * c) / s1);
        }
    }

} // namespace gridsynth
