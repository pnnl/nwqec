#pragma once

#include <optional>
#include <vector>
#include <random>
#include <chrono>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <string>
#include <iostream>
#include <cassert>
#include <gmp.h>

#include "types.hpp"
#include "ring.hpp"

// Port of src/diophantine.py (exact algorithmic structure, no simplifications)
// All timeouts are in milliseconds.

namespace gridsynth
{
    // Thread-local RNG for all probabilistic routines (reduces reseeding overhead)
    inline std::mt19937_64 &global_rng()
    {
        static thread_local std::mt19937_64 rng(std::random_device{}());
        return rng;
    }
    // Sentinel: we use std::nullopt to indicate NO_SOLUTION externally. Internally we use a bool flag.
    struct ZOmegaOrNoSolution
    {
        ZOmega value; // valid only if has_value && !no_solution
        bool has_value{false};
        bool no_solution{false};
    };

    // Helper: modular exponentiation a^e mod m (precise, fast via GMP powm)
    inline Integer mod_pow(const Integer &base, const Integer &exp, const Integer &mod)
    {
        if (mod == 0)
            return Integer(0);
        if (mod == 1)
            return Integer(0);
        Integer b_norm = base; // copy so we can normalize
        // Normalize base into [0, mod-1]
        mpz_mod(b_norm.get_mpz_t(), b_norm.get_mpz_t(), mod.get_mpz_t());
        Integer result;
        mpz_powm(result.get_mpz_t(), b_norm.get_mpz_t(), exp.get_mpz_t(), mod.get_mpz_t());
        return result;
    }

    // Fast digit count in base 10 for GMP-backed Integer
    inline size_t num_decimal_digits(const Integer &n)
    {
        return static_cast<size_t>(mpz_sizeinbase(n.get_mpz_t(), 10));
    }

    inline Integer _rand_between(Integer low_inclusive, Integer high_inclusive, std::mt19937_64 &rng)
    {
        // Convert to long long for standard distribution, assuming the range fits

        long long low_ll = static_cast<long long>(low_inclusive);
        long long high_ll = static_cast<long long>(high_inclusive);
        std::uniform_int_distribution<long long> dist(low_ll, high_ll);
        return Integer(dist(rng));
    }

    // Pollard-Brent like factor finder (probabilistic); returns optional factor
    inline std::optional<Integer> _find_factor(Integer n, int factoring_timeout_ms, Integer M = 128)
    {
        if (!(n & 1LL) && n > 2)
            return 2; // even
        if (n <= 3)
            return std::nullopt;
        auto &rng = global_rng();
        Integer a = _rand_between(1, n - 1, rng);
        Integer y = a;
        Integer r = 1;
        Integer k = 0;
        // L heuristic copied from Python code but avoid floating precision drift:
        // L = int(10 ** (digits/4) * 1.1774 + 10). We'll approximate by using double only on small digit count.
        size_t digits = num_decimal_digits(n);
        double pow_term = std::pow(10.0, static_cast<double>(digits) / 4.0);
        Integer L = Integer(static_cast<long long>(pow_term * 1.1774 + 10.0));

        auto start = std::chrono::steady_clock::now();
        while (true)
        {
            Integer x = y + n; // Python uses x = y + n (not reduced mod n)
            while (k < r)
            {
                Integer q = 1;
                Integer y0 = y;
                for (Integer iter = 0; iter < M; ++iter)
                {
                    // y = (y*y + a) % n
                    Integer mul = y * y;
                    Integer y_new = mul % n;
                    y = (y_new + a) % n;
                    Integer diff = x - y; // always >= n - y >= 0 in Python version
                    Integer qmul = q * (diff % n);
                    q = qmul % n;
                    ++k;
                    if (k == r)
                        break;
                }
                Integer g = gcd(q, n);
                if (g != 1)
                {
                    if (g == n)
                    {
                        y = y0;
                        for (Integer iter = 0; iter < M; ++iter)
                        {
                            Integer mul = y * y;
                            Integer y_new = mul % n;
                            y = (y_new + a) % n;
                            Integer diff = x - y;
                            Integer g2 = gcd(diff % n, n);
                            if (g2 != 1)
                            {
                                if (g2 == n)
                                    return std::nullopt;
                                else
                                    return g2;
                            }
                        }
                        return std::nullopt;
                    }
                    else
                    {
                        return g;
                    }
                }
                auto now = std::chrono::steady_clock::now();
                if (k >= L || std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count() >= factoring_timeout_ms)
                {
                    return std::nullopt;
                }
            }
            r <<= 1;
        }
    }

    // Random square root of -1 mod p for p prime p ≡ 1 (mod 4)
    inline std::optional<Integer> _sqrt_negative_one(Integer p, Integer L = 100)
    {
        if (p <= 2)
            return std::nullopt;
        auto &rng = global_rng();
        for (Integer i = 0; i < L; ++i)
        {
            Integer b = _rand_between(1, p - 1, rng);
            Integer h = mod_pow(b, (p - 1) >> 2, p);
            Integer r = (h * h) % p;
            if (r == p - 1)
                return h;
            else if (r != 1)
                return std::nullopt; // fail fast as python returns None
        }
        return std::nullopt;
    }

    class F_p2
    {
    public:
        static Integer base; // element with base^((p-1)/2) ≡ -1 ?
        static Integer p;
        Integer _a; // a + b * x where x^2 = base
        Integer _b;
        F_p2(Integer a = 0, Integer b = 0) : _a((a % p + p) % p), _b((b % p + p) % p) {}
        Integer a() const { return _a; }
        Integer b() const { return _b; }
        F_p2 operator*(const F_p2 &o) const
        {
            // new_a = a1*a2 + ((b1*b2) % p) * base (mod p)
            Integer prod_a = _a * o._a;
            Integer prod_b = _b * o._b % p;
            Integer new_a = (prod_a + prod_b * base) % p;
            Integer new_b = (_a * o._b + _b * o._a) % p;
            return F_p2(new_a, new_b);
        }
        F_p2 pow(Integer e) const
        {
            if (e < 0)
                throw std::invalid_argument("negative exponent");
            F_p2 res(1, 0), tmp = *this;
            while (e > 0)
            {
                if (e & 1)
                    res = res * tmp;
                tmp = tmp * tmp;
                e >>= 1;
            }
            return res;
        }
    };
    inline Integer F_p2::base = 0;
    inline Integer F_p2::p = 0;

    inline std::optional<Integer> _root_mod(Integer x, Integer p, Integer L = 100)
    {
        x %= p;
        if (x < 0)
            x += p;
        if (p == 2)
            return x;
        if (x == 0)
            return 0LL;
        if (!(p & 1LL) && p > 2)
            return std::nullopt;
        // Euler criterion
        Integer t = mod_pow(x, (p - 1) / 2, p);
        if (t != 1)
            return std::nullopt;
        auto &rng = global_rng();
        for (Integer iter = 0; iter < L; ++iter)
        {
            Integer b = _rand_between(1, p - 1, rng);
            // r = b^{p-1} mod p
            Integer r = mod_pow(b, p - 1, p);
            if (r != 1)
                return std::nullopt; // as python does
            Integer candidate_base = ((b * b + p - x) % p);
            // candidate_base^{(p-1)/2} must !=1
            Integer check = mod_pow(candidate_base, (p - 1) / 2, p);
            if (check != 1)
            {
                F_p2::p = p;
                F_p2::base = candidate_base;
                Integer power = (p + 1) / 2;
                F_p2 elem(b, 1);
                F_p2 rfp = elem.pow(power);
                return rfp.a();
            }
        }
        return std::nullopt;
    }

    inline bool _is_prime(Integer n, Integer L = 4)
    {
        if (n < 0)
            n = -n;
        if (n == 0 || n == 1)
            return false;
        if (!(n & 1LL))
            return n == 2;
        Integer d = n - 1;
        Integer r = 0;
        while (!(d & 1LL))
        {
            d >>= 1;
            ++r;
        }
        auto &rng = global_rng();
        for (Integer i = 0; i < L; ++i)
        {
            Integer a = _rand_between(1, n - 1, rng);
            Integer x = mod_pow(a, d, n);
            if (x == 1)
                return true; // matches python's early return (note: python code seems flawed; we keep behavior)
            bool passed = false;
            for (Integer j = 0; j < r; ++j)
            {
                if (x == n - 1)
                {
                    passed = true;
                    break;
                }
                x = x * x % n;
            }
            if (passed)
                return true;
        }
        return false;
    }

    // Decompose relatively integer prime factors - identical logic
    inline std::pair<Integer, std::vector<std::pair<Integer, Integer>>> _decompose_relatively_int_prime(std::vector<std::pair<Integer, Integer>> partial_facs)
    {
        Integer u = 1;
        std::vector<std::pair<Integer, Integer>> stack(partial_facs.rbegin(), partial_facs.rend());
        std::vector<std::pair<Integer, Integer>> facs;
        facs.reserve(partial_facs.size());
        while (!stack.empty())
        {
            auto [b, k_b] = stack.back();
            stack.pop_back();
            size_t i = 0;
            while (true)
            {
                if (i >= facs.size())
                {
                    if (b == 1 || b == -1)
                    {
                        if (b == -1 && (k_b & 1))
                            u = -u;
                    }
                    else
                        facs.emplace_back(b, k_b);
                    break;
                }
                auto &ak = facs[i];
                Integer a = ak.first;
                Integer k_a = ak.second;
                if (a == b || a == -b)
                {
                    if (a == -b && (k_b & 1))
                        u = -u;
                    ak.second = k_a + k_b;
                    break;
                }
                else
                {
                    Integer g = gcd(a, b);
                    if (g == 1 || g == -1)
                    {
                        ++i;
                        continue;
                    }
                    else
                    {
                        std::vector<std::pair<Integer, Integer>> new_partial = {{a / g, k_a}, {g, k_a + k_b}};
                        auto [u_a, facs_a] = _decompose_relatively_int_prime(new_partial);
                        u *= u_a;
                        facs[i] = facs_a[0];
                        facs.insert(facs.end(), facs_a.begin() + 1, facs_a.end());
                        stack.emplace_back(b / g, k_b);
                        break;
                    }
                }
            }
        }
        return {u, facs};
    }

    inline ZOmegaOrNoSolution _adj_decompose_int_prime(Integer p)
    {
        if (p < 0)
            p = -p;
        if (p == 0 || p == 1)
            return {ZOmega::from_int(p), true, false};
        if (p == 2)
            return {ZOmega(-1, 0, 1, 0), true, false};
        if (_is_prime(p))
        {
            if ((p & 0b11) == 1)
            {
                auto h = _sqrt_negative_one(p);
                if (!h)
                    return {ZOmega(), false, false};
                // Python: t = gcd(h + ZOmega(0,1,0,0), p)
                auto t = ZOmega::gcd(ZOmega(0, 1, 0, 0) + ZOmega(0, 0, 0, *h), ZOmega::from_int(p));
                if (t.conj() * t == ZOmega::from_int(p) || t.conj() * t == ZOmega::from_int(-p))
                    return {t, true, false};
                else
                    return {ZOmega(), false, false};
            }
            else if ((p & 0b111) == 3)
            {
                auto h = _root_mod(-2, p);
                if (!h)
                    return {ZOmega(), false, false};
                // Python: t = gcd(h + ZOmega(1,0,1,0), p)
                auto t = ZOmega::gcd(ZOmega(1, 0, 1, 0) + ZOmega(0, 0, 0, *h), ZOmega::from_int(p));
                if (t.conj() * t == ZOmega::from_int(p) || t.conj() * t == ZOmega::from_int(-p))
                    return {t, true, false};
                else
                    return {ZOmega(), false, false};
            }
            else if ((p & 0b111) == 7)
            {
                auto h = _root_mod(2, p);
                if (h)
                    return {ZOmega(), false, true};
                else
                    return {ZOmega(), false, false};
            }
            else
            {
                return {ZOmega(), false, false};
            }
        }
        else
        {
            if ((p & 0b111) == 7)
            {
                auto h = _root_mod(2, p);
                if (h)
                    return {ZOmega(), false, true};
                else
                    return {ZOmega(), false, false};
            }
            else
            {
                return {ZOmega(), false, false};
            }
        }
    }

    inline ZOmegaOrNoSolution _adj_decompose_int_prime_power(Integer p, Integer k)
    {
        if (!(k & 1))
        { // even
            // return p ** (k//2)
            Integer e = k / 2;
            ZOmega z = ZOmega::from_int(1);
            ZOmega base = ZOmega::from_int(p);
            while (e > 0)
            {
                if (e & 1)
                    z = z * base;
                base = base * base;
                e >>= 1;
            }
            return {z, true, false};
        }
        else
        {
            auto t = _adj_decompose_int_prime(p);
            if (!t.has_value || t.no_solution)
                return t;
            // t ** k
            Integer e = k - 1;
            ZOmega acc = t.value;
            ZOmega base = t.value;
            while (e > 0)
            {
                if (e & 1)
                    acc = acc * base;
                base = base * base;
                e >>= 1;
            }
            return {acc, true, false};
        }
    }

    inline ZOmegaOrNoSolution _adj_decompose_int(Integer n, int diophantine_timeout_ms, int factoring_timeout_ms, std::chrono::steady_clock::time_point start)
    {
        if (n < 0)
            n = -n;
        std::vector<std::pair<Integer, Integer>> facs = {{n, 1}};
        ZOmega t = ZOmega::from_int(1);
        while (!facs.empty())
        {
            auto [p, k] = facs.back();
            facs.pop_back();
            auto t_p = _adj_decompose_int_prime_power(p, k);
            if (t_p.no_solution)
                return {ZOmega(), false, true};
            else if (!t_p.has_value)
            {
                auto fac = _find_factor(p, factoring_timeout_ms);
                if (!fac)
                {
                    facs.emplace_back(p, k);
                    auto now = std::chrono::steady_clock::now();
                    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
                    if (elapsed >= diophantine_timeout_ms)
                    {
                        return {ZOmega(), false, true};
                    }
                }
                else
                {
                    facs.emplace_back(p / *fac, k);
                    facs.emplace_back(*fac, k);
                    auto decomposed = _decompose_relatively_int_prime(facs);
                    facs = decomposed.second;
                }
            }
            else
            {
                t = t * t_p.value;
            }
        }
        return {t, true, false};
    }

    inline ZOmegaOrNoSolution _adj_decompose_selfassociate(const ZRootTwo &xi, int diophantine_timeout_ms, int factoring_timeout_ms, std::chrono::steady_clock::time_point start)
    {
        if (xi == ZRootTwo::from_int(0))
            return {ZOmega::from_int(0), true, false};
        Integer n = gcd(xi.a(), xi.b());
        ZRootTwo r = xi / ZRootTwo::from_int(n);
        auto t1 = _adj_decompose_int(n, diophantine_timeout_ms, factoring_timeout_ms, start);
        ZOmega t2 = ((r % ZRootTwo(0, 1)) == ZRootTwo::from_int(0) ? ZOmega(0, 0, 1, 1) : ZOmega::from_int(1));
        if (!t1.has_value)
            return t1; // includes no_solution flag
        if (t1.no_solution)
            return t1;
        return {t1.value * t2, true, false};
    }

    // Relatively ZOmega prime decomposition (ZRootTwo factors) NOTE: this part may require ZRootTwo.gcd available (it is)
    inline std::pair<ZRootTwo, std::vector<std::pair<ZRootTwo, Integer>>> _decompose_relatively_zomega_prime(std::vector<std::pair<ZRootTwo, Integer>> partial)
    {
        ZRootTwo u = ZRootTwo::from_int(1);
        std::vector<std::pair<ZRootTwo, Integer>> stack(partial.rbegin(), partial.rend());
        std::vector<std::pair<ZRootTwo, Integer>> facs;
        facs.reserve(partial.size());
        while (!stack.empty())
        {
            auto [b, k_b] = stack.back();
            stack.pop_back();
            size_t i = 0;
            while (true)
            {
                if (i >= facs.size())
                {
                    if (ZRootTwo::sim(b, ZRootTwo::from_int(1)))
                    {
                        for (Integer i2 = 0; i2 < k_b; ++i2)
                            u = u * b;
                    }
                    else
                        facs.emplace_back(b, k_b);
                    break;
                }
                auto &ak = facs[i];
                auto a = ak.first;
                auto k_a = ak.second;
                if (ZRootTwo::sim(a, b))
                { // b // a**k_b times
                    // units accumulate: (b // a)^k_b
                    ZRootTwo quotient = b / a; // may be unit
                    for (Integer i2 = 0; i2 < k_b; ++i2)
                        u = u * quotient;
                    ak.second = k_a + k_b;
                    break;
                }
                else
                {
                    ZRootTwo g = ZRootTwo::gcd(a, b);
                    if (ZRootTwo::sim(g, ZRootTwo::from_int(1)))
                    {
                        ++i;
                        continue;
                    }
                    else
                    {
                        std::vector<std::pair<ZRootTwo, Integer>> new_partial = {{a / g, k_a}, {g, k_a + k_b}};
                        auto [u_a, facs_a] = _decompose_relatively_zomega_prime(new_partial);
                        u = u * u_a;
                        facs[i] = facs_a[0];
                        facs.insert(facs.end(), facs_a.begin() + 1, facs_a.end());
                        stack.emplace_back(b / g, k_b);
                        break;
                    }
                }
            }
        }
        return {u, facs};
    }

    inline ZOmegaOrNoSolution _adj_decompose_zomega_prime(const ZRootTwo &eta)
    {
        Integer p = eta.norm();
        if (p < 0)
            p = -p;
        if (p == 0 || p == 1)
            return {ZOmega::from_int(p), true, false};
        if (p == 2)
            return {ZOmega(-1, 0, 1, 0), true, false};
        if (_is_prime(p))
        {
            if ((p & 0b11) == 1)
            {
                auto h = _sqrt_negative_one(p);
                if (!h)
                {
                    return {ZOmega(), false, false};
                }
                // Python: t = gcd(h + ZOmega(0,1,0,0), eta)
                auto t = ZOmega::gcd(ZOmega(0, 1, 0, 0) + ZOmega(0, 0, 0, *h), ZOmega::from_zroottwo(eta));

                // Need sim(t.conj * t, eta)
                if (ZRootTwo::sim(ZRootTwo::from_zomega(t.conj() * t), eta))
                    return {t, true, false};
                else
                    return {ZOmega(), false, false};
            }
            else if ((p & 0b111) == 3)
            {
                auto h = _root_mod(-2, p);
                if (!h)
                    return {ZOmega(), false, false};
                // Python: t = gcd(h + ZOmega(1,0,1,0), eta)
                auto t = ZOmega::gcd(ZOmega(1, 0, 1, 0) + ZOmega(0, 0, 0, *h), ZOmega::from_zroottwo(eta));
                if (ZRootTwo::sim(ZRootTwo::from_zomega(t.conj() * t), eta))
                    return {t, true, false};
                else
                    return {ZOmega(), false, false};
            }
            else if ((p & 0b111) == 7)
            {
                auto h = _root_mod(2, p);
                if (h)
                    return {ZOmega(), false, true};
                else
                    return {ZOmega(), false, false};
            }
            else
                return {ZOmega(), false, false};
        }
        else
        {

            if ((p & 0b111) == 7)
            {
                auto h = _root_mod(2, p);
                if (h)
                    return {ZOmega(), false, true};
                else
                    return {ZOmega(), false, false};
            }
            else
                return {ZOmega(), false, false};
        }
    }

    inline ZOmegaOrNoSolution _adj_decompose_zomega_prime_power(const ZRootTwo &eta, Integer k)
    {
        if (!(k & 1))
        {
            Integer e = k / 2;
            ZRootTwo eta_pow = ZRootTwo::from_int(1);
            ZRootTwo base = eta;
            while (e > 0)
            {
                if (e & 1)
                    eta_pow = eta_pow * base;
                base = base * base;
                e >>= 1;
            }
            // convert to ZOmega
            ZOmega zeta = ZOmega::from_zroottwo(eta_pow);
            return {zeta, true, false};
        }
        else
        {
            auto t = _adj_decompose_zomega_prime(eta);

            if (!t.has_value || t.no_solution)
                return t;
            Integer e = k - 1;
            ZOmega acc = t.value;
            ZOmega base = t.value;
            while (e > 0)
            {
                if (e & 1)
                    acc = acc * base;
                base = base * base;
                e >>= 1;
            }
            return {acc, true, false};
        }
    }

    inline ZOmegaOrNoSolution _adj_decompose_selfcoprime(const ZRootTwo &xi, int diophantine_timeout_ms, int factoring_timeout_ms, std::chrono::steady_clock::time_point start)
    {
        std::vector<std::pair<ZRootTwo, Integer>> facs = {{xi, 1}};
        ZOmega t = ZOmega::from_int(1);
        while (!facs.empty())
        {
            auto [eta, k] = facs.back();
            facs.pop_back();
            auto t_eta = _adj_decompose_zomega_prime_power(eta, k);
            if (t_eta.no_solution)
                return {ZOmega(), false, true};
            else if (!t_eta.has_value)
            {
                Integer n = eta.norm();
                if (n < 0)
                    n = -n;
                auto fac_n = _find_factor(n, factoring_timeout_ms);
                if (!fac_n)
                {
                    facs.emplace_back(eta, k);
                    auto now = std::chrono::steady_clock::now();
                    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
                    if (elapsed >= diophantine_timeout_ms)
                    {
                        return {ZOmega(), false, true};
                    }
                }
                else
                {
                    ZRootTwo fac = ZRootTwo::gcd(xi, ZRootTwo::from_int(*fac_n));
                    facs.emplace_back(eta / fac, k);
                    facs.emplace_back(fac, k);
                    auto decomposed = _decompose_relatively_zomega_prime(facs);
                    facs = decomposed.second;
                }
            }
            else
            {
                t = t * t_eta.value;
            }
        }
        return {t, true, false};
    }

    inline ZOmegaOrNoSolution _adj_decompose(const ZRootTwo &xi, int diophantine_timeout_ms, int factoring_timeout_ms, std::chrono::steady_clock::time_point start)
    {

        if (xi == ZRootTwo::from_int(0))
        {
            return {ZOmega::from_int(0), true, false};
        }

        ZRootTwo xi_conj_sq2 = xi.conj_sq2();

        ZRootTwo d = ZRootTwo::gcd(xi, xi_conj_sq2);

        ZRootTwo eta = xi / d;
        auto t1 = _adj_decompose_selfassociate(d, diophantine_timeout_ms, factoring_timeout_ms, start);
        if (t1.no_solution)
        {
            return t1;
        }

        auto t2 = _adj_decompose_selfcoprime(eta, diophantine_timeout_ms, factoring_timeout_ms, start);
        if (t2.no_solution)
        {
            return t2;
        }

        if (!t1.has_value || !t2.has_value)
            return {ZOmega(), false, t1.no_solution || t2.no_solution};

        return {t1.value * t2.value, true, false};
    }

    inline ZOmegaOrNoSolution _diophantine(const ZRootTwo &xi, int diophantine_timeout_ms, int factoring_timeout_ms, std::chrono::steady_clock::time_point start)
    {
        if (xi == ZRootTwo::from_int(0))
        {
            return {ZOmega::from_int(0), true, false};
        }
        if (xi < ZRootTwo::from_int(0) || xi.conj_sq2() < ZRootTwo::from_int(0))
        {
            return {ZOmega(), false, true};
        }

        auto t = _adj_decompose(xi, diophantine_timeout_ms, factoring_timeout_ms, start);
        if (t.no_solution || !t.has_value)
        {
            return t;
        }

        ZRootTwo xi_associate = ZRootTwo::from_zomega(t.value.conj() * t.value);
        ZRootTwo u = xi / xi_associate;
        auto v_opt = u.sqrt();
        if (!v_opt)
        {
            return {ZOmega(), false, true};
        }
        ZOmega v_zomega = ZOmega::from_zroottwo(*v_opt);
        ZOmega result = v_zomega * t.value;
        return {result, true, false};
    }

    inline std::optional<DOmega> diophantine_dyadic(const DRootTwo &xi, int diophantine_timeout, int factoring_timeout)
    {
        Integer k_div_2 = xi.k() >> 1;
        Integer k_mod_2 = xi.k() & 1;
        ZRootTwo arg = k_mod_2 ? (xi.alpha() * ZRootTwo(1, 1)) : xi.alpha();

        auto start = std::chrono::steady_clock::now();
        auto t = _diophantine(arg, diophantine_timeout, factoring_timeout, start);

        if (!t.has_value || t.no_solution)
        {
            return std::nullopt;
        }
        ZOmega z = t.value;
        if (k_mod_2)
        {
            z = z * ZOmega(0, -1, 1, 0);
        }
        return DOmega(z, k_div_2 + k_mod_2);
    }

} // namespace gridsynth
