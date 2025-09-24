#pragma once
#include <iostream>
#include <string>
#include <type_traits>
#include <gmp.h>

/**
 * gmp_integer.hpp
 *
 * GMP-based arbitrary precision Integer implementation that provides
 * the same interface as the __int128_t-based Integer class.
 */

namespace gridsynth
{
    /**
     * GMP-based Integer class with arbitrary precision arithmetic
     */
    class GMPInteger
    {
    private:
        mpz_t value_;

    public:
        // Constructors
        GMPInteger() { mpz_init(value_); }

        GMPInteger(int val) { mpz_init_set_si(value_, val); }

        GMPInteger(long val) { mpz_init_set_si(value_, val); }

        GMPInteger(long long val) { mpz_init_set_si(value_, val); }

        GMPInteger(double val) { mpz_init_set_d(value_, val); }

        // Copy constructor
        GMPInteger(const GMPInteger &other) { mpz_init_set(value_, other.value_); }

        // Move constructor
        GMPInteger(GMPInteger &&other) noexcept
        {
            mpz_init(value_);
            mpz_swap(value_, other.value_);
        }

        // Destructor
        ~GMPInteger() { mpz_clear(value_); }

        // No global pre-allocation; use GMP growth heuristics

        bool is_odd() const
        {
            return mpz_odd_p(value_);
        }

        // Assignment operators
        GMPInteger &operator=(const GMPInteger &other)
        {
            if (this != &other)
            {
                mpz_set(value_, other.value_);
            }
            return *this;
        }

        GMPInteger &operator=(GMPInteger &&other) noexcept
        {
            if (this != &other)
            {
                mpz_swap(value_, other.value_);
            }
            return *this;
        }

        GMPInteger &operator=(int val)
        {
            mpz_set_si(value_, val);
            return *this;
        }

        GMPInteger &operator=(long val)
        {
            mpz_set_si(value_, val);
            return *this;
        }

        GMPInteger &operator=(long long val)
        {
            mpz_set_si(value_, val);
            return *this;
        }

        // Conversion operators
        explicit operator int() const { return static_cast<int>(mpz_get_si(value_)); }

        explicit operator long long() const { return static_cast<long long>(mpz_get_si(value_)); }

        explicit operator double() const { return mpz_get_d(value_); }

        operator size_t() const { return static_cast<size_t>(mpz_get_ui(value_)); } // Implicit for array indexing

        explicit operator bool() const { return mpz_cmp_si(value_, 0) != 0; }

        // Access to internal representation for efficient operations
        mpz_t &get_mpz_t() { return value_; }
        const mpz_t &get_mpz_t() const { return value_; }

        // Arithmetic operators
        GMPInteger operator+(const GMPInteger &other) const
        {
            GMPInteger result;
            mpz_add(result.value_, value_, other.value_);
            return result;
        }

        GMPInteger operator-(const GMPInteger &other) const
        {
            GMPInteger result;
            mpz_sub(result.value_, value_, other.value_);
            return result;
        }

        GMPInteger operator*(const GMPInteger &other) const
        {
            GMPInteger result;
            mpz_mul(result.value_, value_, other.value_);
            return result;
        }

        GMPInteger operator/(const GMPInteger &other) const
        {
            GMPInteger result;
            mpz_tdiv_q(result.value_, value_, other.value_);
            return result;
        }

        GMPInteger operator%(const GMPInteger &other) const
        {
            GMPInteger result;
            mpz_tdiv_r(result.value_, value_, other.value_);
            return result;
        }

        // Compound assignment operators
        GMPInteger &operator+=(const GMPInteger &other)
        {
            mpz_add(value_, value_, other.value_);
            return *this;
        }

        GMPInteger &operator-=(const GMPInteger &other)
        {
            mpz_sub(value_, value_, other.value_);
            return *this;
        }

        GMPInteger &operator*=(const GMPInteger &other)
        {
            mpz_mul(value_, value_, other.value_);
            return *this;
        }

        GMPInteger &operator/=(const GMPInteger &other)
        {
            mpz_tdiv_q(value_, value_, other.value_);
            return *this;
        }

        GMPInteger &operator%=(const GMPInteger &other)
        {
            mpz_tdiv_r(value_, value_, other.value_);
            return *this;
        }

        // Compound assignment with built-in signed types (avoid temporary objects)
        GMPInteger &operator+=(long long rhs)
        {
            if (rhs >= 0)
                mpz_add_ui(value_, value_, static_cast<unsigned long>(rhs));
            else
                mpz_sub_ui(value_, value_, static_cast<unsigned long>(-rhs));
            return *this;
        }
        GMPInteger &operator-=(long long rhs)
        {
            if (rhs >= 0)
                mpz_sub_ui(value_, value_, static_cast<unsigned long>(rhs));
            else
                mpz_add_ui(value_, value_, static_cast<unsigned long>(-rhs));
            return *this;
        }
        GMPInteger &operator*=(long long rhs)
        {
#if defined(__GNU_MP_VERSION)
            // mpz_mul_si is available in GMP
            mpz_mul_si(value_, value_, rhs);
#else
            if (rhs >= 0)
                mpz_mul_ui(value_, value_, static_cast<unsigned long>(rhs));
            else
            {
                mpz_mul_ui(value_, value_, static_cast<unsigned long>(-rhs));
                mpz_neg(value_, value_);
            }
#endif
            return *this;
        }
        GMPInteger &operator/=(long long rhs)
        {
            if (rhs == 0)
                throw std::runtime_error("Division by zero");
            bool neg = rhs < 0;
            unsigned long mag = static_cast<unsigned long>(neg ? -rhs : rhs);
            mpz_tdiv_q_ui(value_, value_, mag);
            if (neg)
                mpz_neg(value_, value_);
            return *this;
        }
        GMPInteger &operator%=(long long rhs)
        {
            if (rhs == 0)
                throw std::runtime_error("Modulo by zero");
            bool neg = rhs < 0; // Sign of remainder follows dividend per C++ semantics, so ignore sign of rhs
            unsigned long mag = static_cast<unsigned long>(neg ? -rhs : rhs);
            mpz_t r;
            mpz_init(r);
            mpz_tdiv_r_ui(r, value_, mag);
            mpz_set(value_, r); // store remainder
            mpz_clear(r);
            return *this;
        }

        GMPInteger &operator+=(int rhs) { return (*this += static_cast<long long>(rhs)); }
        GMPInteger &operator-=(int rhs) { return (*this -= static_cast<long long>(rhs)); }
        GMPInteger &operator*=(int rhs) { return (*this *= static_cast<long long>(rhs)); }
        GMPInteger &operator/=(int rhs) { return (*this /= static_cast<long long>(rhs)); }
        GMPInteger &operator%=(int rhs) { return (*this %= static_cast<long long>(rhs)); }

        // Increment and decrement operators
        GMPInteger &operator++()
        {
            mpz_add_ui(value_, value_, 1);
            return *this;
        }

        GMPInteger operator++(int)
        {
            GMPInteger temp(*this);
            mpz_add_ui(value_, value_, 1);
            return temp;
        }

        // Comparison operators
        bool operator==(const GMPInteger &other) const { return mpz_cmp(value_, other.value_) == 0; }
        bool operator!=(const GMPInteger &other) const { return mpz_cmp(value_, other.value_) != 0; }
        bool operator<(const GMPInteger &other) const { return mpz_cmp(value_, other.value_) < 0; }
        bool operator<=(const GMPInteger &other) const { return mpz_cmp(value_, other.value_) <= 0; }
        bool operator>(const GMPInteger &other) const { return mpz_cmp(value_, other.value_) > 0; }
        bool operator>=(const GMPInteger &other) const { return mpz_cmp(value_, other.value_) >= 0; }

        // Bitwise operators
        GMPInteger operator<<(int shift) const
        {
            GMPInteger result;
            mpz_mul_2exp(result.value_, value_, shift);
            return result;
        }

        GMPInteger operator>>(int shift) const
        {
            GMPInteger result;
            mpz_tdiv_q_2exp(result.value_, value_, shift);
            return result;
        }

        GMPInteger operator<<(const GMPInteger &shift) const
        {
            GMPInteger result;
            mpz_mul_2exp(result.value_, value_, static_cast<mp_bitcnt_t>(mpz_get_ui(shift.value_)));
            return result;
        }

        GMPInteger operator>>(const GMPInteger &shift) const
        {
            GMPInteger result;
            mpz_tdiv_q_2exp(result.value_, value_, static_cast<mp_bitcnt_t>(mpz_get_ui(shift.value_)));
            return result;
        }

        GMPInteger operator&(const GMPInteger &other) const
        {
            GMPInteger result;
            mpz_and(result.value_, value_, other.value_);
            return result;
        }

        GMPInteger operator|(const GMPInteger &other) const
        {
            GMPInteger result;
            mpz_ior(result.value_, value_, other.value_);
            return result;
        }

        GMPInteger operator^(const GMPInteger &other) const
        {
            GMPInteger result;
            mpz_xor(result.value_, value_, other.value_);
            return result;
        }

        // Compound bitwise assignment operators
        GMPInteger &operator<<=(int shift)
        {
            mpz_mul_2exp(value_, value_, shift);
            return *this;
        }

        GMPInteger &operator>>=(int shift)
        {
            mpz_tdiv_q_2exp(value_, value_, shift);
            return *this;
        }

        GMPInteger &operator&=(const GMPInteger &other)
        {
            mpz_and(value_, value_, other.value_);
            return *this;
        }

        // Unary operators
        GMPInteger operator-() const
        {
            GMPInteger result;
            mpz_neg(result.value_, value_);
            return result;
        }

        bool operator!() const { return mpz_cmp_si(value_, 0) == 0; }

        // Stream output operator
        friend std::ostream &operator<<(std::ostream &os, const GMPInteger &val)
        {
            char *str = mpz_get_str(nullptr, 10, val.value_);
            os << str;
            free(str);
            return os;
        }

        GMPInteger floorsqrt() const
        {
            GMPInteger result;
            mpz_sqrt(result.value_, value_);
            return result;
        }
    };

    // Optimized mixed operations with long long (avoid temporary GMPInteger construction)

    inline GMPInteger operator+(const GMPInteger &lhs, long long rhs)
    {
        GMPInteger result(lhs);
        result += rhs;
        return result;
    }
    inline GMPInteger operator+(long long lhs, const GMPInteger &rhs)
    {
        GMPInteger result(rhs);
        result += lhs;
        return result;
    }
    inline GMPInteger operator-(const GMPInteger &lhs, long long rhs)
    {
        GMPInteger result(lhs);
        result -= rhs;
        return result;
    }
    inline GMPInteger operator-(long long lhs, const GMPInteger &rhs)
    {
        GMPInteger result(lhs);
        result -= static_cast<long long>(rhs.operator long long());
        return result;
    }
    inline GMPInteger operator*(const GMPInteger &lhs, long long rhs)
    {
        GMPInteger result(lhs);
        result *= rhs;
        return result;
    }
    inline GMPInteger operator*(long long lhs, const GMPInteger &rhs)
    {
        GMPInteger result(rhs);
        result *= lhs;
        return result;
    }
    inline GMPInteger operator/(const GMPInteger &lhs, long long rhs)
    {
        GMPInteger result(lhs);
        result /= rhs;
        return result;
    }
    inline GMPInteger operator/(long long lhs, const GMPInteger &rhs)
    {
        GMPInteger result(lhs);
        // Division by large rhs still needs full mpz operation
        mpz_tdiv_q(result.get_mpz_t(), result.get_mpz_t(), rhs.get_mpz_t());
        return result;
    }
    inline GMPInteger operator%(const GMPInteger &lhs, long long rhs)
    {
        GMPInteger result(lhs);
        result %= rhs;
        return result;
    }
    inline GMPInteger operator%(long long lhs, const GMPInteger &rhs)
    {
        GMPInteger result(lhs);
        mpz_tdiv_r(result.get_mpz_t(), result.get_mpz_t(), rhs.get_mpz_t());
        return result;
    }
    inline GMPInteger operator<<(long long lhs, const GMPInteger &rhs)
    {
        GMPInteger temp(lhs);
        return temp << static_cast<int>(rhs); // potential narrowing
    }

    // Comparison operators with long long
    inline bool operator==(const GMPInteger &lhs, long long rhs) { return lhs == GMPInteger(rhs); }
    inline bool operator!=(const GMPInteger &lhs, long long rhs) { return lhs != GMPInteger(rhs); }
    inline bool operator<(const GMPInteger &lhs, long long rhs) { return lhs < GMPInteger(rhs); }
    inline bool operator<=(const GMPInteger &lhs, long long rhs) { return lhs <= GMPInteger(rhs); }
    inline bool operator>(const GMPInteger &lhs, long long rhs) { return lhs > GMPInteger(rhs); }
    inline bool operator>=(const GMPInteger &lhs, long long rhs) { return lhs >= GMPInteger(rhs); }

    inline bool operator==(long long lhs, const GMPInteger &rhs) { return GMPInteger(lhs) == rhs; }
    inline bool operator!=(long long lhs, const GMPInteger &rhs) { return GMPInteger(lhs) != rhs; }
    inline bool operator<(long long lhs, const GMPInteger &rhs) { return GMPInteger(lhs) < rhs; }
    inline bool operator<=(long long lhs, const GMPInteger &rhs) { return GMPInteger(lhs) <= rhs; }
    inline bool operator>(long long lhs, const GMPInteger &rhs) { return GMPInteger(lhs) > rhs; }
    inline bool operator>=(long long lhs, const GMPInteger &rhs) { return GMPInteger(lhs) >= rhs; }

    // Bitwise operators with long long
    inline GMPInteger operator&(const GMPInteger &lhs, long long rhs) { return lhs & GMPInteger(rhs); }
    inline GMPInteger operator|(const GMPInteger &lhs, long long rhs) { return lhs | GMPInteger(rhs); }
    inline GMPInteger operator^(const GMPInteger &lhs, long long rhs) { return lhs ^ GMPInteger(rhs); }

    inline GMPInteger operator&(long long lhs, const GMPInteger &rhs) { return GMPInteger(lhs) & rhs; }
    inline GMPInteger operator|(long long lhs, const GMPInteger &rhs) { return GMPInteger(lhs) | rhs; }
    inline GMPInteger operator^(long long lhs, const GMPInteger &rhs) { return GMPInteger(lhs) ^ rhs; }

    // Mixed operations with int types
    // Optimized mixed operations with int (delegate to long long overloads)
    inline GMPInteger operator+(const GMPInteger &lhs, int rhs) { return lhs + static_cast<long long>(rhs); }
    inline GMPInteger operator+(int lhs, const GMPInteger &rhs) { return static_cast<long long>(lhs) + rhs; }
    inline GMPInteger operator-(const GMPInteger &lhs, int rhs) { return lhs - static_cast<long long>(rhs); }
    inline GMPInteger operator-(int lhs, const GMPInteger &rhs) { return static_cast<long long>(lhs) - rhs; }
    inline GMPInteger operator*(const GMPInteger &lhs, int rhs) { return lhs * static_cast<long long>(rhs); }
    inline GMPInteger operator*(int lhs, const GMPInteger &rhs) { return static_cast<long long>(lhs) * rhs; }
    inline GMPInteger operator/(const GMPInteger &lhs, int rhs) { return lhs / static_cast<long long>(rhs); }
    inline GMPInteger operator/(int lhs, const GMPInteger &rhs) { return static_cast<long long>(lhs) / rhs; }
    inline GMPInteger operator%(const GMPInteger &lhs, int rhs) { return lhs % static_cast<long long>(rhs); }
    inline GMPInteger operator%(int lhs, const GMPInteger &rhs) { return static_cast<long long>(lhs) % rhs; }

    // Comparison operators with int
    inline bool operator==(const GMPInteger &lhs, int rhs) { return lhs == GMPInteger(rhs); }
    inline bool operator!=(const GMPInteger &lhs, int rhs) { return lhs != GMPInteger(rhs); }
    inline bool operator<(const GMPInteger &lhs, int rhs) { return lhs < GMPInteger(rhs); }
    inline bool operator<=(const GMPInteger &lhs, int rhs) { return lhs <= GMPInteger(rhs); }
    inline bool operator>(const GMPInteger &lhs, int rhs) { return lhs > GMPInteger(rhs); }
    inline bool operator>=(const GMPInteger &lhs, int rhs) { return lhs >= GMPInteger(rhs); }

    inline bool operator==(int lhs, const GMPInteger &rhs) { return GMPInteger(lhs) == rhs; }
    inline bool operator!=(int lhs, const GMPInteger &rhs) { return GMPInteger(lhs) != rhs; }
    inline bool operator<(int lhs, const GMPInteger &rhs) { return GMPInteger(lhs) < rhs; }
    inline bool operator<=(int lhs, const GMPInteger &rhs) { return GMPInteger(lhs) <= rhs; }
    inline bool operator>(int lhs, const GMPInteger &rhs) { return GMPInteger(lhs) > rhs; }
    inline bool operator>=(int lhs, const GMPInteger &rhs) { return GMPInteger(lhs) >= rhs; }

    // Bitwise operators with int
    inline GMPInteger operator&(const GMPInteger &lhs, int rhs) { return lhs & GMPInteger(rhs); }
    inline GMPInteger operator|(const GMPInteger &lhs, int rhs) { return lhs | GMPInteger(rhs); }
    inline GMPInteger operator^(const GMPInteger &lhs, int rhs) { return lhs ^ GMPInteger(rhs); }

    inline GMPInteger operator&(int lhs, const GMPInteger &rhs) { return GMPInteger(lhs) & rhs; }
    inline GMPInteger operator|(int lhs, const GMPInteger &rhs) { return GMPInteger(lhs) | rhs; }
    inline GMPInteger operator^(int lhs, const GMPInteger &rhs) { return GMPInteger(lhs) ^ rhs; }

    // Mixed operations with double (result is double to preserve precision)
    inline double operator+(const GMPInteger &lhs, double rhs) { return static_cast<double>(lhs) + rhs; }
    inline double operator-(const GMPInteger &lhs, double rhs) { return static_cast<double>(lhs) - rhs; }
    inline double operator*(const GMPInteger &lhs, double rhs) { return static_cast<double>(lhs) * rhs; }
    inline double operator/(const GMPInteger &lhs, double rhs) { return static_cast<double>(lhs) / rhs; }

    inline double operator+(double lhs, const GMPInteger &rhs) { return lhs + static_cast<double>(rhs); }
    inline double operator-(double lhs, const GMPInteger &rhs) { return lhs - static_cast<double>(rhs); }
    inline double operator*(double lhs, const GMPInteger &rhs) { return lhs * static_cast<double>(rhs); }
    inline double operator/(double lhs, const GMPInteger &rhs) { return lhs / static_cast<double>(rhs); }

    // Comparison operators with double
    inline bool operator==(const GMPInteger &lhs, double rhs) { return static_cast<double>(lhs) == rhs; }
    inline bool operator!=(const GMPInteger &lhs, double rhs) { return static_cast<double>(lhs) != rhs; }
    inline bool operator<(const GMPInteger &lhs, double rhs) { return static_cast<double>(lhs) < rhs; }
    inline bool operator<=(const GMPInteger &lhs, double rhs) { return static_cast<double>(lhs) <= rhs; }
    inline bool operator>(const GMPInteger &lhs, double rhs) { return static_cast<double>(lhs) > rhs; }
    inline bool operator>=(const GMPInteger &lhs, double rhs) { return static_cast<double>(lhs) >= rhs; }

    inline bool operator==(double lhs, const GMPInteger &rhs) { return lhs == static_cast<double>(rhs); }
    inline bool operator!=(double lhs, const GMPInteger &rhs) { return lhs != static_cast<double>(rhs); }
    inline bool operator<(double lhs, const GMPInteger &rhs) { return lhs < static_cast<double>(rhs); }
    inline bool operator<=(double lhs, const GMPInteger &rhs) { return lhs <= static_cast<double>(rhs); }
    inline bool operator>(double lhs, const GMPInteger &rhs) { return lhs > static_cast<double>(rhs); }
    inline bool operator>=(double lhs, const GMPInteger &rhs) { return lhs >= static_cast<double>(rhs); }

} // namespace gridsynth
