#pragma once
#include <iostream>
#include <string>
#include <type_traits>
#include <cmath>
#include <mpfr.h>
#include <cctype>
// <algorithm> not needed after refactor

#include "gmp_integer.hpp" // Include GMPInteger for mixed operations>

/**
 * gmp_float.hpp
 *
 * MPFR-based arbitrary precision floating point implementation that provides
 * the same interface as double with arbitrary precision arithmetic.
 */

namespace gridsynth
{
    /**
     * MPFR-based floating point class with arbitrary precision arithmetic
     */
    class GMPFloat
    {
    private:
        mpfr_t value_;
        static mpfr_prec_t default_precision_;
        // Internal tag to construct directly with a specified precision without re-init overhead
        struct direct_init_tag
        {
        };
        GMPFloat(direct_init_tag, mpfr_prec_t prec)
        {
            mpfr_init2(value_, prec);
            mpfr_set_zero(value_, 1);
        }

    public:
        // Static precision management
        static void set_default_precision(mpfr_prec_t prec) { default_precision_ = prec; }
        static mpfr_prec_t get_default_precision() { return default_precision_; }

        // Static factory method for creating GMPFloat with specific precision
        static GMPFloat with_precision(mpfr_prec_t precision, double val = 0.0)
        {
            GMPFloat result(direct_init_tag{}, precision);
            if (val == 0.0)
            {
                // Already zeroed
            }
            else
            {
                mpfr_set_d(result.value_, val, MPFR_RNDN);
            }
            return result;
        }

        // Constructors
        GMPFloat()
        {
            mpfr_init2(value_, default_precision_);
            mpfr_set_zero(value_, 1);
        }

        GMPFloat(int val)
        {
            mpfr_init2(value_, default_precision_);
            mpfr_set_si(value_, val, MPFR_RNDN);
        }

        GMPFloat(long val)
        {
            mpfr_init2(value_, default_precision_);
            mpfr_set_si(value_, val, MPFR_RNDN);
        }

        GMPFloat(long long val)
        {
            mpfr_init2(value_, default_precision_);
            mpfr_set_si(value_, val, MPFR_RNDN);
        }

        GMPFloat(float val)
        {
            mpfr_init2(value_, default_precision_);
            mpfr_set_flt(value_, val, MPFR_RNDN);
        }

        GMPFloat(double val)
        {
            mpfr_init2(value_, default_precision_);
            mpfr_set_d(value_, val, MPFR_RNDN);
        }

        GMPFloat(const GMPInteger &val)
        {
            mpfr_init2(value_, default_precision_);
            mpfr_set_z(value_, val.get_mpz_t(), MPFR_RNDN);
        }

        GMPFloat(long double val)
        {
            mpfr_init2(value_, default_precision_);
            mpfr_set_ld(value_, val, MPFR_RNDN);
        }

        // Small helper to parse forms like: [sign][coeff][*]?pi[/denom]
        static bool parse_pi_expr(const std::string &in, mpfr_t out)
        {
            std::string s;
            s.reserve(in.size());
            for (char c : in)
                if (!std::isspace(static_cast<unsigned char>(c)))
                    s.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
            if (s.find("pi") == std::string::npos)
                return false;
            size_t i = 0;
            int sign = 1;
            if (i < s.size() && (s[i] == '+' || s[i] == '-'))
            {
                if (s[i] == '-')
                    sign = -1;
                ++i;
            }
            // coefficient
            long long coeff = 0;
            bool have_coeff = false;
            while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i])))
            {
                have_coeff = true;
                coeff = coeff * 10 + (s[i] - '0');
                ++i;
            }
            if (!have_coeff)
                coeff = 1; // implicit 1
            if (i < s.size() && s[i] == '*')
                ++i; // optional '*'
            if (i + 1 >= s.size() || s[i] != 'p' || s[i + 1] != 'i')
                return false;
            i += 2; // consume pi
            long long denom = 1;
            if (i < s.size())
            {
                if (s[i] != '/')
                    return false; // nothing else allowed except /denom
                ++i;
                if (i == s.size())
                    return false; // need digits
                long long d = 0;
                bool have = false;
                while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i])))
                {
                    have = true;
                    d = d * 10 + (s[i] - '0');
                    ++i;
                }
                if (!have || d == 0)
                    return false;
                denom = d;
            }
            if (i != s.size())
                return false; // extra junk
            // value = sign * coeff * pi / denom
            mpfr_t pi_val;
            mpfr_init2(pi_val, mpfr_get_prec(out));
            mpfr_const_pi(pi_val, MPFR_RNDN);
            mpfr_set_si(out, sign * coeff, MPFR_RNDN);
            mpfr_mul(out, out, pi_val, MPFR_RNDN);
            if (denom != 1)
                mpfr_div_si(out, out, denom, MPFR_RNDN);
            mpfr_clear(pi_val);
            return true;
        }

        GMPFloat(const std::string &str)
        {
            mpfr_init2(value_, default_precision_);
            if (!parse_pi_expr(str, value_))
                mpfr_set_str(value_, str.c_str(), 10, MPFR_RNDN);
        }

        // Copy constructor
        GMPFloat(const GMPFloat &other)
        {
            mpfr_init2(value_, mpfr_get_prec(other.value_));
            mpfr_set(value_, other.value_, MPFR_RNDN);
        }

        // Move constructor
        GMPFloat(GMPFloat &&other) noexcept
        {
            mpfr_init2(value_, default_precision_);
            mpfr_swap(value_, other.value_);
        }

        // Destructor
        ~GMPFloat() { mpfr_clear(value_); }

        // Assignment operators
        GMPFloat &operator=(const GMPFloat &other)
        {
            if (this != &other)
            {
                mpfr_set(value_, other.value_, MPFR_RNDN);
            }
            return *this;
        }

        GMPFloat &operator=(GMPFloat &&other) noexcept
        {
            if (this != &other)
            {
                mpfr_swap(value_, other.value_);
            }
            return *this;
        }

        GMPFloat &operator=(int val)
        {
            mpfr_set_si(value_, val, MPFR_RNDN);
            return *this;
        }

        GMPFloat &operator=(long val)
        {
            mpfr_set_si(value_, val, MPFR_RNDN);
            return *this;
        }

        GMPFloat &operator=(long long val)
        {
            mpfr_set_si(value_, val, MPFR_RNDN);
            return *this;
        }

        GMPFloat &operator=(float val)
        {
            mpfr_set_flt(value_, val, MPFR_RNDN);
            return *this;
        }

        GMPFloat &operator=(double val)
        {
            mpfr_set_d(value_, val, MPFR_RNDN);
            return *this;
        }

        GMPFloat &operator=(long double val)
        {
            mpfr_set_ld(value_, val, MPFR_RNDN);
            return *this;
        }

        // Conversion operators
        explicit operator int() const { return static_cast<int>(mpfr_get_si(value_, MPFR_RNDN)); }

        explicit operator long() const { return mpfr_get_si(value_, MPFR_RNDN); }

        explicit operator long long() const { return static_cast<long long>(mpfr_get_si(value_, MPFR_RNDN)); }

        explicit operator float() const { return mpfr_get_flt(value_, MPFR_RNDN); }

        operator double() const { return mpfr_get_d(value_, MPFR_RNDN); } // Implicit for double compatibility

        explicit operator long double() const { return mpfr_get_ld(value_, MPFR_RNDN); }

        explicit operator bool() const { return !mpfr_zero_p(value_); }

        // Arithmetic operators
        GMPFloat operator+(const GMPFloat &other) const
        {
            mpfr_prec_t prec = std::max(mpfr_get_prec(value_), mpfr_get_prec(other.value_));
            GMPFloat result(direct_init_tag{}, prec);
            mpfr_add(result.value_, value_, other.value_, MPFR_RNDN);
            return result;
        }

        GMPFloat operator-(const GMPFloat &other) const
        {
            mpfr_prec_t prec = std::max(mpfr_get_prec(value_), mpfr_get_prec(other.value_));
            GMPFloat result(direct_init_tag{}, prec);
            mpfr_sub(result.value_, value_, other.value_, MPFR_RNDN);
            return result;
        }

        GMPFloat operator*(const GMPFloat &other) const
        {
            mpfr_prec_t prec = std::max(mpfr_get_prec(value_), mpfr_get_prec(other.value_));
            GMPFloat result(direct_init_tag{}, prec);
            mpfr_mul(result.value_, value_, other.value_, MPFR_RNDN);
            return result;
        }

        GMPFloat operator/(const GMPFloat &other) const
        {
            mpfr_prec_t prec = std::max(mpfr_get_prec(value_), mpfr_get_prec(other.value_));
            GMPFloat result(direct_init_tag{}, prec);
            mpfr_div(result.value_, value_, other.value_, MPFR_RNDN);
            return result;
        }

        // Compound assignment operators
        GMPFloat &operator+=(const GMPFloat &other)
        {
            mpfr_add(value_, value_, other.value_, MPFR_RNDN);
            return *this;
        }

        GMPFloat &operator-=(const GMPFloat &other)
        {
            mpfr_sub(value_, value_, other.value_, MPFR_RNDN);
            return *this;
        }

        GMPFloat &operator*=(const GMPFloat &other)
        {
            mpfr_mul(value_, value_, other.value_, MPFR_RNDN);
            return *this;
        }

        GMPFloat &operator/=(const GMPFloat &other)
        {
            mpfr_div(value_, value_, other.value_, MPFR_RNDN);
            return *this;
        }

        // Compound assignment with double (avoids temporary GMPFloat allocations)
        GMPFloat &operator+=(double rhs)
        {
            mpfr_add_d(value_, value_, rhs, MPFR_RNDN);
            return *this;
        }
        GMPFloat &operator-=(double rhs)
        {
            mpfr_sub_d(value_, value_, rhs, MPFR_RNDN);
            return *this;
        }
        GMPFloat &operator*=(double rhs)
        {
            mpfr_mul_d(value_, value_, rhs, MPFR_RNDN);
            return *this;
        }
        GMPFloat &operator/=(double rhs)
        {
            mpfr_div_d(value_, value_, rhs, MPFR_RNDN);
            return *this;
        }

        // Compound assignment with int (uses mpfr_*_si helpers)
        GMPFloat &operator+=(int rhs)
        {
            mpfr_add_si(value_, value_, static_cast<long>(rhs), MPFR_RNDN);
            return *this;
        }
        GMPFloat &operator-=(int rhs)
        {
            mpfr_sub_si(value_, value_, static_cast<long>(rhs), MPFR_RNDN);
            return *this;
        }
        GMPFloat &operator*=(int rhs)
        {
            mpfr_mul_si(value_, value_, static_cast<long>(rhs), MPFR_RNDN);
            return *this;
        }
        GMPFloat &operator/=(int rhs)
        {
            mpfr_div_si(value_, value_, static_cast<long>(rhs), MPFR_RNDN);
            return *this;
        }

        // Increment and decrement operators
        GMPFloat &operator++()
        {
            mpfr_add_ui(value_, value_, 1, MPFR_RNDN);
            return *this;
        }

        GMPFloat operator++(int)
        {
            GMPFloat temp(*this);
            mpfr_add_ui(value_, value_, 1, MPFR_RNDN);
            return temp;
        }

        GMPFloat &operator--()
        {
            mpfr_sub_ui(value_, value_, 1, MPFR_RNDN);
            return *this;
        }

        GMPFloat operator--(int)
        {
            GMPFloat temp(*this);
            mpfr_sub_ui(value_, value_, 1, MPFR_RNDN);
            return temp;
        }

        // Comparison operators
        bool operator==(const GMPFloat &other) const { return mpfr_equal_p(value_, other.value_) != 0; }
        bool operator!=(const GMPFloat &other) const { return mpfr_equal_p(value_, other.value_) == 0; }
        bool operator<(const GMPFloat &other) const { return mpfr_less_p(value_, other.value_) != 0; }
        bool operator<=(const GMPFloat &other) const { return mpfr_lessequal_p(value_, other.value_) != 0; }
        bool operator>(const GMPFloat &other) const { return mpfr_greater_p(value_, other.value_) != 0; }
        bool operator>=(const GMPFloat &other) const { return mpfr_greaterequal_p(value_, other.value_) != 0; }

        // Unary operators
        GMPFloat operator-() const
        {
            GMPFloat result(direct_init_tag{}, mpfr_get_prec(value_));
            mpfr_neg(result.value_, value_, MPFR_RNDN);
            return result;
        }

        GMPFloat operator+() const
        {
            return *this;
        }

        bool operator!() const { return mpfr_zero_p(value_) != 0; }

        // Mathematical functions (mimicking std::cmath functions)
        GMPFloat abs() const
        {
            GMPFloat result(direct_init_tag{}, mpfr_get_prec(value_));
            mpfr_abs(result.value_, value_, MPFR_RNDN);
            return result;
        }

        GMPFloat sqrt() const
        {
            GMPFloat result(direct_init_tag{}, mpfr_get_prec(value_));
            mpfr_sqrt(result.value_, value_, MPFR_RNDN);
            return result;
        }

        GMPFloat pow(const GMPFloat &exp) const
        {
            GMPFloat result(direct_init_tag{}, std::max(mpfr_get_prec(value_), mpfr_get_prec(exp.value_)));
            mpfr_pow(result.value_, value_, exp.value_, MPFR_RNDN);
            return result;
        }

        GMPFloat exp() const
        {
            GMPFloat result(direct_init_tag{}, mpfr_get_prec(value_));
            mpfr_exp(result.value_, value_, MPFR_RNDN);
            return result;
        }

        GMPFloat log() const
        {
            GMPFloat result(direct_init_tag{}, mpfr_get_prec(value_));
            mpfr_log(result.value_, value_, MPFR_RNDN);
            return result;
        }

        GMPFloat log10() const
        {
            GMPFloat result(direct_init_tag{}, mpfr_get_prec(value_));
            mpfr_log10(result.value_, value_, MPFR_RNDN);
            return result;
        }

        GMPFloat sin() const
        {
            GMPFloat result(direct_init_tag{}, mpfr_get_prec(value_));
            mpfr_sin(result.value_, value_, MPFR_RNDN);
            return result;
        }

        GMPFloat cos() const
        {
            GMPFloat result(direct_init_tag{}, mpfr_get_prec(value_));
            mpfr_cos(result.value_, value_, MPFR_RNDN);
            return result;
        }

        GMPFloat tan() const
        {
            GMPFloat result(direct_init_tag{}, mpfr_get_prec(value_));
            mpfr_tan(result.value_, value_, MPFR_RNDN);
            return result;
        }

        GMPFloat asin() const
        {
            GMPFloat result(direct_init_tag{}, mpfr_get_prec(value_));
            mpfr_asin(result.value_, value_, MPFR_RNDN);
            return result;
        }

        GMPFloat acos() const
        {
            GMPFloat result(direct_init_tag{}, mpfr_get_prec(value_));
            mpfr_acos(result.value_, value_, MPFR_RNDN);
            return result;
        }

        GMPFloat atan() const
        {
            GMPFloat result(direct_init_tag{}, mpfr_get_prec(value_));
            mpfr_atan(result.value_, value_, MPFR_RNDN);
            return result;
        }

        GMPFloat atan2(const GMPFloat &x) const
        {
            GMPFloat result(direct_init_tag{}, std::max(mpfr_get_prec(value_), mpfr_get_prec(x.value_)));
            mpfr_atan2(result.value_, value_, x.value_, MPFR_RNDN);
            return result;
        }

        GMPFloat sinh() const
        {
            GMPFloat result(direct_init_tag{}, mpfr_get_prec(value_));
            mpfr_sinh(result.value_, value_, MPFR_RNDN);
            return result;
        }

        GMPFloat cosh() const
        {
            GMPFloat result(direct_init_tag{}, mpfr_get_prec(value_));
            mpfr_cosh(result.value_, value_, MPFR_RNDN);
            return result;
        }

        GMPFloat tanh() const
        {
            GMPFloat result(direct_init_tag{}, mpfr_get_prec(value_));
            mpfr_tanh(result.value_, value_, MPFR_RNDN);
            return result;
        }

        GMPFloat floor() const
        {
            GMPFloat result(direct_init_tag{}, mpfr_get_prec(value_));
            mpfr_floor(result.value_, value_);
            return result;
        }

        GMPFloat ceil() const
        {
            GMPFloat result(direct_init_tag{}, mpfr_get_prec(value_));
            mpfr_ceil(result.value_, value_);
            return result;
        }

        GMPFloat round() const
        {
            GMPFloat result(direct_init_tag{}, mpfr_get_prec(value_));
            mpfr_round(result.value_, value_);
            return result;
        }

        // Utility functions
        bool is_nan() const { return mpfr_nan_p(value_) != 0; }
        bool is_inf() const { return mpfr_inf_p(value_) != 0; }
        bool is_zero() const { return mpfr_zero_p(value_) != 0; }
        bool is_finite() const { return mpfr_number_p(value_) != 0; }

        mpfr_prec_t precision() const { return mpfr_get_prec(value_); }
        void set_precision(mpfr_prec_t prec) { mpfr_prec_round(value_, prec, MPFR_RNDN); }

        std::string to_string(int digits = 2) const
        {
            if (mpfr_zero_p(value_))
                return "0.0";

            bool negative = mpfr_sgn(value_) < 0;

            // Get the string representation and exponent
            char *str;
            mpfr_exp_t exp;
            str = mpfr_get_str(nullptr, &exp, 10, digits, value_, MPFR_RNDN);

            std::string mantissa = str;
            mpfr_free_str(str);

            // Remove negative sign from mantissa if present (we'll handle it separately)
            if (mantissa[0] == '-')
                mantissa = mantissa.substr(1);

            std::string result;

            // Python-style automatic formatting logic:
            // Use scientific notation if exponent is too large or too small
            // Use regular notation for reasonable ranges
            int scientific_exp = exp - 1;

            if (scientific_exp >= -4 && scientific_exp < digits)
            {
                // Use regular notation
                if (exp <= 0)
                {
                    // Number is less than 1, e.g., 0.001234
                    result = "0.";
                    for (int i = 0; i < -exp; i++)
                        result += "0";
                    result += mantissa;
                }
                else if (exp >= (int)mantissa.length())
                {
                    // Number is a whole number or needs trailing zeros
                    result = mantissa;
                    for (int i = mantissa.length(); i < exp; i++)
                        result += "0";
                    result += ".0";
                }
                else
                {
                    // Insert decimal point
                    result = mantissa;
                    result.insert(exp, ".");
                }
            }
            else
            {
                // Use scientific notation
                result = mantissa.substr(0, 1);
                if (mantissa.length() > 1)
                {
                    result += ".";
                    result += mantissa.substr(1);
                }
                else
                {
                    result += ".0";
                }

                // Add exponent
                if (scientific_exp >= 0)
                {
                    result += "e+";
                    result += (scientific_exp < 10 ? "0" : "");
                    result += std::to_string(scientific_exp);
                }
                else
                {
                    int abs_exp = -scientific_exp;
                    result += "e-";
                    result += (abs_exp < 10 ? "0" : "");
                    result += std::to_string(abs_exp);
                }
            }

            return negative ? "-" + result : result;
        }

        // Stream output operator
        friend std::ostream &operator<<(std::ostream &os, const GMPFloat &val)
        {
            os << val.to_string();
            return os;
        }

        // Give access to internal mpfr_t for advanced operations
        const mpfr_t &get_mpfr() const { return value_; }
        mpfr_t &get_mpfr() { return value_; }
    };

    // Static member initialization
    mpfr_prec_t GMPFloat::default_precision_ = 256; // Default precision in bits (sufficient for ~1e-50 accuracy)

    // Global mathematical functions (mimicking std namespace functions)
    inline GMPFloat abs(const GMPFloat &x) { return x.abs(); }
    inline GMPFloat sqrt(const GMPFloat &x) { return x.sqrt(); }
    inline GMPFloat pow(const GMPFloat &base, const GMPFloat &exp) { return base.pow(exp); }
    inline GMPFloat exp(const GMPFloat &x) { return x.exp(); }
    inline GMPFloat log(const GMPFloat &x) { return x.log(); }
    inline GMPFloat log10(const GMPFloat &x) { return x.log10(); }
    inline GMPFloat sin(const GMPFloat &x) { return x.sin(); }
    inline GMPFloat cos(const GMPFloat &x) { return x.cos(); }
    inline GMPFloat tan(const GMPFloat &x) { return x.tan(); }
    inline GMPFloat asin(const GMPFloat &x) { return x.asin(); }
    inline GMPFloat acos(const GMPFloat &x) { return x.acos(); }
    inline GMPFloat atan(const GMPFloat &x) { return x.atan(); }
    inline GMPFloat atan2(const GMPFloat &y, const GMPFloat &x) { return y.atan2(x); }
    inline GMPFloat sinh(const GMPFloat &x) { return x.sinh(); }
    inline GMPFloat cosh(const GMPFloat &x) { return x.cosh(); }
    inline GMPFloat tanh(const GMPFloat &x) { return x.tanh(); }
    inline GMPFloat floor(const GMPFloat &x) { return x.floor(); }
    inline GMPFloat ceil(const GMPFloat &x) { return x.ceil(); }
    inline GMPFloat round(const GMPFloat &x) { return x.round(); }

    // Direct GMPFloat to GMPInteger conversions (much faster than via double)
    inline GMPInteger floor_to_gmpinteger(const GMPFloat &x)
    {
        GMPFloat floored = x.floor();
        GMPInteger result;
        mpfr_get_z(result.get_mpz_t(), floored.get_mpfr(), MPFR_RNDN);
        return result;
    }

    inline GMPInteger ceil_to_gmpinteger(const GMPFloat &x)
    {
        GMPFloat ceiled = x.ceil();
        GMPInteger result;
        mpfr_get_z(result.get_mpz_t(), ceiled.get_mpfr(), MPFR_RNDN);
        return result;
    }

    inline GMPInteger round_to_gmpinteger(const GMPFloat &x)
    {
        GMPFloat rounded = x.round();
        GMPInteger result;
        mpfr_get_z(result.get_mpz_t(), rounded.get_mpfr(), MPFR_RNDN);
        return result;
    }

    // Mixed operations with built-in types
    // inline GMPFloat operator+(double lhs, const GMPFloat &rhs) { return GMPFloat(lhs) + rhs; }
    // inline GMPFloat operator-(double lhs, const GMPFloat &rhs) { return GMPFloat(lhs) - rhs; }
    // inline GMPFloat operator*(double lhs, const GMPFloat &rhs) { return GMPFloat(lhs) * rhs; }
    // inline GMPFloat operator/(double lhs, const GMPFloat &rhs) { return GMPFloat(lhs) / rhs; }

    // inline GMPFloat operator+(const GMPFloat &lhs, double rhs) { return lhs + GMPFloat(rhs); }
    // inline GMPFloat operator-(const GMPFloat &lhs, double rhs) { return lhs - GMPFloat(rhs); }
    // inline GMPFloat operator*(const GMPFloat &lhs, double rhs) { return lhs * GMPFloat(rhs); }
    // inline GMPFloat operator/(const GMPFloat &lhs, double rhs) { return lhs / GMPFloat(rhs); }

    // Comparison operators with double
    inline bool operator==(const GMPFloat &lhs, double rhs) { return lhs == GMPFloat(rhs); }
    inline bool operator!=(const GMPFloat &lhs, double rhs) { return lhs != GMPFloat(rhs); }
    inline bool operator<(const GMPFloat &lhs, double rhs) { return lhs < GMPFloat(rhs); }
    inline bool operator<=(const GMPFloat &lhs, double rhs) { return lhs <= GMPFloat(rhs); }
    inline bool operator>(const GMPFloat &lhs, double rhs) { return lhs > GMPFloat(rhs); }
    inline bool operator>=(const GMPFloat &lhs, double rhs) { return lhs >= GMPFloat(rhs); }

    inline bool operator==(double lhs, const GMPFloat &rhs) { return GMPFloat(lhs) == rhs; }
    inline bool operator!=(double lhs, const GMPFloat &rhs) { return GMPFloat(lhs) != rhs; }
    inline bool operator<(double lhs, const GMPFloat &rhs) { return GMPFloat(lhs) < rhs; }
    inline bool operator<=(double lhs, const GMPFloat &rhs) { return GMPFloat(lhs) <= rhs; }
    inline bool operator>(double lhs, const GMPFloat &rhs) { return GMPFloat(lhs) > rhs; }
    inline bool operator>=(double lhs, const GMPFloat &rhs) { return GMPFloat(lhs) >= rhs; }

    // Optimized mixed operations with double (avoid temporary GMPFloat creation)
    inline GMPFloat operator+(const GMPFloat &lhs, double rhs)
    {
        GMPFloat result = GMPFloat::with_precision(lhs.precision());
        mpfr_add_d(result.get_mpfr(), lhs.get_mpfr(), rhs, MPFR_RNDN);
        return result;
    }
    inline GMPFloat operator+(double lhs, const GMPFloat &rhs)
    {
        GMPFloat result = GMPFloat::with_precision(rhs.precision());
        mpfr_add_d(result.get_mpfr(), rhs.get_mpfr(), lhs, MPFR_RNDN);
        return result;
    }
    inline GMPFloat operator-(const GMPFloat &lhs, double rhs)
    {
        GMPFloat result = GMPFloat::with_precision(lhs.precision());
        mpfr_sub_d(result.get_mpfr(), lhs.get_mpfr(), rhs, MPFR_RNDN);
        return result;
    }
    inline GMPFloat operator-(double lhs, const GMPFloat &rhs)
    {
        GMPFloat result = GMPFloat::with_precision(rhs.precision());
        mpfr_d_sub(result.get_mpfr(), lhs, rhs.get_mpfr(), MPFR_RNDN);
        return result;
    }
    inline GMPFloat operator*(const GMPFloat &lhs, double rhs)
    {
        GMPFloat result = GMPFloat::with_precision(lhs.precision());
        mpfr_mul_d(result.get_mpfr(), lhs.get_mpfr(), rhs, MPFR_RNDN);
        return result;
    }
    inline GMPFloat operator*(double lhs, const GMPFloat &rhs)
    {
        GMPFloat result = GMPFloat::with_precision(rhs.precision());
        mpfr_mul_d(result.get_mpfr(), rhs.get_mpfr(), lhs, MPFR_RNDN);
        return result;
    }
    inline GMPFloat operator/(const GMPFloat &lhs, double rhs)
    {
        GMPFloat result = GMPFloat::with_precision(lhs.precision());
        mpfr_div_d(result.get_mpfr(), lhs.get_mpfr(), rhs, MPFR_RNDN);
        return result;
    }
    inline GMPFloat operator/(double lhs, const GMPFloat &rhs)
    {
        GMPFloat result = GMPFloat::with_precision(rhs.precision());
        mpfr_d_div(result.get_mpfr(), lhs, rhs.get_mpfr(), MPFR_RNDN);
        return result;
    }

    // Optimized mixed operations with int using mpfr_*_si helpers
    inline GMPFloat operator+(const GMPFloat &lhs, int rhs)
    {
        GMPFloat result = GMPFloat::with_precision(lhs.precision());
        mpfr_add_si(result.get_mpfr(), lhs.get_mpfr(), static_cast<long>(rhs), MPFR_RNDN);
        return result;
    }
    inline GMPFloat operator+(int lhs, const GMPFloat &rhs)
    {
        GMPFloat result = GMPFloat::with_precision(rhs.precision());
        mpfr_add_si(result.get_mpfr(), rhs.get_mpfr(), static_cast<long>(lhs), MPFR_RNDN);
        return result;
    }
    inline GMPFloat operator-(const GMPFloat &lhs, int rhs)
    {
        GMPFloat result = GMPFloat::with_precision(lhs.precision());
        mpfr_sub_si(result.get_mpfr(), lhs.get_mpfr(), static_cast<long>(rhs), MPFR_RNDN);
        return result;
    }
    inline GMPFloat operator-(int lhs, const GMPFloat &rhs)
    {
        GMPFloat result = GMPFloat::with_precision(rhs.precision());
        mpfr_si_sub(result.get_mpfr(), static_cast<long>(lhs), rhs.get_mpfr(), MPFR_RNDN);
        return result;
    }
    inline GMPFloat operator*(const GMPFloat &lhs, int rhs)
    {
        GMPFloat result = GMPFloat::with_precision(lhs.precision());
        mpfr_mul_si(result.get_mpfr(), lhs.get_mpfr(), static_cast<long>(rhs), MPFR_RNDN);
        return result;
    }
    inline GMPFloat operator*(int lhs, const GMPFloat &rhs)
    {
        GMPFloat result = GMPFloat::with_precision(rhs.precision());
        mpfr_mul_si(result.get_mpfr(), rhs.get_mpfr(), static_cast<long>(lhs), MPFR_RNDN);
        return result;
    }
    inline GMPFloat operator/(const GMPFloat &lhs, int rhs)
    {
        GMPFloat result = GMPFloat::with_precision(lhs.precision());
        mpfr_div_si(result.get_mpfr(), lhs.get_mpfr(), static_cast<long>(rhs), MPFR_RNDN);
        return result;
    }
    inline GMPFloat operator/(int lhs, const GMPFloat &rhs)
    {
        GMPFloat result = GMPFloat::with_precision(rhs.precision());
        mpfr_si_div(result.get_mpfr(), static_cast<long>(lhs), rhs.get_mpfr(), MPFR_RNDN);
        return result;
    }

    // Comparison operators with int
    inline bool operator==(const GMPFloat &lhs, int rhs) { return lhs == GMPFloat(rhs); }
    inline bool operator!=(const GMPFloat &lhs, int rhs) { return lhs != GMPFloat(rhs); }
    inline bool operator<(const GMPFloat &lhs, int rhs) { return lhs < GMPFloat(rhs); }
    inline bool operator<=(const GMPFloat &lhs, int rhs) { return lhs <= GMPFloat(rhs); }
    inline bool operator>(const GMPFloat &lhs, int rhs) { return lhs > GMPFloat(rhs); }
    inline bool operator>=(const GMPFloat &lhs, int rhs) { return lhs >= GMPFloat(rhs); }

    inline bool operator==(int lhs, const GMPFloat &rhs) { return GMPFloat(lhs) == rhs; }
    inline bool operator!=(int lhs, const GMPFloat &rhs) { return GMPFloat(lhs) != rhs; }
    inline bool operator<(int lhs, const GMPFloat &rhs) { return GMPFloat(lhs) < rhs; }
    inline bool operator<=(int lhs, const GMPFloat &rhs) { return GMPFloat(lhs) <= rhs; }
    inline bool operator>(int lhs, const GMPFloat &rhs) { return GMPFloat(lhs) > rhs; }
    inline bool operator>=(int lhs, const GMPFloat &rhs) { return GMPFloat(lhs) >= rhs; }

    // Mixed operations between GMPInteger and GMPFloat
    inline GMPFloat operator+(const GMPInteger &lhs, const GMPFloat &rhs)
    {
        return GMPFloat(lhs) + rhs;
    }
    inline GMPFloat operator-(const GMPInteger &lhs, const GMPFloat &rhs)
    {
        return GMPFloat(lhs) - rhs;
    }
    inline GMPFloat operator*(const GMPInteger &lhs, const GMPFloat &rhs)
    {
        return GMPFloat(lhs) * rhs;
    }
    inline GMPFloat operator/(const GMPInteger &lhs, const GMPFloat &rhs)
    {
        return GMPFloat(lhs) / rhs;
    }

    inline GMPFloat operator+(const GMPFloat &lhs, const GMPInteger &rhs)
    {
        return lhs + GMPFloat(rhs);
    }
    inline GMPFloat operator-(const GMPFloat &lhs, const GMPInteger &rhs)
    {
        return lhs - GMPFloat(rhs);
    }
    inline GMPFloat operator*(const GMPFloat &lhs, const GMPInteger &rhs)
    {
        return lhs * GMPFloat(rhs);
    }
    inline GMPFloat operator/(const GMPFloat &lhs, const GMPInteger &rhs)
    {
        return lhs / GMPFloat(rhs);
    }

    // Comparison operators between GMPInteger and GMPFloat
    inline bool operator==(const GMPInteger &lhs, const GMPFloat &rhs)
    {
        return GMPFloat(lhs) == rhs;
    }
    inline bool operator!=(const GMPInteger &lhs, const GMPFloat &rhs)
    {
        return GMPFloat(lhs) != rhs;
    }
    inline bool operator<(const GMPInteger &lhs, const GMPFloat &rhs)
    {
        return GMPFloat(lhs) < rhs;
    }
    inline bool operator<=(const GMPInteger &lhs, const GMPFloat &rhs)
    {
        return GMPFloat(lhs) <= rhs;
    }
    inline bool operator>(const GMPInteger &lhs, const GMPFloat &rhs)
    {
        return GMPFloat(lhs) > rhs;
    }
    inline bool operator>=(const GMPInteger &lhs, const GMPFloat &rhs)
    {
        return GMPFloat(lhs) >= rhs;
    }

    inline bool operator==(const GMPFloat &lhs, const GMPInteger &rhs)
    {
        return lhs == GMPFloat(rhs);
    }
    inline bool operator!=(const GMPFloat &lhs, const GMPInteger &rhs)
    {
        return lhs != GMPFloat(rhs);
    }
    inline bool operator<(const GMPFloat &lhs, const GMPInteger &rhs)
    {
        return lhs < GMPFloat(rhs);
    }
    inline bool operator<=(const GMPFloat &lhs, const GMPInteger &rhs)
    {
        return lhs <= GMPFloat(rhs);
    }
    inline bool operator>(const GMPFloat &lhs, const GMPInteger &rhs)
    {
        return lhs > GMPFloat(rhs);
    }
    inline bool operator>=(const GMPFloat &lhs, const GMPInteger &rhs)
    {
        return lhs >= GMPFloat(rhs);
    }

} // namespace gridsynth
