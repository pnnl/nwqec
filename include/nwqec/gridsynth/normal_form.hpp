#pragma once

#include <vector>
#include <string>
#include <array>
#include <stdexcept>
#include <iostream>

#include "types.hpp"
namespace gridsynth
{

    // Enums
    enum class Axis
    {
        I = 0,
        H = 1,
        SH = 2
    };

    enum class Syllable
    {
        I = 0,
        T = 1,
        HT = 2,
        SHT = 3
    };

    // Lookup tables
    inline const std::array<std::pair<int, int>, 8> CONJ2_TABLE = {{{0, 0}, {0, 0}, {1, 0}, {3, 2}, {2, 0}, {2, 4}, {3, 0}, {1, 6}}};

    inline const std::array<std::array<int, 4>, 24> CONJ3_TABLE = {{{0, 0, 0, 0}, {0, 0, 1, 0}, {0, 0, 2, 0}, {0, 0, 3, 0}, {0, 1, 0, 0}, {0, 1, 1, 0}, {0, 1, 2, 0}, {0, 1, 3, 0}, {1, 0, 0, 0}, {2, 0, 3, 6}, {1, 1, 2, 2}, {2, 1, 3, 6}, {1, 0, 2, 0}, {2, 1, 1, 0}, {1, 1, 0, 6}, {2, 0, 1, 4}, {2, 0, 0, 0}, {1, 1, 3, 4}, {2, 1, 0, 0}, {1, 0, 1, 2}, {2, 1, 2, 2}, {1, 1, 1, 0}, {2, 0, 2, 6}, {1, 0, 3, 2}}};

    inline const std::array<std::array<int, 4>, 24> CINV_TABLE = {{{0, 0, 0, 0}, {0, 0, 3, 0}, {0, 0, 2, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 1, 1, 6}, {0, 1, 2, 4}, {0, 1, 3, 2}, {2, 0, 0, 0}, {1, 0, 1, 2}, {2, 1, 0, 0}, {1, 1, 3, 4}, {2, 1, 1, 2}, {1, 1, 1, 6}, {2, 0, 2, 2}, {1, 0, 3, 4}, {1, 0, 0, 0}, {2, 1, 3, 6}, {1, 1, 2, 2}, {2, 0, 3, 6}, {1, 0, 2, 0}, {2, 1, 1, 6}, {1, 1, 0, 2}, {2, 0, 1, 6}}};

    inline const std::array<std::array<int, 3>, 6> TCONJ_TABLE = {{{static_cast<int>(Axis::I), 0, 0}, {static_cast<int>(Axis::I), 1, 7}, {static_cast<int>(Axis::H), 3, 3}, {static_cast<int>(Axis::H), 2, 0}, {static_cast<int>(Axis::SH), 0, 5}, {static_cast<int>(Axis::SH), 1, 4}}};

    /**
     * Clifford group element representation
     */
    class Clifford
    {
    private:
        int _a, _b, _c, _d;

        // Normalize parameters to valid ranges
        void normalize()
        {
            _a = (_a % 3 + 3) % 3; // 0 <= a < 3
            _b = _b & 1;           // 0 <= b < 2
            _c = _c & 0b11;        // 0 <= c < 4
            _d = _d & 0b111;       // 0 <= d < 8
        }

    public:
        Clifford(int a = 0, int b = 0, int c = 0, int d = 0)
            : _a(a), _b(b), _c(c), _d(d)
        {
            normalize();
        }

        // Properties
        int a() const { return _a; }
        int b() const { return _b; }
        int c() const { return _c; }
        int d() const { return _d; }

        // String representation
        friend std::ostream &operator<<(std::ostream &os, const Clifford &cliff)
        {
            os << "E^" << cliff._a << " X^" << cliff._b << " S^" << cliff._c << " ω^" << cliff._d;
            return os;
        }

        // Factory method
        static Clifford from_str(const std::string &g)
        {
            if (g == "H")
            {
                return Clifford(1, 0, 1, 5); // CLIFFORD_H
            }
            else if (g == "S")
            {
                return Clifford(0, 0, 1, 0); // CLIFFORD_S
            }
            else if (g == "X")
            {
                return Clifford(0, 1, 0, 0); // CLIFFORD_X
            }
            else if (g == "W")
            {
                return Clifford(0, 0, 0, 1); // CLIFFORD_W
            }
            else
            {
                throw std::invalid_argument("Unknown gate: " + g);
            }
        }

        // Equality operator
        bool operator==(const Clifford &other) const
        {
            return _a == other._a && _b == other._b && _c == other._c && _d == other._d;
        }

        bool operator!=(const Clifford &other) const
        {
            return !(*this == other);
        }

        // Helper methods for lookup tables
        static std::pair<int, int> _conj2(int c, int b)
        {
            int index = (c << 1) | b;
            return CONJ2_TABLE[index];
        }

        static std::array<int, 4> _conj3(int b, int c, int a)
        {
            int index = (a << 3) | (b << 2) | c;
            return CONJ3_TABLE[index];
        }

        static std::array<int, 4> _cinv(int a, int b, int c)
        {
            int index = (a << 3) | (b << 2) | c;
            return CINV_TABLE[index];
        }

        static std::array<int, 3> _tconj(int a, int b)
        {
            int index = (a << 1) | b;
            return TCONJ_TABLE[index];
        }

        // Multiplication operator
        Clifford operator*(const Clifford &other) const
        {
            auto conj3_result = _conj3(_b, _c, other._a);
            int a1 = conj3_result[0];
            int b1 = conj3_result[1];
            int c1 = conj3_result[2];
            int d1 = conj3_result[3];

            auto conj2_result = _conj2(c1, other._b);
            int c2 = conj2_result.first;
            int d2 = conj2_result.second;

            int new_a = _a + a1;
            int new_b = b1 + other._b;
            int new_c = c2 + other._c;
            int new_d = d2 + d1 + _d + other._d;

            return Clifford(new_a, new_b, new_c, new_d);
        }

        // Inverse
        Clifford inv() const
        {
            auto cinv_result = _cinv(_a, _b, _c);
            int a1 = cinv_result[0];
            int b1 = cinv_result[1];
            int c1 = cinv_result[2];
            int d1 = cinv_result[3];

            return Clifford(a1, b1, c1, d1 - _d);
        }

        // Decompose coset
        std::pair<Axis, Clifford> decompose_coset() const
        {
            if (_a == 0)
            {
                return {Axis::I, *this};
            }
            else if (_a == 1)
            {
                Clifford clifford_h(1, 0, 1, 5);
                return {Axis::H, clifford_h.inv() * (*this)};
            }
            else if (_a == 2)
            {
                Clifford clifford_s(0, 0, 1, 0);
                Clifford clifford_h(1, 0, 1, 5);
                Clifford clifford_sh = clifford_s * clifford_h;
                return {Axis::SH, clifford_sh.inv() * (*this)};
            }
            // Should never reach here
            return {Axis::I, *this};
        }

        // Decompose T conjugation
        std::pair<Axis, Clifford> decompose_tconj() const
        {
            auto tconj_result = _tconj(_a, _b);
            Axis axis = static_cast<Axis>(tconj_result[0]);
            int c1 = tconj_result[1];
            int d1 = tconj_result[2];

            return {axis, Clifford(0, _b, c1 + _c, d1 + _d)};
        }

        // Convert to gate string
        std::string to_gates() const
        {
            auto [axis, c] = decompose_coset();

            std::string result = "";
            if (axis == Axis::H)
            {
                result += "H";
            }
            else if (axis == Axis::SH)
            {
                result += "SH";
            }

            result += std::string(c.b(), 'X');
            result += std::string(c.c(), 'S');
            result += std::string(c.d(), 'W');

            return result;
        }
    };

    /**
     * Normal form representation for gate sequences
     */
    class NormalForm
    {
    private:
        std::vector<Syllable> _syllables;
        Clifford _c;

    public:
        NormalForm(const std::vector<Syllable> &syllables = {}, const Clifford &c = Clifford())
            : _syllables(syllables), _c(c) {}

        // Properties
        const std::vector<Syllable> &syllables() const { return _syllables; }
        const Clifford &c() const { return _c; }
        void set_c(const Clifford &c) { _c = c; }

        // String representation
        friend std::ostream &operator<<(std::ostream &os, const NormalForm &nf)
        {
            os << "NormalForm([";
            for (size_t i = 0; i < nf._syllables.size(); ++i)
            {
                if (i > 0)
                    os << ", ";
                os << static_cast<int>(nf._syllables[i]);
            }
            os << "], " << nf._c << ")";
            return os;
        }

        // Append a gate to the normal form
        void _append_gate(const std::string &g)
        {
            if (g == "H" || g == "S" || g == "X" || g == "W")
            {
                _c = _c * Clifford::from_str(g);
            }
            else if (g == "T")
            {
                auto [axis, new_c] = _c.decompose_tconj();

                if (axis == Axis::I)
                {
                    if (_syllables.empty())
                    {
                        _syllables.push_back(Syllable::T);
                    }
                    else if (_syllables.back() == Syllable::T)
                    {
                        _syllables.pop_back();
                        Clifford clifford_s(0, 0, 1, 0);
                        _c = clifford_s * new_c;
                    }
                    else if (_syllables.back() == Syllable::HT)
                    {
                        _syllables.pop_back();
                        Clifford clifford_h(1, 0, 1, 5);
                        Clifford clifford_s(0, 0, 1, 0);
                        Clifford clifford_hs = clifford_h * clifford_s;
                        _c = clifford_hs * new_c;
                    }
                    else if (_syllables.back() == Syllable::SHT)
                    {
                        _syllables.pop_back();
                        Clifford clifford_s(0, 0, 1, 0);
                        Clifford clifford_h(1, 0, 1, 5);
                        Clifford clifford_shs = clifford_s * clifford_h * clifford_s;
                        _c = clifford_shs * new_c;
                    }
                }
                else if (axis == Axis::H)
                {
                    _syllables.push_back(Syllable::HT);
                    _c = new_c;
                }
                else if (axis == Axis::SH)
                {
                    _syllables.push_back(Syllable::SHT);
                    _c = new_c;
                }
            }
            // Ignore unknown gates silently
        }

        // Factory method
        static NormalForm from_gates(const std::string &gates)
        {
            NormalForm normal_form;
            for (char g : gates)
            {
                normal_form._append_gate(std::string(1, g));
            }
            return normal_form;
        }

        // Convert to gate string
        std::string to_gates() const
        {
            std::string gates = "";

            for (const auto &syllable : _syllables)
            {
                if (syllable == Syllable::T)
                {
                    gates += "T";
                }
                else if (syllable == Syllable::HT)
                {
                    gates += "HT";
                }
                else if (syllable == Syllable::SHT)
                {
                    gates += "SHT";
                }
                // Syllable::I adds nothing
            }

            gates += _c.to_gates();

            return gates.empty() ? "I" : gates;
        }
    };

    // Predefined Clifford elements
    inline const Clifford CLIFFORD_I(0, 0, 0, 0);
    inline const Clifford CLIFFORD_X(0, 1, 0, 0);
    inline const Clifford CLIFFORD_H(1, 0, 1, 5);
    inline const Clifford CLIFFORD_S(0, 0, 1, 0);
    inline const Clifford CLIFFORD_W(0, 0, 0, 1);
    inline const Clifford CLIFFORD_SH = CLIFFORD_S * CLIFFORD_H;
    inline const Clifford CLIFFORD_HS = CLIFFORD_H * CLIFFORD_S;
    inline const Clifford CLIFFORD_SHS = CLIFFORD_S * CLIFFORD_H * CLIFFORD_S;

} // namespace gridsynth
