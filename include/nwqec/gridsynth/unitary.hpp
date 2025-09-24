#pragma once

#include <vector>
#include <string>
#include <array>
#include <iostream>
#include <complex>

#include "ring.hpp"

namespace gridsynth
{

    /**
     * DOmegaUnitary: Represents a 2x2 unitary matrix with entries in D[ω]
     * Used for quantum circuit synthesis with Clifford+T gates
     */
    class DOmegaUnitary
    {
    private:
        DOmega _z, _w;
        int _n;

    public:
        // Constructor
        DOmegaUnitary(const DOmega &z, const DOmega &w, int n, int k = -1)
            : _z(z), _w(w), _n(n & 0b111)
        {

            if (k == -1)
            {
                // Auto-align denominators
                if (_z.k() > _w.k())
                {
                    _w = _w.renew_denomexp(_z.k());
                }
                else if (_z.k() < _w.k())
                {
                    _z = _z.renew_denomexp(_w.k());
                }
            }
            else
            {
                // Use specified denominator exponent
                _z = _z.renew_denomexp(k);
                _w = _w.renew_denomexp(k);
            }
        }

        // Properties
        DOmega z() const { return _z; }
        DOmega w() const { return _w; }
        int n() const { return _n; }
        int k() const { return _w.k(); }

        // Matrix representation as 2x2 array of DOmega
        std::array<std::array<DOmega, 2>, 2> to_matrix() const
        {
            DOmega m00 = _z;
            DOmega m01 = -_w.conj().mul_by_omega_power(_n);
            DOmega m10 = _w;
            DOmega m11 = _z.conj().mul_by_omega_power(_n);

            return {{{{m00, m01}},
                     {{m10, m11}}}};
        }

        // Complex matrix representation (avoid intermediate arrays and extra temporaries)
        std::array<std::array<std::complex<Float>, 2>, 2> to_complex_matrix() const
        {
            // Build entries directly and use fast coord extraction with hoisted invariants.
            DOmega m00 = _z;
            DOmega m01 = -_w.conj().mul_by_omega_power(_n);
            DOmega m10 = _w;
            DOmega m11 = _z.conj().mul_by_omega_power(_n);

            // All entries share the same denominator exponent k
            const Float inv_scale = Float(1.0) / pow_sqrt2(_w.k());
            const Float sqrt2_over_2 = SQRT2 / Float(2.0);

            Float r00, i00, r01, i01, r10, i10, r11, i11;
            m00.coords_into_with(inv_scale, sqrt2_over_2, r00, i00);
            m01.coords_into_with(inv_scale, sqrt2_over_2, r01, i01);
            m10.coords_into_with(inv_scale, sqrt2_over_2, r10, i10);
            m11.coords_into_with(inv_scale, sqrt2_over_2, r11, i11);

            return {{{{std::complex<Float>(r00, i00), std::complex<Float>(r01, i01)}}
                    ,{{std::complex<Float>(r10, i10), std::complex<Float>(r11, i11)}}}};
        }

        // String representation
        friend std::ostream &operator<<(std::ostream &os, const DOmegaUnitary &u)
        {
            auto matrix = u.to_matrix();
            os << "[[" << matrix[0][0] << ", " << matrix[0][1] << "],\n";
            os << " [" << matrix[1][0] << ", " << matrix[1][1] << "]]";
            return os;
        }

        // Equality operator
        bool operator==(const DOmegaUnitary &other) const
        {
            return _z == other._z && _w == other._w && _n == other._n;
        }

        bool operator!=(const DOmegaUnitary &other) const
        {
            return !(*this == other);
        }

        // Gate multiplications from the left
        DOmegaUnitary mul_by_T_from_left() const
        {
            return DOmegaUnitary(_z, _w.mul_by_omega(), _n + 1);
        }

        DOmegaUnitary mul_by_T_inv_from_left() const
        {
            return DOmegaUnitary(_z, _w.mul_by_omega_inv(), _n - 1);
        }

        DOmegaUnitary mul_by_T_power_from_left(int m) const
        {
            m &= 0b111; // mod 8
            return DOmegaUnitary(_z, _w.mul_by_omega_power(m), _n + m);
        }

        DOmegaUnitary mul_by_S_from_left() const
        {
            return DOmegaUnitary(_z, _w.mul_by_omega_power(2), _n + 2);
        }

        DOmegaUnitary mul_by_S_power_from_left(int m) const
        {
            m &= 0b11; // mod 4
            return DOmegaUnitary(_z, _w.mul_by_omega_power(m << 1), _n + (m << 1));
        }

        DOmegaUnitary mul_by_H_from_left() const
        {
            DOmega new_z = (_z + _w).mul_by_inv_sqrt2();
            DOmega new_w = (_z - _w).mul_by_inv_sqrt2();
            return DOmegaUnitary(new_z, new_w, _n + 4);
        }

        DOmegaUnitary mul_by_H_and_T_power_from_left(int m) const
        {
            return mul_by_T_power_from_left(m).mul_by_H_from_left();
        }

        DOmegaUnitary mul_by_X_from_left() const
        {
            return DOmegaUnitary(_w, _z, _n + 4);
        }

        DOmegaUnitary mul_by_W_from_left() const
        {
            return DOmegaUnitary(_z.mul_by_omega(), _w.mul_by_omega(), _n + 2);
        }

        DOmegaUnitary mul_by_W_power_from_left(int m) const
        {
            m &= 0b111; // mod 8
            return DOmegaUnitary(_z.mul_by_omega_power(m), _w.mul_by_omega_power(m), _n + (m << 1));
        }

        // Denominator operations
        DOmegaUnitary renew_denomexp(int new_k) const
        {
            return DOmegaUnitary(_z, _w, _n, new_k);
        }

        DOmegaUnitary reduce_denomexp() const
        {
            DOmega new_z = _z.reduce_denomexp();
            DOmega new_w = _w.reduce_denomexp();
            return DOmegaUnitary(new_z, new_w, _n);
        }

        // Static factory methods
        static DOmegaUnitary identity()
        {
            return DOmegaUnitary(DOmega::from_int(1), DOmega::from_int(0), 0);
        }

        static DOmegaUnitary from_gates(const std::string &gates)
        {
            DOmegaUnitary unitary = identity();

            // Process gates in reverse order (right to left multiplication)
            for (auto it = gates.rbegin(); it != gates.rend(); ++it)
            {
                const char &gate = *it;
                if (gate == 'H')
                {
                    unitary = unitary.renew_denomexp(unitary.k() + 1).mul_by_H_from_left();
                }
                else if (gate == 'T')
                {
                    unitary = unitary.mul_by_T_from_left();
                }
                else if (gate == 'S')
                {
                    unitary = unitary.mul_by_S_from_left();
                }
                else if (gate == 'X')
                {
                    unitary = unitary.mul_by_X_from_left();
                }
                else if (gate == 'W')
                {
                    unitary = unitary.mul_by_W_from_left();
                }
                // Ignore unknown gates silently
            }

            return unitary.reduce_denomexp();
        }
    };

} // namespace gridsynth
