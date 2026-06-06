#pragma once

#include "nwqec/core/circuit.hpp"
#include "nwqec/core/pauli_op.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cassert>
#include <iomanip>
#include <memory>
#include <array>

namespace NWQEC
{
#if defined(_MSC_VER) && !defined(__clang__)
#    define NWQEC_RESTRICT __restrict
#else
#    define NWQEC_RESTRICT __restrict__
#endif

    using packed_t = uint64_t;
    constexpr size_t packed_size = sizeof(packed_t) * 8;
    constexpr packed_t MAX_PACKED = static_cast<packed_t>(~0);

    namespace Utils
    {
        inline size_t calc_elements(size_t rows)
        {
            size_t n = rows / packed_size + 1;
            if (rows % packed_size == 0)
                n--;
            return n;
        }

        inline void set_bit(std::vector<packed_t> &vec, size_t elem, size_t bit, bool val)
        {
            if (val)
                vec[elem] |= (static_cast<packed_t>(1) << bit);
            else
                vec[elem] &= ~(static_cast<packed_t>(1) << bit);
        }

        inline bool get_bit(const std::vector<packed_t> &vec, size_t elem, size_t bit)
        {
            return (vec[elem] & (static_cast<packed_t>(1) << bit)) != 0;
        }
    }

    class VTab
    {
    public:
        VTab(size_t n_qubits, size_t n_gate_stabs)
            : n_qubits(n_qubits)
        {
            init_structure(n_qubits + n_gate_stabs);
            init_identity();
        }

        size_t num_qubits() const { return n_qubits; }
        size_t num_rows() const { return local_rows; }

        void apply_clifford_gate(Operation::Type gate, size_t a, size_t b = SIZE_MAX)
        {
            apply_gate(gate, a, b);
        }

        void add_t_stab(size_t qubit, uint8_t phase)
        {
            if (local_rows % packed_size == 0)
                cur_elements++;
            size_t bit_pos = local_rows % packed_size;
            packed_t bit_mask = static_cast<packed_t>(1) << bit_pos;
            if (phase != 0)
                r[cur_elements - 1] |= bit_mask;
            z_data[qubit * elements_per_qubit + cur_elements - 1] |= bit_mask;
            local_rows++;
        }

        void add_stab(const PauliOp &row)
        {
            if (row.is_small())
            {
                add_stab_bits(row.get_x_bits_small(), row.get_z_bits_small(), row.r());
                return;
            }

            if (local_rows % packed_size == 0)
                cur_elements++;
            size_t bit_pos = local_rows % packed_size;
            packed_t row_mask = static_cast<packed_t>(1) << bit_pos;

            if (row.r())
                r[cur_elements - 1] |= row_mask;

            for (size_t qubit : row.x_indices())
                x_data[qubit * elements_per_qubit + cur_elements - 1] |= row_mask;
            for (size_t qubit : row.z_indices())
                z_data[qubit * elements_per_qubit + cur_elements - 1] |= row_mask;

            local_rows++;
        }

        void add_ccx_stabs(size_t q0, size_t q1, size_t q2)
        {
            if (n_qubits > 64)
            {
                auto ccx_rows = PauliOp::create_ccx_ops(q0, q1, q2, n_qubits);
                for (const auto &stab : ccx_rows)
                    add_stab(stab);
                return;
            }

            const packed_t z0 = static_cast<packed_t>(1) << q0;
            const packed_t z1 = static_cast<packed_t>(1) << q1;
            const packed_t x2 = static_cast<packed_t>(1) << q2;

            add_stab_bits(x2, 0, false);
            add_stab_bits(0, z0 | z1, true);
            add_stab_bits(0, z0, false);
            add_stab_bits(x2, z0, true);
            add_stab_bits(0, z1, false);
            add_stab_bits(x2, z0 | z1, false);
            add_stab_bits(x2, z1, true);
        }

        void add_stab_bits(packed_t x_bits, packed_t z_bits, bool phase)
        {
            if (local_rows % packed_size == 0)
                cur_elements++;
            size_t bit_pos = local_rows % packed_size;
            packed_t row_mask = static_cast<packed_t>(1) << bit_pos;

            if (phase)
                r[cur_elements - 1] |= row_mask;

            while (x_bits)
            {
                int bit = pauli_op_ctz64(x_bits);
                x_data[static_cast<size_t>(bit) * elements_per_qubit + cur_elements - 1] |= row_mask;
                x_bits &= x_bits - 1;
            }
            while (z_bits)
            {
                int bit = pauli_op_ctz64(z_bits);
                z_data[static_cast<size_t>(bit) * elements_per_qubit + cur_elements - 1] |= row_mask;
                z_bits &= z_bits - 1;
            }
            local_rows++;
        }

        std::vector<PauliOp> get_paili_ops()
        {
            std::vector<PauliOp> stabs;
            stabs.reserve(local_rows);

            for (size_t i = 0; i < cur_elements; i++)
            {
                for (size_t j = 0; j < packed_size; j++)
                {
                    if (stabs.size() >= local_rows)
                        break;
                    PauliOp row(n_qubits);
                    row.set_r(Utils::get_bit(r, i, j));

                    // Build sparse representation by checking each qubit
                    for (size_t k = 0; k < n_qubits; k++)
                    {
                        if ((x_data[k * elements_per_qubit + i] & (static_cast<packed_t>(1) << j)) != 0)
                            row.add_x(k);
                        if ((z_data[k * elements_per_qubit + i] & (static_cast<packed_t>(1) << j)) != 0)
                            row.add_z(k);
                    }
                    stabs.push_back(row);
                }
            }
            return stabs;
        }

    private:
        void init_structure(size_t total_rows)
        {
            size_t elements = Utils::calc_elements(total_rows);

            elements_per_qubit = elements;
            x_data.assign(n_qubits * elements_per_qubit, 0);
            z_data.assign(n_qubits * elements_per_qubit, 0);
            r.resize(elements, 0);
            local_rows = 0;
            cur_elements = 0;
        }

        void init_identity()
        {
            for (size_t i = 0; i < n_qubits; i++)
            {
                size_t elem = local_rows / packed_size;
                size_t bit = local_rows % packed_size;
                z_data[i * elements_per_qubit + elem] |= (static_cast<packed_t>(1) << bit);
                local_rows++;
            }
            cur_elements = Utils::calc_elements(local_rows);
        }

        void apply_gate(Operation::Type gate, size_t a, size_t b = SIZE_MAX)
        {
            switch (gate)
            {
            case Operation::Type::H:
                apply_h(a);
                break;
            case Operation::Type::S:
                apply_s(a);
                break;
            case Operation::Type::SDG:
                apply_sdg(a);
                break;
            case Operation::Type::SX:
                apply_sx(a);
                break;
            case Operation::Type::SXDG:
                apply_sxdg(a);
                break;
            case Operation::Type::CX:
                apply_cx(a, b);
                break;
            case Operation::Type::X:
            case Operation::Type::Y:
            case Operation::Type::Z:
                apply_pauli(gate, a);
                break;
            default:
                throw std::runtime_error("Non-Clifford Gate encountered: " + Operation::get_type_name(gate));
            }
        }

        void apply_h(size_t q)
        {
            packed_t *NWQEC_RESTRICT x_ptr = x_row(q);
            packed_t *NWQEC_RESTRICT z_ptr = z_row(q);
            packed_t *NWQEC_RESTRICT r_ptr = r.data();
            for (size_t i = 0; i < cur_elements; i++)
            {
                packed_t x_word = x_ptr[i];
                packed_t z_word = z_ptr[i];
                r_ptr[i] ^= (x_word & z_word);
                x_ptr[i] = z_word;
                z_ptr[i] = x_word;
            }
        }

        void apply_s(size_t q)
        {
            const packed_t *NWQEC_RESTRICT x_ptr = x_row(q);
            packed_t *NWQEC_RESTRICT z_ptr = z_row(q);
            packed_t *NWQEC_RESTRICT r_ptr = r.data();
            for (size_t i = 0; i < cur_elements; i++)
            {
                packed_t x_word = x_ptr[i];
                r_ptr[i] ^= (x_word & z_ptr[i]);
                z_ptr[i] ^= x_word;
            }
        }

        void apply_sdg(size_t q)
        {
            const packed_t *NWQEC_RESTRICT x_ptr = x_row(q);
            packed_t *NWQEC_RESTRICT z_ptr = z_row(q);
            packed_t *NWQEC_RESTRICT r_ptr = r.data();
            for (size_t i = 0; i < cur_elements; i++)
            {
                packed_t x_word = x_ptr[i];
                r_ptr[i] ^= x_word ^ (x_word & z_ptr[i]);
                z_ptr[i] ^= x_word;
            }
        }

        void apply_sx(size_t q)
        {
            packed_t *NWQEC_RESTRICT x_ptr = x_row(q);
            const packed_t *NWQEC_RESTRICT z_ptr = z_row(q);
            packed_t *NWQEC_RESTRICT r_ptr = r.data();
            for (size_t i = 0; i < cur_elements; i++)
            {
                packed_t z_word = z_ptr[i];
                r_ptr[i] ^= (x_ptr[i] & z_word) ^ z_word;
                x_ptr[i] ^= z_word;
            }
        }

        void apply_sxdg(size_t q)
        {
            packed_t *NWQEC_RESTRICT x_ptr = x_row(q);
            const packed_t *NWQEC_RESTRICT z_ptr = z_row(q);
            packed_t *NWQEC_RESTRICT r_ptr = r.data();
            for (size_t i = 0; i < cur_elements; i++)
            {
                packed_t z_word = z_ptr[i];
                r_ptr[i] ^= (x_ptr[i] & z_word);
                x_ptr[i] ^= z_word;
            }
        }

        void apply_cx(size_t ctrl, size_t targ)
        {
            assert(targ != SIZE_MAX);
            const packed_t *NWQEC_RESTRICT x_ctrl_ptr = x_row(ctrl);
            packed_t *NWQEC_RESTRICT z_ctrl_ptr = z_row(ctrl);
            packed_t *NWQEC_RESTRICT x_targ_ptr = x_row(targ);
            const packed_t *NWQEC_RESTRICT z_targ_ptr = z_row(targ);
            packed_t *NWQEC_RESTRICT r_ptr = r.data();
            for (size_t i = 0; i < cur_elements; i++)
            {
                packed_t x_ctrl_word = x_ctrl_ptr[i];
                packed_t z_targ_word = z_targ_ptr[i];
                packed_t x_targ_word = x_targ_ptr[i];
                packed_t z_ctrl_word = z_ctrl_ptr[i];
                r_ptr[i] ^= ((x_ctrl_word & z_targ_word) & ~(x_targ_word ^ z_ctrl_word));
                x_targ_ptr[i] = x_targ_word ^ x_ctrl_word;
                z_ctrl_ptr[i] = z_ctrl_word ^ z_targ_word;
            }
        }

        void apply_pauli(Operation::Type gate, size_t q)
        {
            const packed_t *NWQEC_RESTRICT xq = x_row(q);
            const packed_t *NWQEC_RESTRICT zq = z_row(q);
            if (gate == Operation::Type::X)
                for (size_t i = 0; i < cur_elements; i++)
                    r[i] ^= zq[i];
            else if (gate == Operation::Type::Y)
                for (size_t i = 0; i < cur_elements; i++)
                    r[i] ^= (xq[i] ^ zq[i]);
            else if (gate == Operation::Type::Z)
                for (size_t i = 0; i < cur_elements; i++)
                    r[i] ^= xq[i];
        }

        packed_t *x_row(size_t q)
        {
            return x_data.data() + q * elements_per_qubit;
        }

        packed_t *z_row(size_t q)
        {
            return z_data.data() + q * elements_per_qubit;
        }

        const packed_t *x_row(size_t q) const
        {
            return x_data.data() + q * elements_per_qubit;
        }

        const packed_t *z_row(size_t q) const
        {
            return z_data.data() + q * elements_per_qubit;
        }

        size_t n_qubits, local_rows, cur_elements;
        size_t elements_per_qubit = 0;

        std::vector<packed_t> x_data, z_data;
        std::vector<packed_t> r;
    };

} // namespace NWQEC
