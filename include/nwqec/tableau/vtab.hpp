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
        // Constructor with gates and parameters
        VTab(size_t n_qubits, size_t n_gate_stabs,
             const std::vector<Operation::Type> &gates = {},
             const std::vector<size_t> &qa = {},
             const std::vector<size_t> &qb = {},
             const std::vector<uint8_t> &phases = {},
             const std::vector<PauliOp> &stab_rows = {})
            : n_qubits(n_qubits)
        {
            str_len = static_cast<int>(n_qubits + 1);
            assert(gates.size() == qa.size() && gates.size() == qb.size() && phases.size() == gates.size());

            init_structure(n_qubits + n_gate_stabs);
            init_identity();
            process_gates(gates, qa, qb, phases, stab_rows);
        }

        size_t num_qubits() const { return n_qubits; }
        size_t num_rows() const { return local_rows; }

        void add_t_stab(size_t qubit, uint8_t phase)
        {
            if (local_rows % packed_size == 0)
                cur_elements++;
            size_t bit_pos = local_rows % packed_size;
            Utils::set_bit(r, cur_elements - 1, bit_pos, phase != 0);
            Utils::set_bit(z[qubit], cur_elements - 1, bit_pos, true);
            local_rows++;
        }

        void add_stab(const PauliOp &row)
        {
            if (local_rows % packed_size == 0)
                cur_elements++;
            size_t bit_pos = local_rows % packed_size;
            Utils::set_bit(r, cur_elements - 1, bit_pos, row.r());

            // Only set bits for non-zero entries (sparse representation)
            for (size_t qubit : row.x_indices())
            {
                Utils::set_bit(x[qubit], cur_elements - 1, bit_pos, true);
            }
            for (size_t qubit : row.z_indices())
            {
                Utils::set_bit(z[qubit], cur_elements - 1, bit_pos, true);
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
                        if (Utils::get_bit(x[k], i, j))
                            row.add_x(k);
                        if (Utils::get_bit(z[k], i, j))
                            row.add_z(k);
                    }
                    stabs.push_back(row);
                }
            }
            return stabs;
        }

        void apply_s_from_start(const std::vector<size_t> &qubits)
        {
            size_t start_element = start_row_index / packed_size;

            for (size_t q : qubits)
            {
                if (q < n_qubits)
                {
                    for (size_t i = start_element; i < cur_elements; i++)
                    {
                        r[i] ^= (x[q][i] & z[q][i]);
                        z[q][i] ^= x[q][i];
                    }
                }
            }
        }

        PauliOp pop_front()
        {
            if (is_empty())
                return PauliOp(n_qubits);

            size_t elem = start_row_index / packed_size;
            size_t bit = start_row_index % packed_size;

            PauliOp row(n_qubits);
            row.set_r(Utils::get_bit(r, elem, bit));

            for (size_t k = 0; k < n_qubits; k++)
            {
                if (Utils::get_bit(x[k], elem, bit))
                    row.add_x(k);
                if (Utils::get_bit(z[k], elem, bit))
                    row.add_z(k);
            }

            start_row_index++;
            return row;
        }

        bool is_empty() const
        {
            return start_row_index >= local_rows;
        }

        size_t remaining_rows() const
        {
            return is_empty() ? 0 : local_rows - start_row_index;
        }

    private:
        void init_structure(size_t total_rows)
        {
            size_t elements = Utils::calc_elements(total_rows);

            x.resize(n_qubits, std::vector<packed_t>(elements, 0));
            z.resize(n_qubits, std::vector<packed_t>(elements, 0));
            r.resize(elements, 0);
            local_rows = 0;
            cur_elements = 0;
            start_row_index = 0;
        }

        void init_identity()
        {
            for (size_t i = 0; i < n_qubits; i++)
            {
                size_t elem = local_rows / packed_size;
                size_t bit = local_rows % packed_size;
                Utils::set_bit(z[i], elem, bit, true);
                local_rows++;
            }
            cur_elements = Utils::calc_elements(local_rows);
            next_rank = 0;
        }

        void process_gates(const std::vector<Operation::Type> &gates,
                           const std::vector<size_t> &qa,
                           const std::vector<size_t> &qb,
                           const std::vector<uint8_t> &phases,
                           const std::vector<PauliOp> &stab_rows = {})
        {
            size_t stab_idx = 0;
            for (size_t i = 0; i < gates.size(); i++)
            {
                if (gates[i] == Operation::Type::T || gates[i] == Operation::Type::TDG)
                    add_t_stab(qa[i], phases[i]);
                else if (gates[i] == Operation::Type::T_PAULI || gates[i] == Operation::Type::S_PAULI)
                {
                    assert(stab_idx < stab_rows.size());

                    add_stab(stab_rows[stab_idx++]);
                }
                else
                    apply_gate(gates[i], qa[i], qb[i]);
            }
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
                std::cout << "Non-Clifford Gate!" << std::endl;
            }
        }

        void apply_h(size_t q)
        {
            for (size_t i = 0; i < cur_elements; i++)
            {
                r[i] ^= (x[q][i] & z[q][i]);
                std::swap(x[q][i], z[q][i]);
            }
        }

        void apply_s(size_t q)
        {
            for (size_t i = 0; i < cur_elements; i++)
            {
                r[i] ^= (x[q][i] & z[q][i]);
                z[q][i] ^= x[q][i];
            }
        }

        void apply_sdg(size_t q)
        {
            for (size_t i = 0; i < cur_elements; i++)
            {
                r[i] ^= x[q][i] ^ (x[q][i] & z[q][i]);
                z[q][i] ^= x[q][i];
            }
        }

        void apply_sx(size_t q)
        {
            for (size_t i = 0; i < cur_elements; i++)
            {
                r[i] ^= (x[q][i] & z[q][i]) ^ z[q][i];
                x[q][i] ^= z[q][i];
            }
        }

        void apply_sxdg(size_t q)
        {
            for (size_t i = 0; i < cur_elements; i++)
            {
                r[i] ^= (x[q][i] & z[q][i]);
                x[q][i] ^= z[q][i];
            }
        }

        void apply_cx(size_t ctrl, size_t targ)
        {
            assert(targ != SIZE_MAX);
            for (size_t i = 0; i < cur_elements; i++)
            {
                r[i] ^= ((x[ctrl][i] & z[targ][i]) & (x[targ][i] ^ z[ctrl][i] ^ MAX_PACKED));
                x[targ][i] ^= x[ctrl][i];
                z[ctrl][i] ^= z[targ][i];
            }
        }

        void apply_pauli(Operation::Type gate, size_t q)
        {
            if (gate == Operation::Type::X)
                for (size_t i = 0; i < cur_elements; i++)
                    r[i] ^= z[q][i];
            else if (gate == Operation::Type::Y)
                for (size_t i = 0; i < cur_elements; i++)
                    r[i] ^= (x[q][i] ^ z[q][i]);
            else if (gate == Operation::Type::Z)
                for (size_t i = 0; i < cur_elements; i++)
                    r[i] ^= x[q][i];
        }

        size_t n_qubits, local_rows, cur_elements, start_row_index;
        int next_rank, str_len;

        std::vector<std::vector<packed_t>> x, z;
        std::vector<packed_t> r;
    };

} // namespace NWQEC