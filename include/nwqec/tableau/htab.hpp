#pragma once

#include "nwqec/core/circuit.hpp"
#include "nwqec/core/pauli_op.hpp"

#include "vtab.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cassert>
#include <iomanip>
#include <memory>
#include <array>

#if defined(_MSC_VER) && !defined(__clang__)
#    include <intrin.h>
#endif

namespace NWQEC
{
    class HTab
    {
    public:

        HTab(size_t n_qubits) : n_qubits(n_qubits)
        {
            // rows.reserve(50);
        }

        size_t num_qubits() const { return n_qubits; }
        size_t num_rows() const
        {
            size_t count = 0;
            for (const auto &row : rows)
            {
                if (row.is_valid())
                    count++;
            }
            return count;
        }

        void add_stab(const PauliOp &pauli_op)
        {
            rows.push_back(pauli_op);
        }

        std::vector<PauliOp> get_stabs() const
        {
            std::vector<PauliOp> stabs;
            stabs.reserve(rows.size());

            for (const auto &row : rows)
            {
                if (!row.is_valid())
                    continue;
                stabs.push_back(row);
            }

            return stabs;
        }

        std::vector<std::string> get_str() const
        {
            std::vector<std::string> stabs;
            stabs.reserve(rows.size());

            for (const auto &row : rows)
            {
                if (!row.is_valid())
                    continue;
                stabs.push_back(row.to_string());
            }

            return stabs;
        }

    static inline int popcount64(uint64_t value)
    {
#if defined(_MSC_VER) && !defined(__clang__)
        return static_cast<int>(__popcnt64(static_cast<unsigned __int64>(value)));
#else
        return __builtin_popcountll(static_cast<unsigned long long>(value));
#endif
    }

        bool commutes_with_all(const PauliOp &pauli_op) const
        {
            for (const auto &tableau_row : rows)
            {
                if (!tableau_row.is_valid())
                    continue;

                if (!commutes(pauli_op, tableau_row))
                {
                    return false;
                }
            }
            return true;
        }

        void front_multiply_pauli(const PauliOp &new_pauli)
        {
            for (auto &row : rows)
            {
                if (!row.is_valid())
                    continue;

                // Compute g_function for phase calculation
                int g_val = compute_g_function(new_pauli, row);
                bool anti_commute = (g_val & 1) != 0;

                if (anti_commute)
                {
                    // Adjust g_val for anti-commuting case (matches Python line 517)
                    g_val += 1;
                }

                // Update X and Z bits (XOR operation)
                // Python line 521-523: stab1_expand = np.outer(mask, new_pauli_mtx[:2*self.qubits])
                // stab_new = stab1_expand^stab2
                // This means: if anti_commute, XOR with new_pauli; if commute, XOR with zeros (no change)
                if (row.is_small())
                {
                    if (anti_commute)
                    {
                        row.get_x_bits_small() ^= new_pauli.get_x_bits_small();
                        row.get_z_bits_small() ^= new_pauli.get_z_bits_small();
                    }
                    // If they commute, no change to X and Z bits
                }
                else
                {
                    if (anti_commute)
                    {
                        for (size_t i = 0; i < row.get_x_bits_large().size(); ++i)
                        {
                            row.get_x_bits_large()[i] ^= new_pauli.get_x_bits_large()[i];
                            row.get_z_bits_large()[i] ^= new_pauli.get_z_bits_large()[i];
                        }
                    }
                    // If they commute, no change to X and Z bits
                }

                // Update phase only for anti-commuting pairs (matches Python line 524)
                if (anti_commute)
                {
                    row.set_phase(row.get_phase() ^ new_pauli.get_phase() ^ ((g_val >> 1) & 1));
                }
            }
        }

        bool apply_reduction()
        {
            bool reduced = false;
            for (size_t i = 0; i < rows.size(); ++i)
            {
                if (!rows[i].is_valid())
                    continue;

                for (size_t j = i + 1; j < rows.size(); ++j)
                {
                    if (!rows[j].is_valid())
                        continue;

                    // Only check pairs if they are the same type
                    if (rows[i].get_rowtype() != rows[j].get_rowtype())
                        continue;

                    if (same_pauli_bits(rows[i], rows[j]))
                    {
                        reduced = true;

                        if (rows[i].get_phase() != rows[j].get_phase())
                        {
                            // Opposite phases cancel out
                            rows[i].set_valid(false);
                            rows[j].set_valid(false);
                        }
                        else
                        {
                            // Same phase merge: T merges to S, S merges to Z
                            if (rows[i].get_rowtype() == RowType::T)
                            {
                                rows[i].set_rowtype(RowType::S); // T + T -> S
                            }
                            else if (rows[i].get_rowtype() == RowType::S)
                            {
                                rows[i].set_rowtype(RowType::Z); // S + S -> Z
                            }

                            rows[i].set_valid(true);  // Keep row I valid
                            rows[j].set_valid(false); // Mark row J as invalid
                        }

                        break;
                    }
                }
            }
            return reduced;
        }

        std::vector<PauliOp> get_rows() const
        {
            std::vector<PauliOp> resulting_rows;

            for (const auto &row : rows)
            {
                if (row.is_valid())
                {
                    resulting_rows.push_back(row);
                }
            }

            return resulting_rows;
        }

    private:
        int compute_g_function(const PauliOp &pauli1, const PauliOp &pauli2) const
        {
            int g_val = 0;

            if (pauli1.is_small())
            {
                // Fast path for ≤64 qubits
                for (size_t q = 0; q < n_qubits; ++q)
                {
                    uint64_t mask = 1ULL << q;
                    bool x1 = pauli1.get_x_bits_small() & mask;
                    bool z1 = pauli1.get_z_bits_small() & mask;
                    bool x2 = pauli2.get_x_bits_small() & mask;
                    bool z2 = pauli2.get_z_bits_small() & mask;

                    g_val += (x1 & z1) * (z2 - x2) +
                             (x1 & !z1) * z2 * (2 * x2 - 1) +
                             (!x1 & z1) * x2 * (1 - 2 * z2);
                }
            }
            else
            {
                // Large circuit path
                for (size_t q = 0; q < n_qubits; ++q)
                {
                    size_t word_idx = q / 64;
                    size_t bit_idx = q % 64;
                    uint64_t mask = 1ULL << bit_idx;

                    bool x1 = pauli1.get_x_bits_large()[word_idx] & mask;
                    bool z1 = pauli1.get_z_bits_large()[word_idx] & mask;
                    bool x2 = pauli2.get_x_bits_large()[word_idx] & mask;
                    bool z2 = pauli2.get_z_bits_large()[word_idx] & mask;

                    g_val += (x1 & z1) * (z2 - x2) +
                             (x1 & !z1) * z2 * (2 * x2 - 1) +
                             (!x1 & z1) * x2 * (1 - 2 * z2);
                }
            }

            return g_val & 3;
        }

        bool same_pauli_bits(const PauliOp &row1, const PauliOp &row2) const
        {
            if (row1.is_small())
            {
                return row1.get_x_bits_small() == row2.get_x_bits_small() &&
                       row1.get_z_bits_small() == row2.get_z_bits_small();
            }
            else
            {
                return row1.get_x_bits_large() == row2.get_x_bits_large() &&
                       row1.get_z_bits_large() == row2.get_z_bits_large();
            }
        }

        bool commutes(const PauliOp &row1, const PauliOp &row2) const
        {
            if (row1.is_small())
            {
                // Ultra-fast path for ≤64 qubits - single popcount operation!
                uint64_t anti_commute_word = (row1.get_x_bits_small() & row2.get_z_bits_small()) ^
                                             (row1.get_z_bits_small() & row2.get_x_bits_small());
                return (popcount64(anti_commute_word) & 1) == 0;
            }
            else
            {
                // Large circuit path with loop unrolling
                size_t anti_commuting_pairs = 0;
                size_t num_words = row1.get_x_bits_large().size();
                size_t i = 0;

                // Process 4 words at a time when possible
                for (; i + 3 < num_words; i += 4)
                {
                    // Unroll 4 iterations for better instruction-level parallelism
                    uint64_t ac0 = (row1.get_x_bits_large()[i] & row2.get_z_bits_large()[i]) ^ (row1.get_z_bits_large()[i] & row2.get_x_bits_large()[i]);
                    uint64_t ac1 = (row1.get_x_bits_large()[i + 1] & row2.get_z_bits_large()[i + 1]) ^ (row1.get_z_bits_large()[i + 1] & row2.get_x_bits_large()[i + 1]);
                    uint64_t ac2 = (row1.get_x_bits_large()[i + 2] & row2.get_z_bits_large()[i + 2]) ^ (row1.get_z_bits_large()[i + 2] & row2.get_x_bits_large()[i + 2]);
                    uint64_t ac3 = (row1.get_x_bits_large()[i + 3] & row2.get_z_bits_large()[i + 3]) ^ (row1.get_z_bits_large()[i + 3] & row2.get_x_bits_large()[i + 3]);

                    anti_commuting_pairs += popcount64(ac0) + popcount64(ac1) +
                                            popcount64(ac2) + popcount64(ac3);
                }

                // Handle remaining words
                for (; i < num_words; ++i)
                {
                    uint64_t anti_commute_word = (row1.get_x_bits_large()[i] & row2.get_z_bits_large()[i]) ^
                                                 (row1.get_z_bits_large()[i] & row2.get_x_bits_large()[i]);
                    anti_commuting_pairs += popcount64(anti_commute_word);
                }

                return (anti_commuting_pairs & 1) == 0;
            }
        }

        size_t n_qubits;
        std::vector<PauliOp> rows;
    };

} // namespace NWQEC