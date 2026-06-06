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
        }

        size_t num_qubits() const { return n_qubits; }
        size_t num_rows() const
        {
            if (row_index_enabled)
                return indexed_row_count;

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
            if (row_index_enabled)
            {
                append_indexed_row(pauli_op);
                return;
            }

            rows.push_back(pauli_op);
        }

        std::vector<PauliOp> get_stabs() const
        {
            if (row_index_enabled)
                return get_rows();

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
            if (row_index_enabled)
            {
                std::vector<std::string> stabs;
                auto indexed_rows = get_rows();
                stabs.reserve(indexed_rows.size());
                for (const auto &row : indexed_rows)
                    stabs.push_back(row.to_string());
                return stabs;
            }

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

        static inline int ctz64(uint64_t value)
        {
#if defined(_MSC_VER) && !defined(__clang__)
            unsigned long index = 0;
            _BitScanForward64(&index, static_cast<unsigned __int64>(value));
            return static_cast<int>(index);
#else
            return __builtin_ctzll(static_cast<unsigned long long>(value));
#endif
        }

        bool commutes_with_all(const PauliOp &pauli_op) const
        {
            if (row_index_enabled && pauli_op.is_small())
            {
                const size_t words = indexed_word_count();
                anti_row_scratch.assign(words, 0);
                for_each_set_bit(pauli_op.get_x_bits_small(), [&](size_t q) {
                    if (q >= n_qubits)
                        return;
                    for (size_t word = 0; word < words; ++word)
                        anti_row_scratch[word] ^= z_row_index[q][word];
                });
                for_each_set_bit(pauli_op.get_z_bits_small(), [&](size_t q) {
                    if (q >= n_qubits)
                        return;
                    for (size_t word = 0; word < words; ++word)
                        anti_row_scratch[word] ^= x_row_index[q][word];
                });
                for (uint64_t word : anti_row_scratch)
                {
                    if (word != 0)
                        return false;
                }
                return true;
            }

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
            if (new_pauli.is_small())
            {
                if (!row_index_enabled)
                    enable_row_index();
                front_multiply_pauli_indexed(new_pauli);
                return;
            }

            for (auto &row : rows)
            {
                if (!row.is_valid())
                    continue;

                int g_val = compute_g_function(new_pauli, row);
                bool anti_commute = (g_val & 1) != 0;
                if (!anti_commute)
                    continue;

                // Adjust g_val for anti-commuting case (matches Python line 517)
                g_val += 1;

                // Update X and Z bits (XOR operation)
                // Python line 521-523: stab1_expand = np.outer(mask, new_pauli_mtx[:2*self.qubits])
                // stab_new = stab1_expand^stab2
                // This means: if anti_commute, XOR with new_pauli; if commute, XOR with zeros (no change)
                for (size_t i = 0; i < row.get_x_bits_large().size(); ++i)
                {
                    row.get_x_bits_large()[i] ^= new_pauli.get_x_bits_large()[i];
                    row.get_z_bits_large()[i] ^= new_pauli.get_z_bits_large()[i];
                }

                // Update phase only for anti-commuting pairs (matches Python line 524)
                row.set_phase(row.get_phase() ^ new_pauli.get_phase() ^ ((g_val >> 1) & 1));
            }
        }

        bool apply_reduction()
        {
            materialize_indexed_rows();

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

            if (row_index_enabled)
            {
                resulting_rows.reserve(indexed_row_count);
                for (size_t row_index = 0; row_index < indexed_row_count; ++row_index)
                {
                    PauliOp row(n_qubits);
                    uint64_t x_bits = 0;
                    uint64_t z_bits = 0;
                    for (size_t q = 0; q < n_qubits; ++q)
                    {
                        if (get_index_bit(x_row_index[q], row_index))
                            x_bits |= (1ULL << q);
                        if (get_index_bit(z_row_index[q], row_index))
                            z_bits |= (1ULL << q);
                    }
                    row.get_x_bits_small() = x_bits;
                    row.get_z_bits_small() = z_bits;
                    row.set_phase(get_index_bit(small_phase_rows, row_index));
                    row.set_rowtype(small_rowtype_rows[row_index]);
                    resulting_rows.push_back(std::move(row));
                }
                return resulting_rows;
            }

            for (const auto &row : rows)
            {
                if (row.is_valid())
                {
                    resulting_rows.push_back(row);
                }
            }

            return resulting_rows;
        }

        void enable_row_index()
        {
            if (row_index_enabled || n_qubits > 64)
                return;

            row_index_enabled = true;
            indexed_row_count = 0;
            x_row_index.assign(n_qubits, {});
            z_row_index.assign(n_qubits, {});
            small_phase_rows.clear();
            small_rowtype_rows.reserve(rows.size());

            for (size_t i = 0; i < rows.size(); ++i)
            {
                if (rows[i].is_valid())
                    append_indexed_row(rows[i]);
            }
            rows.clear();
            rows.shrink_to_fit();
        }

    private:
        void ensure_row_index_capacity(size_t row_index)
        {
            const size_t required_words = row_index / 64 + 1;
            if (small_phase_rows.size() < required_words)
                small_phase_rows.resize(required_words, 0);
            for (size_t q = 0; q < n_qubits; ++q)
            {
                if (x_row_index[q].size() < required_words)
                    x_row_index[q].resize(required_words, 0);
                if (z_row_index[q].size() < required_words)
                    z_row_index[q].resize(required_words, 0);
            }
        }

        void set_index_bit(std::vector<uint64_t> &index, size_t row_index)
        {
            index[row_index / 64] |= (1ULL << (row_index % 64));
        }

        bool get_index_bit(const std::vector<uint64_t> &index, size_t row_index) const
        {
            return (index[row_index / 64] & (1ULL << (row_index % 64))) != 0;
        }

        template <typename Fn>
        void for_each_set_bit(uint64_t mask, Fn fn) const
        {
            while (mask)
            {
                int bit = ctz64(mask);
                fn(static_cast<size_t>(bit));
                mask &= mask - 1;
            }
        }

        void add_row_to_index(size_t row_index, uint64_t x_bits, uint64_t z_bits)
        {
            ensure_row_index_capacity(row_index);
            for_each_set_bit(x_bits, [&](size_t q) {
                if (q < n_qubits)
                    set_index_bit(x_row_index[q], row_index);
            });
            for_each_set_bit(z_bits, [&](size_t q) {
                if (q < n_qubits)
                    set_index_bit(z_row_index[q], row_index);
            });
        }

        void append_indexed_row(const PauliOp &row)
        {
            if (!row.is_valid())
                return;

            const size_t row_index = indexed_row_count++;
            uint64_t x_bits = row.get_x_bits_small();
            uint64_t z_bits = row.get_z_bits_small();

            small_rowtype_rows.push_back(row.get_rowtype());

            ensure_row_index_capacity(row_index);
            if (row.get_phase())
                set_index_bit(small_phase_rows, row_index);
            add_row_to_index(row_index, x_bits, z_bits);
        }

        void materialize_indexed_rows()
        {
            if (!row_index_enabled)
                return;

            rows = get_rows();
            row_index_enabled = false;
            indexed_row_count = 0;
            small_phase_rows.clear();
            small_rowtype_rows.clear();
            x_row_index.clear();
            z_row_index.clear();
            anti_row_scratch.clear();
        }

        void update_index_for_xor(uint64_t x_bits, uint64_t z_bits, size_t words)
        {
            for_each_set_bit(x_bits, [&](size_t q) {
                if (q < n_qubits)
                {
                    for (size_t word = 0; word < words; ++word)
                        x_row_index[q][word] ^= anti_row_scratch[word];
                }
            });
            for_each_set_bit(z_bits, [&](size_t q) {
                if (q < n_qubits)
                {
                    for (size_t word = 0; word < words; ++word)
                        z_row_index[q][word] ^= anti_row_scratch[word];
                }
            });
        }

        void front_multiply_pauli_indexed(const PauliOp &new_pauli)
        {
            const size_t words = indexed_word_count();
            anti_row_scratch.assign(words, 0);

            const uint64_t new_x = new_pauli.get_x_bits_small();
            const uint64_t new_z = new_pauli.get_z_bits_small();
            const uint64_t new_x_only = new_x & ~new_z;
            const uint64_t new_z_only = new_z & ~new_x;
            const uint64_t new_y = new_x & new_z;

            for_each_set_bit(new_x, [&](size_t q) {
                if (q >= n_qubits)
                    return;
                for (size_t word = 0; word < words; ++word)
                    anti_row_scratch[word] ^= z_row_index[q][word];
            });
            for_each_set_bit(new_z, [&](size_t q) {
                if (q >= n_qubits)
                    return;
                for (size_t word = 0; word < words; ++word)
                    anti_row_scratch[word] ^= x_row_index[q][word];
            });

            update_indexed_phases(new_pauli.get_phase(), new_x_only, new_z_only, new_y, words);
            update_index_for_xor(new_x, new_z, words);
        }

        size_t indexed_word_count() const
        {
            return (indexed_row_count + 63) / 64;
        }

        void add_mod4_mask(std::vector<uint64_t> &low,
                           std::vector<uint64_t> &high,
                           size_t word,
                           uint64_t mask)
        {
            uint64_t old_low = low[word];
            low[word] ^= mask;
            high[word] ^= old_low & mask;
        }

        void subtract_mod4_mask(std::vector<uint64_t> &low,
                                std::vector<uint64_t> &high,
                                size_t word,
                                uint64_t mask)
        {
            uint64_t old_low = low[word];
            low[word] ^= mask;
            high[word] ^= ~old_low & mask;
        }

        void update_indexed_phases(bool new_phase,
                                   uint64_t new_x_only,
                                   uint64_t new_z_only,
                                   uint64_t new_y,
                                   size_t words)
        {
            phase_low_scratch.assign(words, 0);
            phase_high_scratch.assign(words, 0);

            auto apply_masks = [&](uint64_t qubit_mask, bool positive_uses_x, bool positive_uses_z) {
                for_each_set_bit(qubit_mask, [&](size_t q) {
                    if (q >= n_qubits)
                        return;
                    for (size_t word = 0; word < words; ++word)
                    {
                        uint64_t x_rows = x_row_index[q][word];
                        uint64_t z_rows = z_row_index[q][word];
                        uint64_t y_rows = x_rows & z_rows;
                        uint64_t x_only_rows = x_rows & ~z_rows;
                        uint64_t z_only_rows = z_rows & ~x_rows;
                        uint64_t positive_rows = positive_uses_x ? x_only_rows : positive_uses_z ? z_only_rows
                                                                                                  : y_rows;
                        uint64_t negative_rows = positive_uses_x ? y_rows : positive_uses_z ? x_only_rows
                                                                                             : z_only_rows;
                        add_mod4_mask(phase_low_scratch, phase_high_scratch, word, positive_rows);
                        subtract_mod4_mask(phase_low_scratch, phase_high_scratch, word, negative_rows);
                    }
                });
            };

            apply_masks(new_x_only, false, false);
            apply_masks(new_z_only, true, false);
            apply_masks(new_y, false, true);

            for (size_t word = 0; word < words; ++word)
            {
                uint64_t phase_toggle = phase_low_scratch[word] & ~phase_high_scratch[word] & anti_row_scratch[word];
                if (new_phase)
                    phase_toggle ^= anti_row_scratch[word];
                small_phase_rows[word] ^= phase_toggle;
            }
        }

        int compute_g_function_small(uint64_t x1_only,
                                     uint64_t z1_only,
                                     uint64_t y1,
                                     const PauliOp &pauli2) const
        {
            return compute_g_function_small(x1_only, z1_only, y1,
                                            pauli2.get_x_bits_small(), pauli2.get_z_bits_small());
        }

        int compute_g_function_small(uint64_t x1_only,
                                     uint64_t z1_only,
                                     uint64_t y1,
                                     uint64_t x2,
                                     uint64_t z2) const
        {
            uint64_t x2_only = x2 & ~z2;
            uint64_t z2_only = z2 & ~x2;
            uint64_t y2 = x2 & z2;

            uint64_t positive_mask = (x1_only & y2) |
                                     (z1_only & x2_only) |
                                     (y1 & z2_only);
            uint64_t negative_mask = (x1_only & z2_only) |
                                     (z1_only & y2) |
                                     (y1 & x2_only);
            int positive = popcount64(positive_mask);
            int negative = popcount64(negative_mask);
            return (positive - negative) & 3;
        }

        int compute_g_function(const PauliOp &pauli1, const PauliOp &pauli2) const
        {
            int g_val = 0;

            if (pauli1.is_small())
            {
                uint64_t x1 = pauli1.get_x_bits_small();
                uint64_t z1 = pauli1.get_z_bits_small();

                uint64_t x1_only = x1 & ~z1;
                uint64_t z1_only = z1 & ~x1;
                uint64_t y1 = x1 & z1;
                g_val = compute_g_function_small(x1_only, z1_only, y1, pauli2);
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
                return commutes_small(row1.get_x_bits_small(), row1.get_z_bits_small(),
                                      row2.get_x_bits_small(), row2.get_z_bits_small());
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

        bool commutes_small(uint64_t x1, uint64_t z1, uint64_t x2, uint64_t z2) const
        {
            uint64_t anti_commute_word = (x1 & z2) ^ (z1 & x2);
            return (popcount64(anti_commute_word) & 1) == 0;
        }

        size_t n_qubits;
        std::vector<PauliOp> rows;
        bool row_index_enabled = false;
        size_t indexed_row_count = 0;
        std::vector<uint64_t> small_phase_rows;
        std::vector<RowType> small_rowtype_rows;
        std::vector<std::vector<uint64_t>> x_row_index;
        std::vector<std::vector<uint64_t>> z_row_index;
        mutable std::vector<uint64_t> anti_row_scratch;
        mutable std::vector<uint64_t> phase_low_scratch;
        mutable std::vector<uint64_t> phase_high_scratch;
    };

} // namespace NWQEC
