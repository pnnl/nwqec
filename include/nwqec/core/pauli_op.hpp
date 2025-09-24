#pragma once

#include <vector>
#include <string>
#include <cstdint>
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace NWQEC
{
    enum class RowType
    {
        T,
        S,
        Z
    };

    class PauliOp
    {
    private:
        // Optimized storage: single uint64_t for ≤64 qubits, vector for larger circuits
        uint64_t x_bits_small;
        uint64_t z_bits_small;
        std::vector<uint64_t> x_bits_large;
        std::vector<uint64_t> z_bits_large;
        bool phase;
        size_t num_qubits;
        bool is_small_circuit; // true if ≤64 qubits
        size_t weight;         // Total number of X and Z bits (for sorting)
        bool valid;            // Whether this row is still valid
        RowType rowtype;       // Type of the row: T, S, or Z

    public:
        PauliOp(size_t qubits) : x_bits_small(0), z_bits_small(0), phase(false), num_qubits(qubits), weight(0), valid(true), rowtype(RowType::T)
        {
            is_small_circuit = (qubits <= 64);
            if (!is_small_circuit)
            {
                size_t num_words = (qubits + 63) / 64;
                x_bits_large.resize(num_words, 0);
                z_bits_large.resize(num_words, 0);
            }
        }

        // Default constructor
        PauliOp() : PauliOp(0) {}

        // Getters for TabRow compatibility
        bool get_phase() const { return phase; }
        void set_phase(bool p) { phase = p; }
        size_t get_num_qubits() const { return num_qubits; }
        size_t get_weight() const { return weight; }
        bool is_valid() const { return valid; }
        void set_valid(bool v) { valid = v; }
        RowType get_rowtype() const { return rowtype; }
        void set_rowtype(RowType type) { rowtype = type; }
        bool is_small() const { return is_small_circuit; }

        // Direct bit access for optimizer compatibility
        uint64_t& get_x_bits_small() { return x_bits_small; }
        uint64_t& get_z_bits_small() { return z_bits_small; }
        uint64_t get_x_bits_small() const { return x_bits_small; }
        uint64_t get_z_bits_small() const { return z_bits_small; }
        const std::vector<uint64_t>& get_x_bits_large() const { return x_bits_large; }
        const std::vector<uint64_t>& get_z_bits_large() const { return z_bits_large; }
        std::vector<uint64_t>& get_x_bits_large() { return x_bits_large; }
        std::vector<uint64_t>& get_z_bits_large() { return z_bits_large; }

        void from_string(const std::string &pauli_str)
        {
            phase = (pauli_str[0] == '-');
            weight = 0;

            if (is_small_circuit)
            {
                // Fast path for ≤64 qubits - no vector overhead
                x_bits_small = 0;
                z_bits_small = 0;

                for (size_t i = 1; i < pauli_str.size(); ++i)
                {
                    size_t qubit = i - 1;
                    if (qubit >= num_qubits)
                        continue;

                    uint64_t mask = 1ULL << qubit;

                    switch (pauli_str[i])
                    {
                    case 'X':
                    case 'x':
                        x_bits_small |= mask;
                        weight++;
                        break;
                    case 'Y':
                    case 'y':
                        x_bits_small |= mask;
                        z_bits_small |= mask;
                        weight++; // Y = X + Z
                        break;
                    case 'Z':
                    case 'z':
                        z_bits_small |= mask;
                        weight++;
                        break;
                    default:
                        break;
                    }
                }
            }
            else
            {
                // Large circuit path
                std::fill(x_bits_large.begin(), x_bits_large.end(), 0);
                std::fill(z_bits_large.begin(), z_bits_large.end(), 0);

                for (size_t i = 1; i < pauli_str.size(); ++i)
                {
                    size_t qubit = i - 1;
                    if (qubit >= num_qubits)
                        continue;

                    size_t word_idx = qubit / 64;
                    size_t bit_idx = qubit % 64;
                    uint64_t mask = 1ULL << bit_idx;

                    switch (pauli_str[i])
                    {
                    case 'X':
                    case 'x':
                        x_bits_large[word_idx] |= mask;
                        weight++;
                        break;
                    case 'Y':
                    case 'y':
                        x_bits_large[word_idx] |= mask;
                        z_bits_large[word_idx] |= mask;
                        weight++; // Y = X + Z
                        break;
                    case 'Z':
                    case 'z':
                        z_bits_large[word_idx] |= mask;
                        weight++;
                        break;
                    default:
                        break;
                    }
                }
            }
        }

        std::string to_string() const
        {
            std::string result;
            result += phase ? '-' : '+';

            if (is_small_circuit)
            {
                for (size_t qubit = 0; qubit < num_qubits; ++qubit)
                {
                    uint64_t mask = 1ULL << qubit;
                    bool x_bit = (x_bits_small & mask) != 0;
                    bool z_bit = (z_bits_small & mask) != 0;

                    if (x_bit && z_bit)
                        result += 'Y';
                    else if (x_bit)
                        result += 'X';
                    else if (z_bit)
                        result += 'Z';
                    else
                        result += 'I';
                }
            }
            else
            {
                for (size_t qubit = 0; qubit < num_qubits; ++qubit)
                {
                    size_t word_idx = qubit / 64;
                    size_t bit_idx = qubit % 64;
                    uint64_t mask = 1ULL << bit_idx;

                    bool x_bit = (x_bits_large[word_idx] & mask) != 0;
                    bool z_bit = (z_bits_large[word_idx] & mask) != 0;

                    if (x_bit && z_bit)
                        result += 'Y';
                    else if (x_bit)
                        result += 'X';
                    else if (z_bit)
                        result += 'Z';
                    else
                        result += 'I';
                }
            }

            return result;
        }

        std::string to_string(size_t n_qubits) const
        {
            std::string result(n_qubits + 1, 'I');
            result[0] = phase ? '-' : '+';

            if (is_small_circuit)
            {
                for (size_t i = 0; i < std::min(n_qubits, num_qubits); ++i)
                {
                    uint64_t mask = 1ULL << i;
                    bool x_bit = (x_bits_small & mask) != 0;
                    bool z_bit = (z_bits_small & mask) != 0;

                    if (x_bit && z_bit)
                        result[i + 1] = 'Y';
                    else if (x_bit)
                        result[i + 1] = 'X';
                    else if (z_bit)
                        result[i + 1] = 'Z';
                }
            }
            else
            {
                for (size_t i = 0; i < std::min(n_qubits, num_qubits); ++i)
                {
                    size_t word_idx = i / 64;
                    size_t bit_idx = i % 64;
                    uint64_t mask = 1ULL << bit_idx;

                    bool x_bit = (x_bits_large[word_idx] & mask) != 0;
                    bool z_bit = (z_bits_large[word_idx] & mask) != 0;

                    if (x_bit && z_bit)
                        result[i + 1] = 'Y';
                    else if (x_bit)
                        result[i + 1] = 'X';
                    else if (z_bit)
                        result[i + 1] = 'Z';
                }
            }

            return result;
        }

        // PauliStr compatibility methods
        std::vector<size_t> get_x_indices() const
        {
            std::vector<size_t> indices;
            
            if (is_small_circuit)
            {
                for (size_t i = 0; i < num_qubits; ++i)
                {
                    if (x_bits_small & (1ULL << i))
                        indices.push_back(i);
                }
            }
            else
            {
                for (size_t i = 0; i < num_qubits; ++i)
                {
                    size_t word_idx = i / 64;
                    size_t bit_idx = i % 64;
                    if (x_bits_large[word_idx] & (1ULL << bit_idx))
                        indices.push_back(i);
                }
            }
            
            return indices;
        }

        std::vector<size_t> get_z_indices() const
        {
            std::vector<size_t> indices;
            
            if (is_small_circuit)
            {
                for (size_t i = 0; i < num_qubits; ++i)
                {
                    if (z_bits_small & (1ULL << i))
                        indices.push_back(i);
                }
            }
            else
            {
                for (size_t i = 0; i < num_qubits; ++i)
                {
                    size_t word_idx = i / 64;
                    size_t bit_idx = i % 64;
                    if (z_bits_large[word_idx] & (1ULL << bit_idx))
                        indices.push_back(i);
                }
            }
            
            return indices;
        }

        bool has_x(size_t qubit) const
        {
            if (qubit >= num_qubits) return false;
            
            if (is_small_circuit)
            {
                return (x_bits_small & (1ULL << qubit)) != 0;
            }
            else
            {
                size_t word_idx = qubit / 64;
                size_t bit_idx = qubit % 64;
                return (x_bits_large[word_idx] & (1ULL << bit_idx)) != 0;
            }
        }

        bool has_z(size_t qubit) const
        {
            if (qubit >= num_qubits) return false;
            
            if (is_small_circuit)
            {
                return (z_bits_small & (1ULL << qubit)) != 0;
            }
            else
            {
                size_t word_idx = qubit / 64;
                size_t bit_idx = qubit % 64;
                return (z_bits_large[word_idx] & (1ULL << bit_idx)) != 0;
            }
        }

        void add_x(size_t qubit)
        {
            if (qubit >= num_qubits) return;
            if (has_x(qubit)) return;
            
            if (is_small_circuit)
            {
                x_bits_small |= (1ULL << qubit);
            }
            else
            {
                size_t word_idx = qubit / 64;
                size_t bit_idx = qubit % 64;
                x_bits_large[word_idx] |= (1ULL << bit_idx);
            }
            weight++;
        }

        void add_z(size_t qubit)
        {
            if (qubit >= num_qubits) return;
            if (has_z(qubit)) return;
            
            if (is_small_circuit)
            {
                z_bits_small |= (1ULL << qubit);
            }
            else
            {
                size_t word_idx = qubit / 64;
                size_t bit_idx = qubit % 64;
                z_bits_large[word_idx] |= (1ULL << bit_idx);
            }
            weight++;
        }

        // Legacy property access for backward compatibility
        bool r() const { return phase; }
        void set_r(bool value) { phase = value; }
        std::vector<size_t> x_indices() const { return get_x_indices(); }
        std::vector<size_t> z_indices() const { return get_z_indices(); }

        // CCX stabilizer creation (moved from PBC passes)
        static std::vector<PauliOp> create_ccx_ops(size_t q0, size_t q1, size_t q2, size_t n_qubits)
        {
            std::vector<PauliOp> stabs;
            stabs.reserve(7);

            // Stabilizer 0: + IIX  (r=false, x2)
            PauliOp s0(n_qubits);
            s0.set_r(false);
            s0.add_x(q2);
            stabs.push_back(s0);

            // Stabilizer 1: - ZZI  (r=true, z0, z1)
            PauliOp s1(n_qubits);
            s1.set_r(true);
            s1.add_z(q0);
            s1.add_z(q1);
            stabs.push_back(s1);

            // Stabilizer 2: + ZII  (r=false, z0)
            PauliOp s2(n_qubits);
            s2.set_r(false);
            s2.add_z(q0);
            stabs.push_back(s2);

            // Stabilizer 3: - ZIX  (r=true, z0, x2)
            PauliOp s3(n_qubits);
            s3.set_r(true);
            s3.add_z(q0);
            s3.add_x(q2);
            stabs.push_back(s3);

            // Stabilizer 4: + IZI  (r=false, z1)
            PauliOp s4(n_qubits);
            s4.set_r(false);
            s4.add_z(q1);
            stabs.push_back(s4);

            // Stabilizer 5: + ZZX  (r=false, z0, z1, x2)
            PauliOp s5(n_qubits);
            s5.set_r(false);
            s5.add_z(q0);
            s5.add_z(q1);
            s5.add_x(q2);
            stabs.push_back(s5);

            // Stabilizer 6: - IZX  (r=true, z1, x2)
            PauliOp s6(n_qubits);
            s6.set_r(true);
            s6.add_z(q1);
            s6.add_x(q2);
            stabs.push_back(s6);

            return stabs;
        }
    };

} // namespace NWQEC
