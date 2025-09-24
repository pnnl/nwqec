#pragma once

#include "nwqec/core/dag_circuit.hpp"
#include "nwqec/core/pauli_op.hpp"
#include "pass_template.hpp"
#include <vector>
#include <cmath>
#include <map>
#include <set>
#include <algorithm>
#include <optional>
#include <cassert>

namespace NWQEC
{
    /**
     * @brief Specialized pass for single-qubit gate optimization with specific rules
     *
     * This pass applies a specific sequence of optimizations to single-qubit gate sequences:
     * 1. Replace CCX gate with a sequence of T_PAULI gates
     * 2. General optimization (combining rotation gates, canceling pairs)
     * 3. Convert T S sequences to TDG Z
     * 4. Commute all Hadamard gates to the end using rewriting rules:
     *    - HXH = Z, HZH = X, HYH = -Y
     *    - HSH = SX H, H SDG H = SXDG H
     *    - H SX H = SDG H, H SXDG H = S H
     *    - HTH = RX(π/4) H, H TDG H = RX(-π/4) H
     *    - HRY(θ)H = RY(-θ) H, HRZ(θ)H = RZ(-θ) H
     * 5. Final cleanup (T/TDG are already converted to RX(±π/4) in step 4)
     */
    class CRPass : public Pass
    {
    public:
        // Constructor that accepts TranspilationConfig
        CRPass() = default;

        std::string get_name() const override
        {
            return "Clifford Reduction Pass";
        }

        bool run(Circuit &circuit) override
        {

            DAGCircuit dag(circuit);

            auto consecutive_sequences = get_1q_sequences(dag);

            bool modified = false;
            std::set<size_t> optimized_indices;
            std::map<size_t, std::vector<Operation>> replacements;

            // First, handle CCX gates (Step 1 - Replace CCX gates with T_PAULI operations)
            const auto &original_ops = circuit.get_operations();
            for (size_t i = 0; i < original_ops.size(); i++)
            {
                if (original_ops[i].get_type() == Operation::Type::CCX)
                {
                    auto qubits = original_ops[i].get_qubits();
                    auto ccx_ops = create_ccx_t_ops(qubits[0], qubits[1], qubits[2], circuit.get_num_qubits());
                    replacements[i] = ccx_ops;
                    optimized_indices.insert(i);
                    modified = true;
                }
            }

            for (const auto &sequence : consecutive_sequences)
            {
                if (sequence.empty())
                    continue;

                size_t first_idx = sequence[0].first;

                std::vector<Operation> ops_only;
                for (const auto &[idx, op] : sequence)
                {
                    ops_only.push_back(op);
                    optimized_indices.insert(idx);
                }

                // Apply the four-step optimization process
                auto optimized = gate_merging(ops_only); // Step 1, general optimization

                optimized = convert_t_s_to_tdg_z(optimized); // Step 2, remove S

                optimized = commute_hadamards_to_end(optimized); // Step 3, remove H

                // Check if any modification occurred
                if (optimized.size() != ops_only.size() ||
                    !std::equal(optimized.begin(), optimized.end(), ops_only.begin(),
                                [](const Operation &a, const Operation &b)
                                {
                                    return a.get_type() == b.get_type();
                                }))
                {
                    modified = true;
                }

                replacements[first_idx] = optimized;
            }

            // Rebuild circuit with optimizations
            if (modified)
            {
                const auto &original_ops = circuit.get_operations();
                std::vector<Operation> final_operations;

                for (size_t i = 0; i < original_ops.size(); i++)
                {
                    if (replacements.count(i))
                    {
                        for (const auto &opt_op : replacements[i])
                        {
                            final_operations.push_back(opt_op);
                        }
                    }
                    else if (optimized_indices.count(i) == 0)
                    {
                        final_operations.push_back(original_ops[i]);
                    }
                }

                // Create new circuit
                Circuit new_circuit;

                // Copy register information
                new_circuit.add_qreg("q", circuit.get_num_qubits());
                new_circuit.add_creg("c", circuit.get_num_bits());

                // Add operations
                for (const auto &op : final_operations)
                {
                    new_circuit.add_operation(op);
                }

                circuit = std::move(new_circuit);
            }

            return modified;
        }

    private:
        // Create T_PAULI operations for CCX gate replacement
        std::vector<Operation> create_ccx_t_ops(size_t q0, size_t q1, size_t q2, size_t total_qubits) const
        {
            // Get PauliOp stabilizers from the new unified method
            auto pauliops = PauliOp::create_ccx_ops(q0, q1, q2, total_qubits);

            std::vector<Operation> t_ops;
            t_ops.reserve(pauliops.size());

            // Convert each PauliOp to a T_PAULI Operation
            for (const auto &op : pauliops)
            {
                t_ops.push_back(Operation(Operation::Type::T_PAULI, {}, {}, {}, op));
            }

            return t_ops;
        }

        // Check if an operation is a single-qubit gate
        bool is_single_qubit_gate(const Operation &op) const
        {
            return op.get_qubits().size() == 1 &&
                   op.get_type() != Operation::Type::MEASURE &&
                   op.get_type() != Operation::Type::RESET &&
                   op.get_type() != Operation::Type::BARRIER;
        }

        // Check if two gates can cancel each other
        bool gates_cancel(Operation::Type gate1, Operation::Type gate2) const
        {
            // Self-inverse gates
            if (gate1 == gate2 && (gate1 == Operation::Type::X ||
                                   gate1 == Operation::Type::Y ||
                                   gate1 == Operation::Type::Z ||
                                   gate1 == Operation::Type::H))
            {
                return true;
            }

            // Inverse pairs
            if ((gate1 == Operation::Type::S && gate2 == Operation::Type::SDG) ||
                (gate1 == Operation::Type::SDG && gate2 == Operation::Type::S) ||
                (gate1 == Operation::Type::T && gate2 == Operation::Type::TDG) ||
                (gate1 == Operation::Type::TDG && gate2 == Operation::Type::T) ||
                (gate1 == Operation::Type::SX && gate2 == Operation::Type::SXDG) ||
                (gate1 == Operation::Type::SXDG && gate2 == Operation::Type::SX) ||
                // Note: P4 with x_rotation inverse pairs would need Operation objects, not just types
                false) // Placeholder - P4 inverse pairs need more complex logic
            {
                return true;
            }

            return false;
        }

        // Step 1: General optimization - combine rotation gates and remove canceling pairs
        std::vector<Operation> gate_merging(const std::vector<Operation> &sequence) const
        {
            std::vector<Operation> result;

            for (size_t i = 0; i < sequence.size(); i++)
            {
                const Operation &current = sequence[i];

                // Check for canceling pairs first
                if (i + 1 < sequence.size() &&
                    gates_cancel(current.get_type(), sequence[i + 1].get_type()))
                {
                    i++; // Skip both operations
                    continue;
                }

                // Check for merging pairs (gates that combine into other gates)
                if (i + 1 < sequence.size())
                {
                    auto merged_gate = try_merge_gates(current, sequence[i + 1]);
                    if (merged_gate.has_value())
                    {
                        result.push_back(merged_gate.value());
                        i++; // Skip the second gate as it's been merged
                        continue;
                    }
                }

                result.push_back(current);
            }

            // Apply iterative merging until no more changes occur
            bool changed = true;
            while (changed)
            {
                changed = false;
                std::vector<Operation> new_result;

                for (size_t i = 0; i < result.size(); i++)
                {
                    const Operation &current = result[i];

                    // Check for canceling pairs first
                    if (i + 1 < result.size() &&
                        gates_cancel(current.get_type(), result[i + 1].get_type()))
                    {
                        i++; // Skip both operations
                        changed = true;
                        continue;
                    }

                    // Check for merging pairs (gates that combine into other gates)
                    if (i + 1 < result.size())
                    {
                        auto merged_gate = try_merge_gates(current, result[i + 1]);
                        if (merged_gate.has_value())
                        {
                            new_result.push_back(merged_gate.value());
                            i++; // Skip the second gate as it's been merged
                            changed = true;
                            continue;
                        }
                    }

                    new_result.push_back(current);
                }

                result = std::move(new_result);
            }

            return result;
        }

        // Try to merge two consecutive gates into a single gate
        std::optional<Operation> try_merge_gates(const Operation &gate1, const Operation &gate2) const
        {
            // Only merge gates on the same qubit
            if (gate1.get_qubits() != gate2.get_qubits())
            {
                return std::nullopt;
            }

            Operation::Type type1 = gate1.get_type();
            Operation::Type type2 = gate2.get_type();

            // T T = S
            if (type1 == Operation::Type::T && type2 == Operation::Type::T)
            {
                return Operation(Operation::Type::S, gate1.get_qubits());
            }

            // TDG TDG = SDG
            if (type1 == Operation::Type::TDG && type2 == Operation::Type::TDG)
            {
                return Operation(Operation::Type::SDG, gate1.get_qubits());
            }

            // S S = Z
            if (type1 == Operation::Type::S && type2 == Operation::Type::S)
            {
                return Operation(Operation::Type::Z, gate1.get_qubits());
            }

            // SDG SDG = Z
            if (type1 == Operation::Type::SDG && type2 == Operation::Type::SDG)
            {
                return Operation(Operation::Type::Z, gate1.get_qubits());
            }

            // T TDG T = T (T TDG cancels, leaving T)
            // This is handled by the cancellation logic, not merging

            // S TDG = T (since S = T T, so S TDG = T T TDG = T)
            if (type1 == Operation::Type::S && type2 == Operation::Type::TDG)
            {
                return Operation(Operation::Type::T, gate1.get_qubits());
            }

            // TDG S = TDG (since S = T T, so TDG S = TDG T T = T T = S, wait that's wrong)
            // Actually: TDG S = TDG T T = (TDG T) T = I T = T
            if (type1 == Operation::Type::TDG && type2 == Operation::Type::S)
            {
                return Operation(Operation::Type::T, gate1.get_qubits());
            }

            // SDG T = TDG (since SDG = TDG TDG, so SDG T = TDG TDG T = TDG (TDG T) = TDG I = TDG)
            if (type1 == Operation::Type::SDG && type2 == Operation::Type::T)
            {
                return Operation(Operation::Type::TDG, gate1.get_qubits());
            }

            // T SDG = TDG (since SDG = TDG TDG, so T SDG = T TDG TDG = (T TDG) TDG = I TDG = TDG)
            if (type1 == Operation::Type::T && type2 == Operation::Type::SDG)
            {
                return Operation(Operation::Type::TDG, gate1.get_qubits());
            }

            // No merging rule found
            return std::nullopt;
        }

        // Step 2: Convert T S sequences to TDG Z
        std::vector<Operation> convert_t_s_to_tdg_z(const std::vector<Operation> &sequence) const
        {
            std::vector<Operation> result;

            for (size_t i = 0; i < sequence.size(); i++)
            {
                const Operation &current = sequence[i];

                // Look for T S pattern
                if (current.get_type() == Operation::Type::T &&
                    i + 1 < sequence.size() &&
                    sequence[i + 1].get_type() == Operation::Type::S)
                {
                    // Replace T S with TDG Z
                    result.push_back(Operation(Operation::Type::TDG, current.get_qubits()));
                    result.push_back(Operation(Operation::Type::Z, current.get_qubits()));
                    i++; // Skip the S gate
                }
                else
                {
                    result.push_back(current);
                }
            }

            return result;
        }

        // Step 3: Commute Hadamard gates to the end using gate rewriting rules
        std::vector<Operation> commute_hadamards_to_end(const std::vector<Operation> &sequence) const
        {
            std::vector<Operation> result;
            size_t hadamard_count = 0;

            // Process each gate, applying Hadamard commutation rules
            for (const auto &op : sequence)
            {
                if (op.get_type() == Operation::Type::H)
                {
                    hadamard_count++;
                }
                else
                {
                    // Apply Hadamard rewriting rules based on the number of Hadamards
                    auto rewritten_gate = apply_hadamard_rewriting_rules(op, hadamard_count);
                    result.push_back(rewritten_gate);
                }
            }

            // Add final Hadamard gates (combined if even number, single if odd)
            if (hadamard_count % 2 == 1)
            {
                // Odd number of Hadamards - add one H at the end
                if (!sequence.empty())
                {
                    result.push_back(Operation(Operation::Type::H, sequence[0].get_qubits()));
                }
            }
            // Even number of Hadamards cancel out - add nothing

            return result;
        }

        // Apply Hadamard rewriting rules to move H through other gates
        Operation apply_hadamard_rewriting_rules(const Operation &op, size_t num_hadamards) const
        {
            if (num_hadamards == 0)
                return op;

            Operation::Type new_type = op.get_type();
            bool odd_hadamards = (num_hadamards % 2 == 1);

            // Apply rewriting rules only for odd number of Hadamards
            // Even number of Hadamards cancel out
            if (odd_hadamards)
            {
                switch (op.get_type())
                {
                // Pauli gate rewriting: HXH = Z, HZH = X, HYH = -Y (ignore phase)
                case Operation::Type::X:
                    new_type = Operation::Type::Z;
                    break;
                case Operation::Type::Z:
                    new_type = Operation::Type::X;
                    break;
                case Operation::Type::Y:
                    new_type = Operation::Type::Y; // Global phase -1 ignored
                    break;

                // S gate rewriting: HSH = SX H → SX (H moved to end)
                case Operation::Type::S:
                    new_type = Operation::Type::SX;
                    break;
                case Operation::Type::SDG:
                    new_type = Operation::Type::SXDG;
                    break;

                // SX gate rewriting: H SX H = SDG H → SDG (H moved to end)
                case Operation::Type::SX:
                    new_type = Operation::Type::S;
                    break;
                case Operation::Type::SXDG:
                    new_type = Operation::Type::SDG;
                    break;

                // T gate rewriting: HTH = RX(π/4), H TDG H = RX(-π/4)
                case Operation::Type::T:
                    // T = RZ(π/4) -> H T H = RX(π/4) = P4 with x_rotation=true, dagger=false
                    return Operation(Operation::Type::P4, op.get_qubits(), op.get_parameters(),
                                     op.get_bits(), op.get_pauli_op(), false, true);
                case Operation::Type::TDG:
                    // TDG = RZ(-π/4) -> H TDG H = RX(-π/4) = P4 with x_rotation=true, dagger=true
                    return Operation(Operation::Type::P4, op.get_qubits(), op.get_parameters(),
                                     op.get_bits(), op.get_pauli_op(), true, true);

                // P4 gate rewriting: Handle various P4 configurations
                case Operation::Type::P4:
                {
                    // P4 with x_rotation=false is RZ(±π/4) -> becomes RX(±π/4)
                    // P4 with x_rotation=true is RX(±π/4) -> becomes RZ(±π/4)
                    bool new_x_rotation = !op.get_x_rotation(); // Flip the x_rotation flag
                    return Operation(Operation::Type::P4, op.get_qubits(), op.get_parameters(),
                                     op.get_bits(), op.get_pauli_op(), op.get_dagger(), new_x_rotation);
                }

                // P8 gate rewriting: Handle various P8 configurations
                case Operation::Type::P8:
                {
                    // P8 with x_rotation=false is RZ(±π/8) -> becomes RX(±π/8)
                    // P8 with x_rotation=true is RX(±π/8) -> becomes RZ(±π/8)
                    bool new_x_rotation = !op.get_x_rotation(); // Flip the x_rotation flag
                    return Operation(Operation::Type::P8, op.get_qubits(), op.get_parameters(),
                                     op.get_bits(), op.get_pauli_op(), op.get_dagger(), new_x_rotation);
                }

                // P16 gate rewriting: Handle various P16 configurations
                case Operation::Type::P16:
                {
                    // P16 with x_rotation=false is RZ(±π/16) -> becomes RX(±π/16)
                    // P16 with x_rotation=true is RX(±π/16) -> becomes RZ(±π/16)
                    bool new_x_rotation = !op.get_x_rotation(); // Flip the x_rotation flag
                    return Operation(Operation::Type::P16, op.get_qubits(), op.get_parameters(),
                                     op.get_bits(), op.get_pauli_op(), op.get_dagger(), new_x_rotation);
                }

                // Rotation gates: HRX(θ)H = RZ(θ), HRY(θ)H = RY(-θ), HRZ(θ)H = RX(θ)
                case Operation::Type::RX:
                {
                    // RX(θ) -> RZ(θ) under Hadamard conjugation
                    auto params = op.get_parameters();
                    assert(params.size() == 1 && "RX gate must have exactly one parameter");
                    return Operation(Operation::Type::RZ, op.get_qubits(), params);
                }
                case Operation::Type::RY:
                {
                    // RY(θ) -> RY(-θ) under Hadamard conjugation
                    auto params = op.get_parameters();
                    assert(params.size() == 1 && "RY gate must have exactly one parameter");
                    params[0] = -params[0];
                    return Operation(Operation::Type::RY, op.get_qubits(), params);
                }
                case Operation::Type::RZ:
                {
                    // RZ(θ) -> RX(θ) under Hadamard conjugation
                    auto params = op.get_parameters();
                    assert(params.size() == 1 && "RZ gate must have exactly one parameter");
                    return Operation(Operation::Type::RX, op.get_qubits(), params);
                }

                default:
                    // No rewriting rule defined - keep original gate
                    break;
                }
            }

            return Operation(new_type, op.get_qubits(), op.get_parameters());
        }

        // Extract consecutive single-qubit gate sequences using DAG traversal
        std::vector<std::vector<std::pair<size_t, Operation>>> get_1q_sequences(const DAGCircuit &dag) const
        {
            std::vector<std::vector<std::pair<size_t, Operation>>> sequences;
            std::set<size_t> visited;
            const auto &operations = dag.get_operations();

            for (size_t start_idx = 0; start_idx < operations.size(); start_idx++)
            {
                if (visited.count(start_idx) || !is_single_qubit_gate(operations[start_idx]))
                {
                    continue;
                }

                std::vector<std::pair<size_t, Operation>> sequence;
                size_t current_qubit = operations[start_idx].get_qubits()[0];
                size_t current_idx = start_idx;

                while (current_idx < operations.size() &&
                       !visited.count(current_idx) &&
                       is_single_qubit_gate(operations[current_idx]) &&
                       operations[current_idx].get_qubits()[0] == current_qubit)
                {
                    sequence.push_back({current_idx, operations[current_idx]});
                    visited.insert(current_idx);

                    const auto &successors = dag.get_successors(current_idx);
                    size_t next_idx = std::numeric_limits<size_t>::max();

                    for (const auto &succ : successors)
                    {
                        if (succ.qubit == current_qubit)
                        {
                            next_idx = succ.node;
                            break;
                        }
                    }

                    if (next_idx == std::numeric_limits<size_t>::max() ||
                        next_idx >= operations.size() ||
                        !is_single_qubit_gate(operations[next_idx]) ||
                        operations[next_idx].get_qubits()[0] != current_qubit)
                    {
                        break;
                    }

                    current_idx = next_idx;
                }

                // Add sequences with at least one gate (apply optimization even to single gates)
                if (!sequence.empty())
                {
                    sequences.push_back(sequence);
                }
            }

            return sequences;
        }
    };
}
