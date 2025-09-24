#pragma once

#include "nwqec/core/dag_circuit.hpp"
#include "pass_template.hpp"
#include <vector>
#include <cmath>
#include <map>
#include <set>

namespace NWQEC
{
    /**
     * @brief Pass to optimize consecutive single-qubit gate sequences
     *
     * This pass identifies consecutive single-qubit gates on the same qubit and optimizes them
     * by combining, canceling, or replacing with more efficient gate sequences.
     */
    class GateFusionPass : public Pass
    {
    public:
        GateFusionPass() = default;

    private:
        const double TOLERANCE = 1e-10;

        // Check if an operation is a single-qubit gate
        bool is_single_qubit_gate(const Operation &op) const
        {
            return op.get_qubits().size() == 1 &&
                   op.get_type() != Operation::Type::MEASURE &&
                   op.get_type() != Operation::Type::RESET &&
                   op.get_type() != Operation::Type::BARRIER;
        }

        // Check if two single-qubit gates can be combined
        bool gates_commute_and_combine(Operation::Type gate1, Operation::Type gate2) const
        {
            // Pauli gates commute with themselves
            if (gate1 == gate2 && (gate1 == Operation::Type::X ||
                                   gate1 == Operation::Type::Y ||
                                   gate1 == Operation::Type::Z))
            {
                return true;
            }

            // S and SDG cancel out
            if ((gate1 == Operation::Type::S && gate2 == Operation::Type::SDG) ||
                (gate1 == Operation::Type::SDG && gate2 == Operation::Type::S))
            {
                return true;
            }

            // T and TDG cancel out
            if ((gate1 == Operation::Type::T && gate2 == Operation::Type::TDG) ||
                (gate1 == Operation::Type::TDG && gate2 == Operation::Type::T))
            {
                return true;
            }

            // SX and SXDG cancel out
            if ((gate1 == Operation::Type::SX && gate2 == Operation::Type::SXDG) ||
                (gate1 == Operation::Type::SXDG && gate2 == Operation::Type::SX))
            {
                return true;
            }

            return false;
        }

        // Check if a gate is its own inverse (self-canceling)
        bool is_self_inverse(Operation::Type gate) const
        {
            return gate == Operation::Type::X || gate == Operation::Type::Y ||
                   gate == Operation::Type::Z || gate == Operation::Type::H;
        }

        // Combine rotation gates with the same axis
        std::vector<Operation> combine_rotation_gates(const std::vector<Operation> &sequence) const
        {
            std::vector<Operation> result;

            for (size_t i = 0; i < sequence.size();)
            {
                const Operation &current = sequence[i];

                // Handle RZ gates
                if (current.get_type() == Operation::Type::RZ)
                {
                    double total_angle = current.get_parameters()[0];
                    size_t j = i + 1;

                    // Combine consecutive RZ gates
                    while (j < sequence.size() && sequence[j].get_type() == Operation::Type::RZ)
                    {
                        total_angle += sequence[j].get_parameters()[0];
                        j++;
                    }

                    // Normalize angle to [0, 2π)
                    total_angle = std::fmod(total_angle, 2 * M_PI);
                    if (total_angle < 0)
                        total_angle += 2 * M_PI;

                    // Skip if close to identity
                    if (std::abs(total_angle) > TOLERANCE && std::abs(total_angle - 2 * M_PI) > TOLERANCE)
                    {
                        result.push_back(Operation(Operation::Type::RZ, current.get_qubits(), {total_angle}));
                    }

                    i = j;
                }
                // Handle RX gates
                else if (current.get_type() == Operation::Type::RX)
                {
                    double total_angle = current.get_parameters()[0];
                    size_t j = i + 1;

                    // Combine consecutive RX gates
                    while (j < sequence.size() && sequence[j].get_type() == Operation::Type::RX)
                    {
                        total_angle += sequence[j].get_parameters()[0];
                        j++;
                    }

                    // Normalize angle
                    total_angle = std::fmod(total_angle, 2 * M_PI);
                    if (total_angle < 0)
                        total_angle += 2 * M_PI;

                    // Skip if close to identity
                    if (std::abs(total_angle) > TOLERANCE && std::abs(total_angle - 2 * M_PI) > TOLERANCE)
                    {
                        result.push_back(Operation(Operation::Type::RX, current.get_qubits(), {total_angle}));
                    }

                    i = j;
                }
                // Handle RY gates
                else if (current.get_type() == Operation::Type::RY)
                {
                    double total_angle = current.get_parameters()[0];
                    size_t j = i + 1;

                    // Combine consecutive RY gates
                    while (j < sequence.size() && sequence[j].get_type() == Operation::Type::RY)
                    {
                        total_angle += sequence[j].get_parameters()[0];
                        j++;
                    }

                    // Normalize angle
                    total_angle = std::fmod(total_angle, 2 * M_PI);
                    if (total_angle < 0)
                        total_angle += 2 * M_PI;

                    // Skip if close to identity
                    if (std::abs(total_angle) > TOLERANCE && std::abs(total_angle - 2 * M_PI) > TOLERANCE)
                    {
                        result.push_back(Operation(Operation::Type::RY, current.get_qubits(), {total_angle}));
                    }

                    i = j;
                }
                else
                {
                    result.push_back(current);
                    i++;
                }
            }

            return result;
        }

        // Remove canceling gate pairs
        std::vector<Operation> remove_canceling_gates(const std::vector<Operation> &sequence) const
        {
            std::vector<Operation> result;

            for (size_t i = 0; i < sequence.size(); i++)
            {
                const Operation &current = sequence[i];

                // Check for self-canceling pairs (e.g., X X = I)
                if (i + 1 < sequence.size() && is_self_inverse(current.get_type()) &&
                    current.get_type() == sequence[i + 1].get_type())
                {
                    i++; // Skip both operations
                    continue;
                }

                // Check for inverse pairs (e.g., S SDG = I)
                if (i + 1 < sequence.size() &&
                    gates_commute_and_combine(current.get_type(), sequence[i + 1].get_type()))
                {
                    i++; // Skip both operations
                    continue;
                }

                result.push_back(current);
            }

            return result;
        }

        // Extract consecutive single-qubit gate sequences using DAG traversal
        std::vector<std::vector<std::pair<size_t, Operation>>> extract_consecutive_single_qubit_sequences(const DAGCircuit &dag) const
        {
            std::vector<std::vector<std::pair<size_t, Operation>>> sequences;
            std::set<size_t> visited;
            const auto &operations = dag.get_operations();

            // For each operation, try to build a consecutive single-qubit sequence starting from it
            for (size_t start_idx = 0; start_idx < operations.size(); start_idx++)
            {
                if (visited.count(start_idx) || !is_single_qubit_gate(operations[start_idx]))
                {
                    continue;
                }

                // Start building a sequence from this operation
                std::vector<std::pair<size_t, Operation>> sequence;
                size_t current_qubit = operations[start_idx].get_qubits()[0];
                size_t current_idx = start_idx;

                // Follow the chain of consecutive single-qubit gates on the same qubit
                while (current_idx < operations.size() &&
                       !visited.count(current_idx) &&
                       is_single_qubit_gate(operations[current_idx]) &&
                       operations[current_idx].get_qubits()[0] == current_qubit)
                {
                    sequence.push_back({current_idx, operations[current_idx]});
                    visited.insert(current_idx);

                    // Find the next operation on this qubit using DAG successors
                    const auto &successors = dag.get_successors(current_idx);
                    size_t next_idx = std::numeric_limits<size_t>::max();

                    // Look for the immediate successor that operates on the same qubit
                    for (const auto &succ : successors)
                    {
                        if (succ.qubit == current_qubit)
                        {
                            next_idx = succ.node;
                            break;
                        }
                    }

                    // If no valid successor found, or it's not a single-qubit gate on same qubit, break
                    if (next_idx == std::numeric_limits<size_t>::max() ||
                        next_idx >= operations.size() ||
                        !is_single_qubit_gate(operations[next_idx]) ||
                        operations[next_idx].get_qubits()[0] != current_qubit)
                    {
                        break;
                    }

                    current_idx = next_idx;
                }

                // Only add sequences with more than one gate (worth optimizing)
                if (sequence.size() > 1)
                {
                    sequences.push_back(sequence);
                }
            }

            return sequences;
        }

    public:
        std::string get_name() const override
        {
            return "Optimize Single Qubit Pass";
        }

        bool run(Circuit &circuit) override
        {
            // Convert to DAG for dependency analysis
            DAGCircuit dag_circuit(circuit);

            // Find consecutive single-qubit gate sequences
            auto gate_sequences = extract_consecutive_single_qubit_sequences(dag_circuit);

            if (gate_sequences.empty())
                return false;

            bool circuit_modified = false;
            auto optimization_plan = create_optimization_plan(gate_sequences);

            if (optimization_plan.has_optimizations)
            {
                circuit_modified = true;
                apply_optimizations(circuit, optimization_plan);
            }

            return circuit_modified;
        }

    private:
        struct OptimizationPlan
        {
            bool has_optimizations = false;
            std::set<size_t> optimized_indices;
            std::map<size_t, std::vector<Operation>> replacements;
        };

        /**
         * @brief Create optimization plan for all gate sequences
         */
        OptimizationPlan create_optimization_plan(const std::vector<std::vector<std::pair<size_t, Operation>>> &sequences)
        {
            OptimizationPlan plan;

            for (const auto &sequence : sequences)
            {
                if (sequence.size() <= 1)
                    continue;

                size_t first_index = sequence[0].first;

                // Extract operations for optimization
                std::vector<Operation> operations;
                for (const auto &[index, operation] : sequence)
                {
                    operations.push_back(operation);
                    plan.optimized_indices.insert(index);
                }

                // Apply optimization transformations
                auto optimized_ops = optimize_gate_sequence(operations);

                // Check if optimization actually changed anything
                if (optimized_ops.size() != operations.size() || !sequences_equal(optimized_ops, operations))
                {
                    plan.has_optimizations = true;
                }

                plan.replacements[first_index] = std::move(optimized_ops);
            }

            return plan;
        }

        /**
         * @brief Apply optimizations to the circuit
         */
        void apply_optimizations(Circuit &circuit, const OptimizationPlan &plan)
        {
            const auto &original_operations = circuit.get_operations();

            // Create new circuit and copy register information
            Circuit new_circuit;
            new_circuit.add_qreg("q", circuit.get_num_qubits());
            new_circuit.add_creg("c", circuit.get_num_bits());

            for (size_t i = 0; i < original_operations.size(); ++i)
            {
                if (plan.replacements.count(i))
                {
                    // Add optimized sequence
                    for (const auto &optimized_op : plan.replacements.at(i))
                    {
                        new_circuit.add_operation(optimized_op);
                    }
                }
                else if (plan.optimized_indices.count(i) == 0)
                {
                    // Keep original operation if not part of optimized sequence
                    new_circuit.add_operation(original_operations[i]);
                }
                // Skip operations that were part of optimized sequences but not the first index
            }

            circuit = std::move(new_circuit);
        }

        /**
         * @brief Optimize a sequence of single-qubit gates
         */
        std::vector<Operation> optimize_gate_sequence(const std::vector<Operation> &operations)
        {
            auto optimized = combine_rotation_gates(operations);
            optimized = remove_canceling_gates(optimized);
            return optimized;
        }

        /**
         * @brief Check if two operation sequences are equivalent
         */
        bool sequences_equal(const std::vector<Operation> &seq1, const std::vector<Operation> &seq2) const
        {
            if (seq1.size() != seq2.size())
                return false;

            for (size_t i = 0; i < seq1.size(); ++i)
            {
                if (seq1[i].get_type() != seq2[i].get_type())
                    return false;
            }

            return true;
        }
    };
}
