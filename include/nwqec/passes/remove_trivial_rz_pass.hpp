#pragma once

#include "nwqec/core/circuit.hpp"
#include "pass_template.hpp"
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <set>
#include <optional>
#include <algorithm>

namespace NWQEC
{
    /**
     * @brief Pass to optimize RZ gates
     *
     * This pass combines consecutive RZ gates on the same qubit and removes RZ(0) gates.
     * Also tracks RZ angles for analysis purposes.
     */
    class RemoveTrivialRzPass : public Pass
    {
    public:
        RemoveTrivialRzPass() {}

    private:
        const double TOLERANCE = 1e-4;

        /**
         * @brief Convert angle to string with specified number of significant digits
         * @param angle The angle to convert
         * @param sig_digits Number of significant digits (default: 4)
         */
        std::string angle_to_string(double angle, int sig_digits = 4) const
        {
            if (angle == 0.0)
            {
                std::string result = "0.";
                result.append(std::max(0, sig_digits - 1), '0');
                return result;
            }

            // Get the order of magnitude
            double abs_angle = std::abs(angle);
            int order = static_cast<int>(std::floor(std::log10(abs_angle)));

            // Scale to get specified significant digits
            double scaled = abs_angle / std::pow(10, order - (sig_digits - 1));
            int scaled_int = static_cast<int>(std::round(scaled));

            // Convert back to the original scale
            double rounded = scaled_int * std::pow(10, order - (sig_digits - 1));
            if (angle < 0)
                rounded = -rounded;

            std::stringstream ss;
            ss << std::fixed << std::setprecision(std::max(0, (sig_digits - 1) - order)) << rounded;
            return ss.str();
        }

    public:
        std::string get_name() const override
        {
            return "Remove Trivial RZ Pass";
        }

        bool run(Circuit &circuit) override
        {
            bool circuit_modified = false;

            // Create new circuit and copy register information
            Circuit new_circuit;
            new_circuit.add_qreg("q", circuit.get_num_qubits());
            new_circuit.add_creg("c", circuit.get_num_bits());

            // Track distinct RZ angles for grouping
            std::vector<std::string> distinct_angle_strings;

            // Process each operation
            for (const auto &operation : circuit.get_operations())
            {
                if (operation.get_type() != Operation::Type::RZ)
                {
                    new_circuit.add_operation(operation);
                    continue;
                }

                circuit_modified = true;
                process_rz_gate(operation, new_circuit, distinct_angle_strings);
            }

            // Replace circuit if modifications were made
            if (circuit_modified)
            {
                circuit = std::move(new_circuit);
            }

            return circuit_modified;
        }

    private:
        /**
         * @brief Process an RZ gate and add optimized operations to the circuit
         */
        void process_rz_gate(const Operation &operation, Circuit &new_circuit,
                             std::vector<std::string> &distinct_angle_strings)
        {
            auto params = operation.get_parameters();
            auto qubits = operation.get_qubits();

            if (params.empty())
            {
                throw std::runtime_error("RZ gate has no parameters");
            }

            double angle = normalize_angle(params[0]);

            // Check for identity operations
            if (is_identity_angle(angle))
                return;

            // Try to replace with standard gates
            if (auto standard_gate = get_standard_gate_replacement(angle, qubits))
            {
                new_circuit.add_operation(*standard_gate);
                return;
            }

            // Try exact decomposition into Z, S, T gates
            if (auto decomposed_ops = try_exact_decomposition(angle, qubits); !decomposed_ops.empty())
            {
                for (const auto &op : decomposed_ops)
                {
                    new_circuit.add_operation(op);
                }
                return;
            }

            // Keep as RZ gate with angle grouping
            add_grouped_rz_gate(angle, qubits, new_circuit, distinct_angle_strings);
        }

        /**
         * @brief Normalize angle to [0, 2π)
         */
        double normalize_angle(double angle) const
        {
            angle = std::fmod(angle, 2 * M_PI);
            if (angle < 0)
                angle += 2 * M_PI;
            return angle;
        }

        /**
         * @brief Check if angle represents identity operation
         */
        bool is_identity_angle(double angle) const
        {
            return std::abs(angle) < TOLERANCE || std::abs(angle - 2 * M_PI) < TOLERANCE;
        }

        /**
         * @brief Get standard gate replacement for common angles
         */
        std::optional<Operation> get_standard_gate_replacement(double angle, const std::vector<size_t> &qubits) const
        {
            if (std::abs(angle - M_PI) < TOLERANCE)
                return Operation(Operation::Type::Z, qubits);

            if (std::abs(angle - M_PI_2) < TOLERANCE)
                return Operation(Operation::Type::S, qubits);

            if (std::abs(angle - 3 * M_PI_2) < TOLERANCE)
                return Operation(Operation::Type::SDG, qubits);

            if (std::abs(angle - M_PI_4) < TOLERANCE)
                return Operation(Operation::Type::T, qubits);

            if (std::abs(angle - 7 * M_PI_4) < TOLERANCE)
                return Operation(Operation::Type::TDG, qubits);

            return std::nullopt;
        }

        /**
         * @brief Try to decompose angle into exact Z, S, T gate sequence
         */
        std::vector<Operation> try_exact_decomposition(double angle, const std::vector<size_t> &qubits) const
        {
            std::vector<Operation> decomposed;
            double remaining = angle;

            // Add Z gates (π)
            if (remaining >= M_PI - TOLERANCE)
            {
                decomposed.push_back(Operation(Operation::Type::Z, qubits));
                remaining -= M_PI;
            }

            // Add S gates (π/2)
            if (remaining >= M_PI_2 - TOLERANCE)
            {
                decomposed.push_back(Operation(Operation::Type::S, qubits));
                remaining -= M_PI_2;
            }

            // Add T gates (π/4)
            if (remaining >= M_PI_4 - TOLERANCE)
            {
                decomposed.push_back(Operation(Operation::Type::T, qubits));
                remaining -= M_PI_4;
            }

            // Return decomposition only if exact (within tolerance)
            if (std::abs(remaining) < TOLERANCE)
                return decomposed;

            return {};
        }

        /**
         * @brief Add RZ gate with angle grouping for synthesis
         */
        void add_grouped_rz_gate(double angle, const std::vector<size_t> &qubits,
                                 Circuit &new_circuit, std::vector<std::string> &distinct_angle_strings)
        {
            std::string angle_str = angle_to_string(angle);

            auto it = std::find(distinct_angle_strings.begin(), distinct_angle_strings.end(), angle_str);
            size_t angle_index;

            if (it == distinct_angle_strings.end())
            {
                angle_index = distinct_angle_strings.size();
                distinct_angle_strings.push_back(angle_str);
                new_circuit.distinct_rz_angles.push_back(angle);
            }
            else
            {
                angle_index = std::distance(distinct_angle_strings.begin(), it);
            }

            new_circuit.add_operation(Operation(Operation::Type::RZ, qubits, {angle}));
            new_circuit.rz_angle_map[new_circuit.get_operations().size() - 1] = angle_index;
        }
    };

} // namespace NWQEC