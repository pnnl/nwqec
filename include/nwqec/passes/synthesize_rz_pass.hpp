#ifndef NWQEC_WITH_GRIDSYNTH_CPP
#define NWQEC_WITH_GRIDSYNTH_CPP 0
#endif

#pragma once

#include "nwqec/core/circuit.hpp"
#if NWQEC_WITH_GRIDSYNTH_CPP
#include "nwqec/gridsynth/gridsynth.hpp"
#endif
#include "nwqec/core/constants.hpp"

#include "pass_template.hpp"
#include <vector>
#include <cmath>
#include <chrono>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <map>
#include <string>
#include <sstream>
#include <iostream>

namespace NWQEC
{
    /**
     * @brief Pass to optimize RZ gates
     *
     * This pass combines consecutive RZ gates on the same qubit and removes RZ(0) gates.
     */
    class SynthesizeRzPass : public Pass
    {
    private:
        const double synthesis_error_ = NWQEC::DEFAULT_EPSILON_MULTIPLIER; // Default synthesis error tolerance multiplier
        double epsilon_override_ = -1.0;                                   // If >=0, use this epsilon for all angles

    public:
        SynthesizeRzPass() = default;
        explicit SynthesizeRzPass(double epsilon_override) : epsilon_override_(epsilon_override) {}

        std::string get_name() const override
        {
            return "Synthesize RZ Pass";
        }

        bool run(Circuit &circuit) override
        {

            bool circuit_modified = false;

            // Create new circuit and copy register information
            Circuit new_circuit;
            new_circuit.add_qreg("q", circuit.get_num_qubits());
            new_circuit.add_creg("c", circuit.get_num_bits());

            // Ensure RZ angle grouping is done
            ensure_rz_angle_grouping(circuit);

            // Pre-synthesize all distinct RZ angles
            auto pre_synthesized_gates = synthesize_all_angles(circuit.distinct_rz_angles);

            // Process each operation
            const auto &operations = circuit.get_operations();
            for (size_t i = 0; i < operations.size(); ++i)
            {
                const auto &operation = operations[i];

                if (operation.get_type() != Operation::Type::RZ)
                {
                    new_circuit.add_operation(operation);
                    continue;
                }

                circuit_modified = true;
                synthesize_rz_operation(operation, i, circuit, pre_synthesized_gates, new_circuit);
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
         * @brief Ensure RZ angle grouping is performed if not already done
         */
        void ensure_rz_angle_grouping(Circuit &circuit)
        {
            if (!circuit.distinct_rz_angles.empty())
                return;

            std::vector<double> distinct_angles;
            std::vector<std::string> distinct_angle_strings;
            std::map<size_t, size_t> rz_angle_map;

            const auto &operations = circuit.get_operations();
            for (size_t i = 0; i < operations.size(); ++i)
            {
                const auto &operation = operations[i];
                if (operation.get_type() == Operation::Type::RZ)
                {
                    auto params = operation.get_parameters();
                    if (!params.empty())
                    {
                        double angle = params[0];
                        std::string angle_str = angle_to_string(angle);

                        auto it = std::find(distinct_angle_strings.begin(), distinct_angle_strings.end(), angle_str);
                        size_t angle_index;

                        if (it == distinct_angle_strings.end())
                        {
                            angle_index = distinct_angle_strings.size();
                            distinct_angle_strings.push_back(angle_str);
                            distinct_angles.push_back(angle);
                        }
                        else
                        {
                            angle_index = std::distance(distinct_angle_strings.begin(), it);
                        }

                        rz_angle_map[i] = angle_index;
                    }
                }
            }

            circuit.distinct_rz_angles = std::move(distinct_angles);
            circuit.rz_angle_map = std::move(rz_angle_map);
        }

        /**
         * @brief Synthesize all distinct RZ angles
         */
        std::vector<std::string> synthesize_all_angles(const std::vector<double> &angles)
        {
            std::vector<std::string> synthesized_gates;
            synthesized_gates.reserve(angles.size());

            for (const auto &angle : angles)
            {
                synthesized_gates.push_back(synthesize_angle(angle));
            }

            return synthesized_gates;
        }

        /**
         * @brief Synthesize a single RZ operation
         */
        void synthesize_rz_operation(const Operation &operation, size_t operation_index,
                                     const Circuit &circuit, const std::vector<std::string> &pre_synthesized_gates,
                                     Circuit &new_circuit)
        {
            auto qubits = operation.get_qubits();

            auto angle_map_it = circuit.rz_angle_map.find(operation_index);
            if (angle_map_it == circuit.rz_angle_map.end() || pre_synthesized_gates.empty())
            {
                throw std::runtime_error("RZ gate found without corresponding pre-synthesized gate");
            }

            const std::string &gate_sequence = pre_synthesized_gates[angle_map_it->second];
            add_gate_sequence_to_circuit(gate_sequence, qubits, new_circuit);
        }

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

        std::string synthesize_angle(double angle)
        {
            try
            {
                std::string theta = std::to_string(angle);
                // Choose epsilon: explicit override or default multiplier rule
                double epsilon_abs = (epsilon_override_ >= 0.0) ? epsilon_override_ : synthesis_error_ * std::abs(angle);
                // Call the gridsynth function with epsilon in scientific notation to avoid truncation
                std::ostringstream eps_ss;
                eps_ss.setf(std::ios::scientific);
                eps_ss << std::setprecision(16) << epsilon_abs;
                return gridsynth::gridsynth_gates(theta, eps_ss.str());
            }
            catch (const std::exception &e)
            {
                std::cerr << "Warning: Failed to synthesize RZ angle " << angle
                          << ": " << e.what() << std::endl;
                std::cerr << "Falling back to identity (no gates)" << std::endl;
                return ""; // Return empty string as fallback
            }
        }

        void add_gate_sequence_to_circuit(const std::string &gate_sequence, const std::vector<size_t> &qubits, Circuit &circuit)
        {
            for (const char &gate : gate_sequence)
            {
                switch (gate)
                {
                case 'X':
                    circuit.add_operation(Operation(Operation::Type::X, qubits));
                    break;
                case 'Y':
                    circuit.add_operation(Operation(Operation::Type::Y, qubits));
                    break;
                case 'Z':
                    circuit.add_operation(Operation(Operation::Type::Z, qubits));
                    break;
                case 'H':
                    circuit.add_operation(Operation(Operation::Type::H, qubits));
                    break;
                case 'S':
                    circuit.add_operation(Operation(Operation::Type::S, qubits));
                    break;
                case 'T':
                    circuit.add_operation(Operation(Operation::Type::T, qubits));
                    break;
                case 'W':
                    // Skip global phase
                    break;
                default:
                    std::cerr << "Unknown gate type: " << gate << std::endl;
                }
            }
        }
    };

} // namespace NWQEC
