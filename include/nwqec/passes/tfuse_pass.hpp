#pragma once

#include "nwqec/core/circuit.hpp"
#include "pass_template.hpp"

#include "nwqec/tableau/htab.hpp"

#include <vector>
#include <unordered_map>
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>

namespace NWQEC
{
    // PauliOp is now directly available from core/pauli_op.hpp
    class TfusePass : public Pass
    {
    private:
        mutable size_t num_qubits_ = 0;

    public:
        TfusePass() {}

        bool run(Circuit &circuit) override
        {
            num_qubits_ = circuit.get_num_qubits();

            const auto &operations = circuit.get_operations();

            // Write original circuit before optimization
            // write_operations_to_file(operations, num_qubits_, "original_circuit.qasm");

            if (!verify_pure_t_pauli_circuit(operations))
                return false;

            std::pair<std::vector<PauliOp>, std::vector<PauliOp>> pauli_rows = get_pauli_rows(operations);

            std::vector<PauliOp> t_pauli_rows = std::move(pauli_rows.first);
            std::vector<PauliOp> m_pauli_rows = std::move(pauli_rows.second);

            // // Quick check for layering benefit
            // auto greedy_layers = create_layers_greedy(t_pauli_rows);
            // auto earilest_fit_layers = create_layers(t_pauli_rows);

            // size_t average_greedy_layer_size = 0;
            // for (const auto &layer : greedy_layers)
            // {
            //     average_greedy_layer_size += layer.num_rows();
            // }
            // average_greedy_layer_size /= greedy_layers.size();

            // size_t average_earliest_fit_layer_size = 0;
            // for (const auto &layer : earilest_fit_layers)
            // {
            //     average_earliest_fit_layer_size += layer.num_rows();
            // }
            // average_earliest_fit_layer_size /= earilest_fit_layers.size();

            // std::cout << "Average greedy layer size: " << average_greedy_layer_size << ", Average earliest fit layer size: " << average_earliest_fit_layer_size << std::endl;

            // std::cout << "Greedy layers: " << greedy_layers.size() << ", Earliest fit layers: " << earilest_fit_layers.size() << std::endl;

            // exit(0);

            // std::cout << "Total M-Pauli gates: " << m_pauli_rows.size() << std::endl;

            HTab m_tab(num_qubits_);
            for (const auto &row : m_pauli_rows)
            {
                // Add measurement Pauli rows to tableau
                m_tab.add_stab(row);
            }

            auto optimized_rows = optimize(t_pauli_rows);

            std::vector<PauliOp> final_t_rows = std::move(optimized_rows.first);
            std::vector<PauliOp> final_s_rows = std::move(optimized_rows.second);

            while (true)
            {
                optimized_rows = optimize(final_t_rows);

                final_t_rows = std::move(optimized_rows.first);
                final_s_rows.insert(final_s_rows.end(), optimized_rows.second.begin(), optimized_rows.second.end());

                if (optimized_rows.second.empty()) // No more S-Pauli rows produced
                    break;
            }

            for (const auto &s_row : final_s_rows)
            {
                // Multiply S-Pauli rows with measurement tableau
                m_tab.front_multiply_pauli(s_row);
            }
            std::vector<PauliOp> m_tab_rows = m_tab.get_rows();

            update_circuit(circuit, final_t_rows, m_tab_rows);

            // Write optimized circuit after optimization
            std::vector<Operation> optimized_operations = circuit.get_operations();
            write_operations_to_file(optimized_operations, num_qubits_, "optimized_circuit.qasm", final_s_rows);

            return true;
        }

        std::string get_name() const override
        {
            return "Tfuse Pass";
        }

    private:
        bool verify_pure_t_pauli_circuit(const std::vector<Operation> &operations) const
        {
            for (const auto &op : operations)
            {
                if (!(op.get_type() == Operation::Type::T_PAULI ||
                      op.get_type() == Operation::Type::M_PAULI))
                {
                    return false;
                }
            }
            return true;
        }

        std::pair<std::vector<PauliOp>, std::vector<PauliOp>> get_pauli_rows(const std::vector<Operation> &operations) const
        {
            std::vector<PauliOp> t_pauli_rows;
            std::vector<PauliOp> m_pauli_rows;

            t_pauli_rows.reserve(operations.size());

            for (auto it = operations.rbegin(); it != operations.rend(); ++it)
            {
                if (it->get_type() == Operation::Type::T_PAULI)
                {
                    PauliOp row(num_qubits_);
                    row.from_string(it->get_pauli_string());
                    t_pauli_rows.push_back(row);
                }
                else if (it->get_type() == Operation::Type::M_PAULI)
                {
                    PauliOp row(num_qubits_);
                    row.from_string(it->get_pauli_string());
                    m_pauli_rows.push_back(row);
                }
                else
                {
                    throw std::runtime_error("Unsupported operation type in T-Pauli optimization: " + it->get_type_name());
                }
            }

            return {t_pauli_rows, m_pauli_rows};
        }

        std::vector<HTab> create_layers(const std::vector<PauliOp> &t_pauli_rows)
        {
            if (t_pauli_rows.empty())
                return {};

            std::vector<HTab> layers;

            // Add first T_PAULI operation to create the first layer
            HTab first_tableau(t_pauli_rows[0].get_num_qubits());
            first_tableau.add_stab(t_pauli_rows[0]);
            layers.push_back(first_tableau);

            // Process remaining T_PAULI operations
            for (size_t i = 1; i < t_pauli_rows.size(); i++)
            {
                bool placed = false;

                const PauliOp &pauli_row = t_pauli_rows[i];

                // Check from most recent layer (back) to oldest layer (front)
                for (int layer_idx = static_cast<int>(layers.size()) - 1; layer_idx >= 0; --layer_idx)
                {
                    bool commutes_with_layer = layers[layer_idx].commutes_with_all(pauli_row);

                    if (!commutes_with_layer)
                    {
                        // Found a non-commuting layer, place in the newer layer (layer_idx + 1)
                        if (layer_idx == static_cast<int>(layers.size()) - 1)
                        {
                            // Already at the newest layer, create a new one
                            HTab new_tableau(num_qubits_);
                            new_tableau.add_stab(pauli_row);
                            layers.push_back(new_tableau);
                        }
                        else
                        {
                            // Add to the newer layer (layer_idx + 1)
                            layers[layer_idx + 1].add_stab(pauli_row);
                        }
                        placed = true;
                        break;
                    }
                }

                if (!placed)
                {
                    // Commutes with all layers, add to the oldest layer (front)
                    layers.front().add_stab(pauli_row);
                }
            }

            return layers;
        }

        std::vector<HTab> create_layers_greedy(const std::vector<PauliOp> &t_pauli_rows)
        {
            if (t_pauli_rows.empty())
                return {};

            std::vector<HTab> layers;

            // Add first T_PAULI operation to create the first layer
            HTab first_tableau(t_pauli_rows[0].get_num_qubits());
            first_tableau.add_stab(t_pauli_rows[0]);
            layers.push_back(first_tableau);

            // Process remaining T_PAULI operations
            for (size_t i = 1; i < t_pauli_rows.size(); i++)
            {
                const PauliOp &pauli_row = t_pauli_rows[i];

                bool commutes_with_layer = layers[layers.size() - 1].commutes_with_all(pauli_row); // Check only the most recent layer

                if (!commutes_with_layer)
                {
                    // Already at the newest layer, create a new one
                    HTab new_tableau(num_qubits_);
                    new_tableau.add_stab(pauli_row);
                    layers.push_back(new_tableau);
                }
                else
                {
                    // Add to the most recent layer
                    layers.back().add_stab(pauli_row);
                }
            }

            return layers;
        }

        std::pair<std::vector<PauliOp>, std::vector<PauliOp>> optimize(std::vector<PauliOp> &t_pauli_rows)
        {
            std::vector<HTab> layers = create_layers(t_pauli_rows);

            std::vector<PauliOp> result_s_rows;

            HTab result_tab(num_qubits_);

            for (auto &layer : layers)
            {
                layer.apply_reduction();
                auto rows = layer.get_rows();

                std::vector<PauliOp> cur_layer_t_rows;

                for (const auto &row : rows)
                {
                    if (row.get_rowtype() == NWQEC::RowType::S)
                    {
                        result_tab.front_multiply_pauli(row);

                        result_s_rows.push_back(row);
                    }
                    else
                    {
                        cur_layer_t_rows.push_back(row);
                    }
                }

                for (const auto &row : cur_layer_t_rows)
                {
                    result_tab.add_stab(row);
                }
            }

            std::vector<PauliOp> result_t_rows = result_tab.get_rows();

            return {result_t_rows, result_s_rows};
        }

        void update_circuit(Circuit &circuit,
                            std::vector<PauliOp> &t_pauli_rows,
                            std::vector<PauliOp> &m_pauli_rows,
                            const std::vector<PauliOp> &s_pauli_rows = {}) const
        {

            // Rebuild circuit
            Circuit new_circuit;
            new_circuit.add_qreg("q", circuit.get_num_qubits());
            new_circuit.add_creg("c", circuit.get_num_bits());

            // Add T paulis first
            for (auto it = t_pauli_rows.rbegin(); it != t_pauli_rows.rend(); ++it)
            {
                const PauliOp &row = *it;

                assert(row.get_rowtype() != NWQEC::RowType::S && row.is_valid());
                new_circuit.add_operation(Operation(Operation::Type::T_PAULI, {}, {}, {}, row));
            }

            // Add S paulis if provided
            for (auto it = s_pauli_rows.rbegin(); it != s_pauli_rows.rend(); ++it)
            {
                const PauliOp &row = *it;

                assert(row.get_rowtype() == NWQEC::RowType::S && row.is_valid());
                new_circuit.add_operation(Operation(Operation::Type::S_PAULI, {}, {}, {}, row));
            }

            // Finally, add measurement Pauli strings
            for (const auto &row : m_pauli_rows)
            {
                Operation m_pauli_op(Operation::Type::M_PAULI, {}, {}, {}, row);
                new_circuit.add_operation(m_pauli_op);
            }

            circuit = std::move(new_circuit);
        }

        void write_operations_to_file(const std::vector<Operation> &operations,
                                      size_t num_qubits,
                                      const std::string &filename,
                                      const std::vector<PauliOp> &s_pauli_rows = {}) const
        {
            std::ofstream file(filename);
            if (!file.is_open())
                return;

            file << "OPENQASM 2.0;\n";
            file << "include \"qelib1.inc\";\n";
            file << "\n";
            file << "qreg q[" << num_qubits << "];\n";
            file << "\n";

            for (const auto &op : operations)
            {
                if (op.get_type() == Operation::Type::T_PAULI ||
                    op.get_type() == Operation::Type::S_PAULI ||
                    op.get_type() == Operation::Type::M_PAULI)
                {
                    op.print(file);
                    file << "\n";
                }
            }

            for (auto it = s_pauli_rows.rbegin(); it != s_pauli_rows.rend(); ++it)
            {
                const PauliOp &row = *it;
                if (row.is_valid() && row.get_rowtype() == NWQEC::RowType::S)
                {
                    file << "s_pauli " << row.to_string() << ";\n";
                }
            }

            file.close();
        }
    };

} // namespace NWQEC