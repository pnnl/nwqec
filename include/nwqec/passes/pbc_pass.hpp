#pragma once

#include "pass_template.hpp"
#include "nwqec/tableau/vtab.hpp"
#include "nwqec/core/pauli_op.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cassert>
#include <algorithm>

namespace NWQEC
{
    class PbcPass : public Pass
    {
    private:
        bool keep_cx;

    public:
        PbcPass(bool keep_cx = false) : keep_cx(keep_cx) {}

        bool run(Circuit &circuit) override
        {
            std::vector<Operation> operations = circuit.get_operations();
            
            // Check if circuit already contains PBC operations
            bool has_pbc_operations = false;
            for (const auto &op : operations) {
                if (op.get_type() == Operation::Type::T_PAULI ||
                    op.get_type() == Operation::Type::S_PAULI ||
                    op.get_type() == Operation::Type::M_PAULI ||
                    op.get_type() == Operation::Type::Z_PAULI) {
                    has_pbc_operations = true;
                    break;
                }
            }
            
            // If circuit already contains PBC operations, skip processing
            if (has_pbc_operations) {
                return false; // No modification needed
            }
            
            size_t n_qubits = circuit.get_num_qubits();
            size_t n_gate_stabs = 0;
            std::vector<Operation::Type> gate_types;
            std::vector<size_t> qubit_a_list;
            std::vector<size_t> qubit_b_list;
            std::vector<uint8_t> phase_bits;
            std::vector<PauliOp> pbc_stabs;
            std::vector<bool> is_t_stab;

            // Process operations in reverse order - simpler approach without sequence optimization for now
            for (auto it = operations.rbegin(); it != operations.rend(); ++it)
            {
                if (it->get_type() == Operation::Type::MEASURE ||
                    it->get_type() == Operation::Type::RESET ||
                    it->get_type() == Operation::Type::BARRIER)
                {
                    continue; // Skip measurement and reset operations
                }

                if (it->get_type() == Operation::Type::CCX)
                {
                    // Expand CCX to 7 individual T_PAULI operations
                    auto qubits = it->get_qubits();
                    auto ccx_rows = PauliOp::create_ccx_ops(qubits[0], qubits[1], qubits[2], n_qubits);

                    // Add each CCX stabilizer as a separate gate entry
                    for (const auto &stab : ccx_rows)
                    {
                        gate_types.push_back(Operation::Type::T_PAULI);
                        qubit_a_list.push_back(0); // Dummy values for T_PAULI
                        qubit_b_list.push_back(SIZE_MAX);
                        phase_bits.push_back(0);
                        pbc_stabs.push_back(stab);
                        n_gate_stabs++;

                        is_t_stab.push_back(true);
                    }
                }
                else if (keep_cx && it->get_type() == Operation::Type::CX)
                {
                    // Add CX gate
                    auto qubits = it->get_qubits();

                    // Add Sdg q[0]
                    gate_types.push_back(Operation::Type::SDG);
                    qubit_a_list.push_back(qubits[0]);
                    qubit_b_list.push_back(SIZE_MAX);
                    phase_bits.push_back(0);

                    // Add Sxdg q[1]
                    gate_types.push_back(Operation::Type::SXDG);
                    qubit_a_list.push_back(qubits[1]);
                    qubit_b_list.push_back(SIZE_MAX);
                    phase_bits.push_back(0);

                    // Add S_PAULI Z_q[0] X_q[1]
                    gate_types.push_back(Operation::Type::S_PAULI);
                    qubit_a_list.push_back(0); // Dummy values for S_PAULI
                    qubit_b_list.push_back(SIZE_MAX);
                    phase_bits.push_back(0);

                    PauliOp stab(n_qubits);
                    stab.set_r(false); // S_PAULI stabilizer
                    stab.add_z(qubits[0]);
                    stab.add_x(qubits[1]);
                    pbc_stabs.push_back(stab);
                    n_gate_stabs++;
                    is_t_stab.push_back(false);
                }
                else
                {
                    // Add regular gate
                    auto qubits = it->get_qubits();
                    gate_types.push_back(it->get_type());
                    qubit_a_list.push_back(qubits[0]);
                    qubit_b_list.push_back(qubits.size() > 1 ? qubits[1] : SIZE_MAX);
                    uint8_t phase = (it->get_type() == Operation::Type::T) ? 0 : 1;
                    phase_bits.push_back(phase);

                    if (it->get_type() == Operation::Type::T || it->get_type() == Operation::Type::TDG)
                    {
                        n_gate_stabs++;
                        is_t_stab.push_back(true);
                    }
                }
            }

            assert(gate_types.size() == qubit_a_list.size());
            assert(gate_types.size() == qubit_b_list.size());
            assert(gate_types.size() == phase_bits.size());
            // Create tableau with expanded gate list
            VTab tableau(n_qubits, n_gate_stabs, gate_types, qubit_a_list, qubit_b_list, phase_bits, pbc_stabs);
            std::vector<PauliOp> stabilizers = tableau.get_paili_ops();

            update_circuit(stabilizers, circuit, is_t_stab);
            return true;
        }

        std::string get_name() const override
        {
            return "PBC Pass";
        }

    private:
        void update_circuit(std::vector<PauliOp> &stabilizers, Circuit &circuit, std::vector<bool> &is_t_stab)
        {
            size_t n_qubits = circuit.get_num_qubits();
            // First n stabilizers are measurement Pauli strings
            std::vector<PauliOp> measurement_stabilizers(stabilizers.begin(), stabilizers.begin() + n_qubits);

            // Remaining stabilizers are rotation Pauli strings (reverse order)
            std::vector<PauliOp> rotation_stabilizers(stabilizers.begin() + n_qubits, stabilizers.end());

            // Build new circuit from stabilizers
            Circuit new_circuit;
            new_circuit.add_qreg("q", n_qubits);

            // Add rotation operations from back to front
            std::reverse(is_t_stab.begin(), is_t_stab.end());
            size_t stab_idx = 0;

            assert(rotation_stabilizers.size() == is_t_stab.size());
            for (auto it = rotation_stabilizers.rbegin(); it != rotation_stabilizers.rend(); ++it)
            {
                const PauliOp &pauli_op = *it;
                
                if (is_t_stab[stab_idx++])
                {
                    // Add T_PAULI operation
                    new_circuit.add_operation(Operation(Operation::Type::T_PAULI, {}, {}, {}, pauli_op));
                }
                else
                {
                    // Add S_PAULI operation
                    new_circuit.add_operation(Operation(Operation::Type::S_PAULI, {}, {}, {}, pauli_op));
                }
            }

            // Add measurement operations
            for (const auto &pauli_op : measurement_stabilizers)
            {
                new_circuit.add_operation(Operation(Operation::Type::M_PAULI, {}, {}, {}, pauli_op));
            }

            circuit = std::move(new_circuit);
        }
    };

} // namespace NWQEC