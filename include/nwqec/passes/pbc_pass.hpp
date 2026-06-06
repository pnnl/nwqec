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
            const std::vector<Operation> &operations = circuit.get_operations();

            std::vector<const Operation *> pbc_operations;
            pbc_operations.reserve(operations.size());
            for (const auto &op : operations)
            {
                if (op.get_type() == Operation::Type::T_PAULI ||
                    op.get_type() == Operation::Type::S_PAULI ||
                    op.get_type() == Operation::Type::M_PAULI ||
                    op.get_type() == Operation::Type::Z_PAULI)
                {
                    return false; // Already a PBC circuit
                }

                if (!pbc_operations.empty() &&
                    is_self_inverse_gate(op.get_type()) &&
                    same_gate_target(*pbc_operations.back(), op))
                {
                    pbc_operations.pop_back();
                    continue;
                }
                pbc_operations.push_back(&op);
            }

            size_t n_qubits = circuit.get_num_qubits();
            size_t n_gate_stabs = 0;
            std::vector<bool> is_t_stab;
            is_t_stab.reserve(pbc_operations.size());

            for (const Operation *op_ptr : pbc_operations)
            {
                const Operation &op = *op_ptr;
                if (op.get_type() == Operation::Type::CCX)
                {
                    n_gate_stabs += 7;
                    for (size_t i = 0; i < 7; ++i)
                        is_t_stab.push_back(true);
                }
                else if (keep_cx && op.get_type() == Operation::Type::CX)
                {
                    ++n_gate_stabs;
                    is_t_stab.push_back(false);
                }
                else if (op.get_type() == Operation::Type::T || op.get_type() == Operation::Type::TDG)
                {
                    ++n_gate_stabs;
                    is_t_stab.push_back(true);
                }
            }

            VTab tableau(n_qubits, n_gate_stabs);

            for (auto it = pbc_operations.rbegin(); it != pbc_operations.rend(); ++it)
            {
                const Operation &op = **it;
                if (op.get_type() == Operation::Type::MEASURE ||
                    op.get_type() == Operation::Type::RESET ||
                    op.get_type() == Operation::Type::BARRIER)
                {
                    continue; // Skip measurement and reset operations
                }

                if (op.get_type() == Operation::Type::CCX)
                {
                    const auto &qubits = op.get_qubits();
                    tableau.add_ccx_stabs(qubits[0], qubits[1], qubits[2]);
                }
                else if (keep_cx && op.get_type() == Operation::Type::CX)
                {
                    // Add CX gate
                    const auto &qubits = op.get_qubits();

                    // Add Sdg q[0]
                    tableau.apply_clifford_gate(Operation::Type::SDG, qubits[0]);

                    // Add Sxdg q[1]
                    tableau.apply_clifford_gate(Operation::Type::SXDG, qubits[1]);

                    // Add S_PAULI Z_q[0] X_q[1]
                    PauliOp stab(n_qubits);
                    stab.set_r(false); // S_PAULI stabilizer
                    stab.add_z(qubits[0]);
                    stab.add_x(qubits[1]);
                    tableau.add_stab(stab);
                }
                else
                {
                    // Add regular gate
                    const auto &qubits = op.get_qubits();
                    if (op.get_type() == Operation::Type::T || op.get_type() == Operation::Type::TDG)
                    {
                        uint8_t phase = (op.get_type() == Operation::Type::T) ? 0 : 1;
                        tableau.add_t_stab(qubits[0], phase);
                    }
                    else
                    {
                        tableau.apply_clifford_gate(op.get_type(), qubits[0], qubits.size() > 1 ? qubits[1] : SIZE_MAX);
                    }
                }
            }
            std::vector<PauliOp> stabilizers = tableau.get_paili_ops();

            update_circuit(stabilizers, circuit, is_t_stab);
            return true;
        }

        std::string get_name() const override
        {
            return "PBC Pass";
        }

    private:
        bool is_self_inverse_gate(Operation::Type type) const
        {
            return type == Operation::Type::H ||
                   type == Operation::Type::X ||
                   type == Operation::Type::Y ||
                   type == Operation::Type::Z ||
                   type == Operation::Type::CX;
        }

        bool same_gate_target(const Operation &lhs, const Operation &rhs) const
        {
            return lhs.get_type() == rhs.get_type() &&
                   lhs.get_qubits() == rhs.get_qubits();
        }

        void update_circuit(std::vector<PauliOp> &stabilizers, Circuit &circuit, std::vector<bool> &is_t_stab)
        {
            size_t n_qubits = circuit.get_num_qubits();

            // Build new circuit from stabilizers
            Circuit new_circuit;
            new_circuit.add_qreg("q", n_qubits);

            // Add rotation operations from back to front
            std::reverse(is_t_stab.begin(), is_t_stab.end());
            size_t stab_idx = 0;

            assert(stabilizers.size() == n_qubits + is_t_stab.size());
            for (size_t pos = stabilizers.size(); pos-- > n_qubits;)
            {
                const PauliOp &pauli_op = stabilizers[pos];
                
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
            for (size_t i = 0; i < n_qubits; ++i)
            {
                const PauliOp &pauli_op = stabilizers[i];
                new_circuit.add_operation(Operation(Operation::Type::M_PAULI, {}, {}, {}, pauli_op));
            }

            circuit = std::move(new_circuit);
        }
    };

} // namespace NWQEC
