#pragma once

#include "nwqec/core/circuit.hpp"
#include "pass_template.hpp"
#include <set>
#include <string>
#include <cmath>
#include <cassert>

namespace NWQEC
{
    /**
     * @brief Pass to decompose gates into basis gates
     *
     * This pass decomposes complex gates into a sequence of basis gates.
     */
    class DecomposePass : public Pass
    {

    public:
        /**
         * @brief Construct a decomposition pass with specified basis gates
         *
         * If no basis gates are provided, defaults to {CX, U, ID, RZ, SX}.
         *
         * @param basis_gates Set of gate types to use as basis
         */
        // DecompositionPass() {}

        DecomposePass(bool keep_ccx = false) : keep_ccx_(keep_ccx) {}

        std::string get_name() const override
        {
            return "Decompose Gates Pass";
        }

        bool run(Circuit &circuit) override
        {

            bool circuit_modified = false;

            // Create new circuit and copy register information
            Circuit new_circuit;
            new_circuit.add_qreg("q", circuit.get_num_qubits());
            new_circuit.add_creg("c", circuit.get_num_bits());

            // Process each operation
            for (const auto &operation : circuit.get_operations())
            {
                if (should_keep_operation(operation, circuit))
                {
                    new_circuit.add_operation(operation);
                }
                else if (auto decomposed_ops = decompose_gate(operation); !decomposed_ops.empty())
                {
                    for (const auto &decomposed_op : decomposed_ops)
                    {
                        new_circuit.add_operation(decomposed_op);
                    }
                    circuit_modified = true;
                }
                else
                {
                    // Keep original operation if decomposition not available
                    new_circuit.add_operation(operation);
                }
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
         * @brief Check if an operation should be kept as-is without decomposition
         */
        bool should_keep_operation(const Operation &operation, const Circuit &circuit) const
        {
            return circuit.is_clifford_t_operation(operation.get_type()) ||
                   (keep_ccx_ && operation.get_type() == Operation::Type::CCX);
        }

        /**
         * @brief Decompose a gate into a sequence of basis gates
         *
         * @param op The operation to decompose
         * @return Vector of decomposed operations
         */
        std::vector<Operation> decompose_gate(const Operation &op)
        {
            std::vector<Operation> result;

            // Get the operation parameters and qubits
            auto type = op.get_type();
            auto qubits = op.get_qubits();
            auto params = op.get_parameters();

            switch (type)
            {

            case Operation::Type::RX:
            {
                result.push_back(Operation(Operation::Type::H, {qubits[0]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[0]}, {params[0]}));
                result.push_back(Operation(Operation::Type::H, {qubits[0]}, {}));
                break;
            }
            case Operation::Type::RY:
            {
                result.push_back(Operation(Operation::Type::SDG, {qubits[0]}, {}));
                result.push_back(Operation(Operation::Type::H, {qubits[0]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[0]}, {params[0]}));
                result.push_back(Operation(Operation::Type::H, {qubits[0]}, {}));
                result.push_back(Operation(Operation::Type::S, {qubits[0]}, {}));
                break;
            }
            case Operation::Type::P:
            {
                result.push_back(Operation(Operation::Type::RZ, {qubits[0]}, {params[0]}));
                break;
            }
            case Operation::Type::U:
            case Operation::Type::U3:
            {
                result.push_back(Operation(Operation::Type::RZ, {qubits[0]}, {params[2]})); // λ
                result.push_back(Operation(Operation::Type::SX, {qubits[0]}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[0]}, {params[0] + M_PI})); // θ
                result.push_back(Operation(Operation::Type::SX, {qubits[0]}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[0]}, {params[1] + 3 * M_PI})); // φ

                break;
            }
            case Operation::Type::U1: // u1(θ) = rz(θ)
            {
                result.push_back(Operation(Operation::Type::RZ, {qubits[0]}, {params[0]}));
                break;
            }
            case Operation::Type::U2: // u2(θ, φ) = u(pi/2, φ, θ)
            {
                result.push_back(Operation(Operation::Type::RZ, {qubits[0]}, {params[1]})); // λ
                result.push_back(Operation(Operation::Type::SX, {qubits[0]}));
                result.push_back(Operation(Operation::Type::S, {qubits[0]}));
                result.push_back(Operation(Operation::Type::Z, {qubits[0]}));
                result.push_back(Operation(Operation::Type::SX, {qubits[0]}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[0]}, {params[0] - 3 * M_PI})); // φ
                break;
            }
            case Operation::Type::CY:
            {
                result.push_back(Operation(Operation::Type::SDG, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::S, {qubits[1]}, {}));
                break;
            }
            case Operation::Type::CZ:
            {
                result.push_back(Operation(Operation::Type::H, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::H, {qubits[1]}, {}));
                break;
            }
            case Operation::Type::CH:
            {
                result.push_back(Operation(Operation::Type::S, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::H, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::T, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::TDG, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::H, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::SDG, {qubits[1]}, {}));
                break;
            }
            case Operation::Type::CS:
            {
                result.push_back(Operation(Operation::Type::S, {qubits[0]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::SDG, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::S, {qubits[1]}, {}));
                break;
            }
            case Operation::Type::CSDG:
            {
                result.push_back(Operation(Operation::Type::TDG, {qubits[0]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::T, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::TDG, {qubits[1]}, {}));
                break;
            }
            case Operation::Type::CT:
            {
                result.push_back(Operation(Operation::Type::RZ, {qubits[0]}, {M_PI / 8}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {-M_PI / 8}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {M_PI / 8}));

                break;
            }
            case Operation::Type::CTDG:
            {
                result.push_back(Operation(Operation::Type::RZ, {qubits[0]}, {-M_PI / 8}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {M_PI / 8}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {-M_PI / 8}));
                break;
            }
            case Operation::Type::CSX:
            {

                result.push_back(Operation(Operation::Type::T, {qubits[0]}, {}));
                result.push_back(Operation(Operation::Type::H, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::TDG, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::T, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::H, {qubits[1]}, {}));
                break;
            }
            case Operation::Type::CRX:
            {
                result.push_back(Operation(Operation::Type::H, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {params[0] / 2}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {-params[0] / 2}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::H, {qubits[1]}, {}));
                break;
            }
            case Operation::Type::CRY:
            {

                result.push_back(Operation(Operation::Type::SX, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {params[0] / 2}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {-params[0] / 2}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::SXDG, {qubits[1]}, {}));
                break;
            }
            case Operation::Type::CRZ:
            {
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {params[0] / 2}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {-params[0] / 2}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                break;
            }
            case Operation::Type::CP:
            case Operation::Type::CU1:
            {
                result.push_back(Operation(Operation::Type::RZ, {qubits[0]}, {params[0] / 2}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {-params[0] / 2}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {params[0] / 2}));
                break;
            }
            case Operation::Type::CU:
            case Operation::Type::CU3:
            {

                result.push_back(Operation(Operation::Type::RZ, {qubits[0]}, {params[3] + params[2] / 2 + params[1] / 2}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {params[2] / 2 - params[1] / 2}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {-params[2] / 2 - params[1] / 2}));
                result.push_back(Operation(Operation::Type::SX, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {M_PI - params[0] / 2}));
                result.push_back(Operation(Operation::Type::SX, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::Z, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::SX, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {params[0] / 2 + M_PI}));
                result.push_back(Operation(Operation::Type::SX, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {params[1] + 3 * M_PI}));
                break;
            }

            case Operation::Type::RXX:
            {

                result.push_back(Operation(Operation::Type::H, {qubits[0]}, {}));
                result.push_back(Operation(Operation::Type::H, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {params[0]}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::H, {qubits[0]}, {}));
                result.push_back(Operation(Operation::Type::H, {qubits[1]}, {}));
                break;
            }
            case Operation::Type::RYY:
            {

                result.push_back(Operation(Operation::Type::SX, {qubits[0]}, {}));
                result.push_back(Operation(Operation::Type::SX, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {params[0]}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::SXDG, {qubits[0]}, {}));
                result.push_back(Operation(Operation::Type::SXDG, {qubits[1]}, {}));
                break;
            }
            case Operation::Type::RZZ:
            {
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::RZ, {qubits[1]}, {params[0]}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                break;
            }
            case Operation::Type::SWAP:
            {
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}));
                result.push_back(Operation(Operation::Type::CX, {qubits[1], qubits[0]}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}));
                break;
            }
            case Operation::Type::CCX:
            {
                result.push_back(Operation(Operation::Type::H, {qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[1], qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::TDG, {qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::T, {qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[1], qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::T, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::TDG, {qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::T, {qubits[0]}, {}));
                result.push_back(Operation(Operation::Type::TDG, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::T, {qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::H, {qubits[2]}, {}));

                break;
            }
            case Operation::Type::CSWAP:
            {
                result.push_back(Operation(Operation::Type::CX, {qubits[2], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::H, {qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[1], qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::TDG, {qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::T, {qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[1], qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::T, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::TDG, {qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::T, {qubits[0]}, {}));
                result.push_back(Operation(Operation::Type::TDG, {qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[1]}, {}));
                result.push_back(Operation(Operation::Type::T, {qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::H, {qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[2], qubits[1]}, {}));
                break;
            }
            case Operation::Type::RCCX:
            {
                result.push_back(Operation(Operation::Type::H, {qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::T, {qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[1], qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::TDG, {qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[0], qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::T, {qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::CX, {qubits[1], qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::TDG, {qubits[2]}, {}));
                result.push_back(Operation(Operation::Type::H, {qubits[2]}, {}));
                break;
            }
            default:
                // For unsupported operations, return empty vector
                return {};
            }

            return result;
        }

    private:
        bool keep_ccx_ = false;
    };

} // namespace NWQEC