#pragma once

#include "nwqec/core/circuit.hpp"
#include "pass_template.hpp"
#include <vector>

namespace NWQEC
{
    /**
     * @brief Pass to remove all Pauli gates (X, Y, Z) from the circuit
     *
     * This pass removes all Pauli X, Y, and Z gates from the circuit.
     * Useful for analysis purposes where Pauli gates are considered "free"
     * operations that don't contribute to circuit complexity.
     */
    class RemovePauliPass : public Pass
    {
    public:
        RemovePauliPass() = default;

        std::string get_name() const override
        {
            return "Remove Pauli Gates Pass";
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
                // Skip Pauli gates (X, Y, Z)
                if (operation.get_type() == Operation::Type::X ||
                    operation.get_type() == Operation::Type::Y ||
                    operation.get_type() == Operation::Type::Z)
                {
                    circuit_modified = true;
                    continue; // Skip this operation
                }

                // Keep all other operations
                new_circuit.add_operation(operation);
            }

            // Replace circuit if modifications were made
            if (circuit_modified)
            {
                circuit = std::move(new_circuit);
            }

            return circuit_modified;
        }
    };

} // namespace NWQEC