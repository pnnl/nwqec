#pragma once

#include <memory>
#include "nwqec/core/circuit.hpp"
#include "nwqec/parser/ast.hpp"
#include <unordered_map>
#include <algorithm>
#include <cctype>

namespace NWQEC
{

    /**
     * CircuitBuilder converts a parsed AST into a flattened circuit representation
     * with all user-defined gates expanded to fundamental gates
     */
    class ASTCircuitConverter
    {
    private:
        struct gate_defination
        {
            std::vector<std::string> params;
            std::vector<std::string> qubits;
            std::vector<std::unique_ptr<Stmt>> body;
        };

        // Maps gate names to their definitions
        std::unordered_map<std::string, gate_defination> gate_definitions;
        Circuit circuit;

        // Stack of currently active parameter bindings for gate expansion
        std::vector<std::unordered_map<std::string, std::vector<size_t>>> qubit_binding_stack;

        // Push/pop qubit bindings during gate expansion
        void push_qubit_bindings(const std::vector<std::string> &formal_params, const std::vector<std::vector<size_t>> &actual_indices)
        {
            std::unordered_map<std::string, std::vector<size_t>> bindings;
            for (size_t i = 0; i < formal_params.size() && i < actual_indices.size(); ++i)
            {
                bindings[formal_params[i]] = actual_indices[i];
            }
            qubit_binding_stack.push_back(std::move(bindings));
        }

        void pop_qubit_bindings()
        {
            if (!qubit_binding_stack.empty())
            {
                qubit_binding_stack.pop_back();
            }
        }

        // Get the actual qubit index for a variable, considering current bindings
        std::vector<size_t> resolve_qubit_index(const std::string &reg_name, int local_index)
        {
            // Check if it's a bound parameter (from innermost scope outward)
            for (auto it = qubit_binding_stack.rbegin(); it != qubit_binding_stack.rend(); ++it)
            {
                auto binding = it->find(reg_name);
                if (binding != it->end())
                {
                    return binding->second;
                }
            }
            // Not a binding, so it's a register reference
            if (local_index >= 0)
            {
                return {circuit.get_qubit_index(reg_name, local_index)};
            }
            else
            {
                size_t reg_size = circuit.get_qubit_reg_size(reg_name);
                std::vector<size_t> indices;

                for (size_t i = 0; i < reg_size; i++)
                {
                    indices.push_back(circuit.get_qubit_index(reg_name, i));
                }
                return indices;
            }
        }

        // Evaluate an expression to get a numeric value
        double evaluate_expr(const Expr *expr)
        {
            Value value = expr->evaluate();
            return value.as_double();
        }

        // Process an expression that refers to a qubit and return its global index
        std::vector<size_t> process_qubit_expr(const Expr *expr)
        {
            if (auto *index_expr = dynamic_cast<const IndexExpr *>(expr))
            {
                const std::string &reg_name = index_expr->get_name();
                int local_index = static_cast<int>(evaluate_expr(index_expr->get_index()));
                return resolve_qubit_index(reg_name, local_index);
            }
            else if (auto *var_expr = dynamic_cast<const VariableExpr *>(expr))
            {
                // This is a variable reference, check if it's a binding
                for (auto it = qubit_binding_stack.rbegin(); it != qubit_binding_stack.rend(); ++it)
                {
                    auto binding = it->find(var_expr->get_name());
                    if (binding != it->end())
                    {
                        return binding->second;
                    }
                }
                // Not found, this is a whole register reference
                return resolve_qubit_index(var_expr->get_name(), -1);
            }

            throw std::runtime_error("Invalid qubit expression");
        }

        // Process an expression that refers to a classical bit and return its global index
        std::vector<size_t> process_bit_expr(const Expr *expr)
        {
            if (auto *index_expr = dynamic_cast<const IndexExpr *>(expr))
            {
                const std::string &reg_name = index_expr->get_name();
                size_t local_index = static_cast<size_t>(evaluate_expr(index_expr->get_index()));
                return {circuit.get_bit_index(reg_name, local_index)};
            }
            else if (auto *var_expr = dynamic_cast<const VariableExpr *>(expr))
            {
                // This is a variable reference, likely a whole register
                std::vector<size_t> indices;
                size_t creg_size = circuit.get_bit_reg_size(var_expr->get_name());

                for (size_t i = 0; i < creg_size; i++)
                {
                    indices.push_back(circuit.get_bit_index(var_expr->get_name(), i));
                }
                return indices;
            }

            throw std::runtime_error("Invalid classical bit expression");
        }

        // Process a built-in gate statement
        void process_builtin_gate(const GateStmt *gate_stmt)
        {
            std::string gate_name = gate_stmt->get_name();
            std::transform(gate_name.begin(), gate_name.end(), gate_name.begin(),
                           [](unsigned char c)
                           { return std::tolower(c); });

            // Get qubit indices
            std::vector<std::vector<size_t>> qubit_indices;
            size_t max_qubit_reg_size = 0;
            for (const auto &qubit : gate_stmt->get_qubits())
            {
                qubit_indices.push_back(process_qubit_expr(qubit.get()));
                max_qubit_reg_size = std::max(max_qubit_reg_size, qubit_indices.back().size());
            }

            // Get parameters if any
            std::vector<double> parameters;
            for (const auto &param : gate_stmt->get_parameters())
            {
                parameters.push_back(evaluate_expr(param.get()));
            }

            // Create the operation
            auto op_type = Operation::name_to_type(gate_name);

            for (size_t i = 0; i < max_qubit_reg_size; i++)
            {
                std::vector<size_t> sub_qubit_indices;
                for (size_t j = 0; j < qubit_indices.size(); j++)
                {
                    if (qubit_indices[j].size() > 1)
                    {
                        sub_qubit_indices.push_back(qubit_indices[j][i]);
                    }
                    else
                    {
                        sub_qubit_indices.push_back(qubit_indices[j][0]);
                    }
                }

                if (parameters.empty())
                {
                    circuit.add_operation(Operation(op_type, sub_qubit_indices));
                }
                else
                {
                    circuit.add_operation(Operation(op_type, sub_qubit_indices, parameters));
                }
            }
        }

        // Process a user-defined gate statement (expand it)
        void process_user_defined_gate(const std::string &gate_name, const GateStmt *gate_stmt)
        {
            auto it = gate_definitions.find(gate_name);
            if (it == gate_definitions.end())
            {
                throw std::runtime_error("Unknown gate: " + gate_name);
            }

            const gate_defination &gate_def = it->second;

            // Get actual qubit indices
            std::vector<std::vector<size_t>> qubit_indices;
            for (const auto &qubit : gate_stmt->get_qubits())
            {
                qubit_indices.push_back(process_qubit_expr(qubit.get()));
            }

            // Bind the formal parameters to actual qubits
            push_qubit_bindings(gate_def.qubits, qubit_indices);

            // Process the gate body
            for (const auto &stmt : gate_def.body)
            {
                process_stmt(stmt.get());
            }

            // Clean up bindings after gate expansion
            pop_qubit_bindings();
        }

        // Process a gate application statement
        void process_gate_stmt(const GateStmt *gate_stmt)
        {
            std::string gate_name = gate_stmt->get_name();
            std::transform(gate_name.begin(), gate_name.end(), gate_name.begin(),
                           [](unsigned char c)
                           { return std::tolower(c); });

            // Check if it's a built-in gate
            if (Operation::is_builtin_gate(gate_name))
            {
                process_builtin_gate(gate_stmt);
            }
            else
            {
                // It's a user-defined gate, expand it
                process_user_defined_gate(gate_name, gate_stmt);
            }
        }

        // Process a measurement statement
        void process_measure_stmt(const MeasureStmt *measure_stmt)
        {
            std::vector<size_t> qubit_indices = process_qubit_expr(measure_stmt->get_qubit());
            std::vector<size_t> bit_indices = process_bit_expr(measure_stmt->get_bit());

            size_t num_measurements = std::max(qubit_indices.size(), bit_indices.size());

            for (size_t i = 0; i < num_measurements; i++)
            {
                size_t qubit_index = qubit_indices.size() > 1 ? qubit_indices[i] : qubit_indices[0];
                size_t bit_index = bit_indices.size() > 1 ? bit_indices[i] : bit_indices[0];

                circuit.add_operation(Operation(Operation::Type::MEASURE, {qubit_index}, {}, {bit_index}));
            }
        }

        // Process a reset statement
        void process_reset_stmt(const ResetStmt *reset_stmt)
        {
            std::vector<size_t> qubit_indices = process_qubit_expr(reset_stmt->get_qubit());

            for (size_t i = 0; i < qubit_indices.size(); i++)
            {
                circuit.add_operation(Operation(Operation::Type::RESET, {qubit_indices[i]}));
            }
        }

        // Process a barrier statement
        void process_barrier_stmt(const BarrierStmt *barrier_stmt)
        {
            std::vector<size_t> qubit_indices;

            for (const auto &qubit : barrier_stmt->get_qubits())
            {
                std::vector<size_t> sub_qubit_indices = process_qubit_expr(qubit.get());
                qubit_indices.insert(qubit_indices.end(), sub_qubit_indices.begin(), sub_qubit_indices.end());
            }

            // Create barrier operation
            circuit.add_operation(Operation(Operation::Type::BARRIER, qubit_indices));
        }

        // Process a Pauli gate statement
        void process_pauli_stmt(const PauliStmt *pauli_stmt)
        {
            std::string gate_name = pauli_stmt->get_gate_name();
            std::string pauli_string = pauli_stmt->get_pauli_string();

            // Determine the operation type
            Operation::Type op_type;
            if (gate_name == "t_pauli")
            {
                op_type = Operation::Type::T_PAULI;
            }
            else if (gate_name == "m_pauli")
            {
                op_type = Operation::Type::M_PAULI;
            }
            else if (gate_name == "s_pauli")
            {
                op_type = Operation::Type::S_PAULI;
            }
            else if (gate_name == "z_pauli")
            {
                op_type = Operation::Type::Z_PAULI;
            }
            else
            {
                throw std::runtime_error("Unknown Pauli gate type: " + gate_name);
            }

            // Create PauliOp from string and add operation
            PauliOp pauli_op(circuit.get_num_qubits());
            pauli_op.from_string(pauli_string);
            circuit.add_operation(Operation(op_type, {}, {}, {}, pauli_op));
        }

        // Process a statement in the AST
        void process_stmt(const Stmt *stmt)
        {
            if (auto *version_decl = dynamic_cast<const VersionDecl *>(stmt))
            {
                // Version declaration - nothing to do for circuit building
                (void)version_decl;
            }
            else if (auto *include_stmt = dynamic_cast<const IncludeStmt *>(stmt))
            {
                // Include statement - for now, we assume standard includes are already processed
                (void)include_stmt;
            }
            else if (auto *qreg = dynamic_cast<const QRegDecl *>(stmt))
            {
                // Quantum register declaration
                circuit.add_qreg(qreg->get_name(), qreg->get_size());
            }
            else if (auto *creg = dynamic_cast<const CRegDecl *>(stmt))
            {
                // Classical register declaration
                circuit.add_creg(creg->get_name(), creg->get_size());
            }
            else if (auto *gate = dynamic_cast<const GateStmt *>(stmt))
            {
                // Gate application
                process_gate_stmt(gate);
            }
            else if (auto *measure = dynamic_cast<const MeasureStmt *>(stmt))
            {
                // Measurement
                process_measure_stmt(measure);
            }
            else if (auto *reset = dynamic_cast<const ResetStmt *>(stmt))
            {
                // Reset
                process_reset_stmt(reset);
            }
            else if (auto *barrier = dynamic_cast<const BarrierStmt *>(stmt))
            {
                // Barrier
                process_barrier_stmt(barrier);
            }
            else if (auto *pauli = dynamic_cast<const PauliStmt *>(stmt))
            {
                // Pauli gate (t_pauli, m_pauli, s_pauli)
                process_pauli_stmt(pauli);
            }
            else if (auto *gate_decl = dynamic_cast<const GateDeclStmt *>(stmt))
            {
                // Gate declaration
                gate_defination gate_def;
                gate_def.qubits = gate_decl->get_qubits();
                gate_def.params = gate_decl->get_params();

                // Deep copy the gate body
                for (const auto &body_stmt : gate_decl->get_body())
                {
                    gate_def.body.push_back(std::unique_ptr<Stmt>(body_stmt->clone()));
                }

                // Store the gate definition as lowercase
                std::string gate_name = gate_decl->get_name();
                std::transform(gate_name.begin(), gate_name.end(), gate_name.begin(),
                               [](unsigned char c)
                               { return std::tolower(c); });

                gate_definitions[gate_name] = std::move(gate_def);
            }
            else if (auto *block = dynamic_cast<const BlockStmt *>(stmt))
            {
                // Block statement
                for (const auto &s : block->get_statements())
                {
                    process_stmt(s.get());
                }
            }
            else if (auto *if_stmt = dynamic_cast<const IfStmt *>(stmt))
            {
                // If statement - not handling conditional execution in flattened circuit
                // For a complete implementation, this would need special handling
                (void)if_stmt;
                throw std::runtime_error("Conditional statements not yet supported in circuit flattening");
            }
        }

    public:
        // Build a circuit from a parsed program
        Circuit build(const ASTProgram *program)
        {
            circuit = Circuit(); // Reset circuit
            gate_definitions.clear();
            qubit_binding_stack.clear();

            // Process all statements in the program
            for (const auto &stmt : program->get_statements())
            {
                process_stmt(stmt.get());
            }

            return circuit;
        }
    };

} // namespace NWQEC