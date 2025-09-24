#pragma once

#include "nwqec/core/circuit.hpp" // Parent class
#include <vector>
#include <string>    // For std::string, std::to_string
#include <algorithm> // For std::find
#include <stdexcept> // For std::out_of_range
#include <limits>    // For std::numeric_limits
#include <numeric>   // For std::iota (if used, not currently)
#include <memory>
#include <fstream>

namespace NWQEC
{
    // A struct to hold the predecessors and successors of an operation that contains the qubit index and the operation index
    struct OperationDependency
    {
        size_t qubit; // The qubit index this operation depends on
        size_t node;  // The operation index that this operation depends on

        OperationDependency(size_t q, size_t n) : qubit(q), node(n) {}
    };

    class DAGCircuit : public Circuit
    {

    public:
        DAGCircuit() : Circuit()
        {
            // last_op_on_qubit is initially empty (default constructed).
            // It will be sized by add_qreg or rebuild_dag (via clear_dag_structure).
        }

        // New constructor to create a DAGCircuit from an existing Circuit object
        DAGCircuit(const Circuit &existing_circuit) : Circuit(existing_circuit)
        {
            // The base Circuit part of this DAGCircuit is now a copy of existing_circuit.
            // This includes num_qubits, operations, register maps, etc.
            // Now, build the DAG structures based on these copied operations.
            // rebuild_dag() will call clear_dag_structure(), which correctly initializes
            // last_op_on_qubit based on get_num_qubits() from the copied circuit.

            // check if it's already a DAGCircuit
            if (dynamic_cast<const DAGCircuit *>(&existing_circuit) != nullptr)
            {
                // If it's already a DAGCircuit, we can just copy the DAG structure
                const DAGCircuit &dag_circuit = static_cast<const DAGCircuit &>(existing_circuit);
                successors = dag_circuit.successors;
                predecessors = dag_circuit.predecessors;
                last_op_on_qubit = dag_circuit.last_op_on_qubit;
            }
            else
            {
                // Otherwise, we need to rebuild the DAG structure from the operations
                rebuild_dag();
            }
        }

        // Override add_qreg to initialize/resize last_op_on_qubit
        void add_qreg(const std::string &name, size_t size) override
        {
            Circuit::add_qreg(name, size);
            // Resize last_op_on_qubit to the new total num_qubits.
            // New entries for newly added qubits are initialized to SIZE_MAX.
            // Existing entries for older qubits retain their last_op_idx.
            last_op_on_qubit.resize(get_num_qubits(), std::numeric_limits<size_t>::max());
        }

        // Override add_operation to build/update the DAG
        void add_operation(Operation operation) override
        {
            // Capture qubit data before moving the operation object
            const std::vector<size_t> &op_qubits = operation.get_qubits();
            size_t new_op_idx = get_operations().size(); // Index this operation will have

            // Validate qubit indices against declared qubits
            for (size_t qubit_idx : op_qubits)
            {
                if (qubit_idx >= get_num_qubits())
                {
                    throw std::out_of_range("DAGCircuit::add_operation: Qubit index " + std::to_string(qubit_idx) +
                                            " (in op " + operation.get_type_name() + ")" +
                                            " is out of range for declared qubits (" + std::to_string(get_num_qubits()) +
                                            "). Ensure qregs cover all used qubits before adding operations.");
                }
                // Also ensure last_op_on_qubit is large enough (should be if add_qreg was called)
                if (qubit_idx >= last_op_on_qubit.size())
                {
                    // This might happen if add_operation is called before add_qreg for those qubits
                    // Or if num_qubits was manipulated in a way that last_op_on_qubit didn't catch up.
                    // Forcing a resize here could be an option, but strict adherence to add_qreg first is safer.
                    throw std::logic_error("DAGCircuit::add_operation: Qubit index " + std::to_string(qubit_idx) +
                                           " encountered but last_op_on_qubit not properly sized. Call add_qreg first.");
                }
            }

            // Resize graph structures for the new operation
            if (new_op_idx >= successors.size())
            {
                successors.resize(new_op_idx + 1);
            }
            else
            {
                successors[new_op_idx].clear(); // Clear any old data if resizing didn't occur
            }
            if (new_op_idx >= predecessors.size())
            {
                predecessors.resize(new_op_idx + 1);
            }
            else
            {
                predecessors[new_op_idx].clear();
            }

            // Determine dependencies based on qubit usage
            for (size_t qubit_idx : op_qubits)
            {
                size_t prev_op_idx_on_this_qubit = last_op_on_qubit[qubit_idx];

                if (prev_op_idx_on_this_qubit != std::numeric_limits<size_t>::max())
                {
                    successors[prev_op_idx_on_this_qubit].emplace_back(qubit_idx, new_op_idx);
                    predecessors[new_op_idx].emplace_back(qubit_idx, prev_op_idx_on_this_qubit);
                }
                // Update this qubit to now be last touched by the new operation
                last_op_on_qubit[qubit_idx] = new_op_idx;
            }

            // Add the operation to the parent class's storage
            Circuit::add_operation(std::move(operation));
        }

        // Accessor for successors of an operation
        const std::vector<OperationDependency> &get_successors(size_t op_idx) const
        {
            if (op_idx >= successors.size())
            {
                throw std::out_of_range("DAGCircuit::get_successors: op_idx " + std::to_string(op_idx) + " out of range.");
            }
            return successors[op_idx];
        }

        // Accessor for predecessors of an operation
        const std::vector<OperationDependency> &get_predecessors(size_t op_idx) const
        {
            if (op_idx >= predecessors.size())
            {
                throw std::out_of_range("DAGCircuit::get_predecessors: op_idx " + std::to_string(op_idx) + " out of range.");
            }
            return predecessors[op_idx];
        }

        // Get all root nodes (operations with no predecessors)
        std::vector<size_t> get_root_nodes() const
        {
            std::vector<size_t> roots;
            size_t num_actual_ops = get_operations().size();
            for (size_t i = 0; i < predecessors.size() && i < num_actual_ops; ++i)
            {
                if (predecessors[i].empty())
                {
                    roots.push_back(i);
                }
            }
            return roots;
        }

        // Get all leaf nodes (operations with no successors)
        std::vector<size_t> get_leaf_nodes() const
        {
            std::vector<size_t> leaves;
            size_t num_actual_ops = get_operations().size();
            for (size_t i = 0; i < successors.size() && i < num_actual_ops; ++i)
            {
                if (successors[i].empty())
                {
                    leaves.push_back(i);
                }
            }
            return leaves;
        }

        // Clears the DAG structure. Operations in the base Circuit remain.
        void clear_dag_structure()
        {
            successors.clear();
            predecessors.clear();
            if (get_num_qubits() > 0)
            {
                last_op_on_qubit.assign(get_num_qubits(), std::numeric_limits<size_t>::max());
            }
            else
            {
                last_op_on_qubit.clear();
            }
        }

        // Rebuilds the DAG from operations currently in the base Circuit.
        // Useful if operations were added/modified without going through DAGCircuit's add_operation,
        // or after loading a circuit.
        void rebuild_dag()
        {
            clear_dag_structure(); // Reset DAG state and last_op_on_qubit

            const auto &ops = get_operations();
            if (ops.empty())
            {
                return; // Nothing to build
            }

            successors.resize(ops.size());
            predecessors.resize(ops.size());
            // last_op_on_qubit should be correctly sized by clear_dag_structure if get_num_qubits() > 0
            // If get_num_qubits() is 0 but ops exist with qubits, it's an inconsistent state.
            // We rely on get_num_qubits() being the authority for qubit count.

            for (size_t op_idx = 0; op_idx < ops.size(); ++op_idx)
            {
                const auto &current_op_qubits = ops[op_idx].get_qubits();

                // Ensure internal consistency for qubit indices during rebuild
                for (size_t qubit_idx : current_op_qubits)
                {
                    if (qubit_idx >= get_num_qubits())
                    {
                        throw std::out_of_range("DAGCircuit::rebuild_dag: Qubit index " + std::to_string(qubit_idx) +
                                                " in op " + std::to_string(op_idx) + " (" + ops[op_idx].get_type_name() + ")" +
                                                " is out of range for declared qubits (" + std::to_string(get_num_qubits()) + ").");
                    }
                    if (qubit_idx >= last_op_on_qubit.size())
                    {
                        // Should not happen if last_op_on_qubit is sized by get_num_qubits()
                        throw std::logic_error("DAGCircuit::rebuild_dag: Qubit index " + std::to_string(qubit_idx) +
                                               " encountered but last_op_on_qubit not properly sized. Inconsistent state.");
                    }
                }

                // Process dependencies for the current operation
                for (size_t qubit_idx : current_op_qubits)
                {
                    size_t prev_op_idx = last_op_on_qubit[qubit_idx];
                    if (prev_op_idx != std::numeric_limits<size_t>::max())
                    {
                        // There is a previous operation on this qubit
                        successors[prev_op_idx].emplace_back(qubit_idx, op_idx);
                        predecessors[op_idx].emplace_back(qubit_idx, prev_op_idx);
                    }
                    last_op_on_qubit[qubit_idx] = op_idx;
                }
            }
        }

        /**
         * @brief Generate a visual representation of the DAG
         *
         * @param filename Base filename for output (without extension)
         * @throws std::invalid_argument if filename is empty
         * @throws std::runtime_error if file operations fail
         */
        void draw(const std::string &filename) const
        {
            if (filename.empty())
            {
                throw std::invalid_argument("Filename cannot be empty");
            }

            const std::string dot_filename = filename + ".dot";
            const std::string png_filename = filename + ".png";

            if (!write_dot_file(dot_filename))
            {
                throw std::runtime_error("Failed to write DOT file: " + dot_filename);
            }

            if (!generate_graph_image(dot_filename, png_filename))
            {
                throw std::runtime_error("Failed to generate graph image. Make sure Graphviz is installed.");
            }
        }

    private:
        // Adjacency lists for the DAG. Indices correspond to operation indices in 'operations' vector.
        std::vector<std::vector<OperationDependency>> successors;   // successors[i] = list of ops that depend on op i
        std::vector<std::vector<OperationDependency>> predecessors; // predecessors[i] = list of ops that op i depends on

        // Tracks the last operation index that affected each qubit (global qubit index to op index)
        std::vector<size_t> last_op_on_qubit;

        /**
         * @brief Write the DAG structure to a DOT file for visualization
         *
         * @param filename Output DOT filename
         * @return bool True if successful, false otherwise
         */
        bool write_dot_file(const std::string &filename) const
        {
            std::ofstream dot_file(filename);
            if (!dot_file)
            {
                std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
                return false;
            }

            dot_file << "digraph G {\n"
                     << "    node [style=filled, shape=ellipse];\n";

            auto ops = get_operations();
            size_t num_ops = ops.size();
            size_t num_qubits = get_num_qubits();

            // Add all operation nodes
            for (size_t i = 0; i < num_ops; ++i)
            {
                if (ops[i].get_type() == Operation::Type::MEASURE)
                {
                    dot_file << "    " << i << " [fillcolor=\"red\", label=\"measure q["
                             << ops[i].get_qubits()[0] << "]\"];\n";
                }
                else
                {
                    auto qubits = ops[i].get_qubits();
                    std::string qubit_label = " ";
                    for (size_t j = 0; j < qubits.size(); ++j)
                    {
                        qubit_label += "q[" + std::to_string(qubits[j]) + "]";
                        if (j < qubits.size() - 1)
                        {
                            qubit_label += ",";
                        }
                    }
                    dot_file << "    " << i << " [fillcolor=\"lightblue\", label=\"" << ops[i].get_type_name()
                             << qubit_label << "\"];\n";
                }
            }

            // Manually add input nodes with node index starting from num_ops
            // This is to ensure that input nodes are indexed after all operation nodes
            for (size_t i = 0; i < num_qubits; ++i)
            {
                dot_file << "    " << num_ops + i << " [fillcolor=\"lightgreen\", label=\"q[" << i << "]\"];\n";
            }

            // Add edges
            for (size_t i = 0; i < num_ops; ++i)
            {
                std::vector<OperationDependency> node_predecessors = predecessors[i];

                std::vector<size_t> node_predecessors_qubits;
                for (const OperationDependency &dep : node_predecessors)
                {
                    node_predecessors_qubits.push_back(dep.qubit);
                }
                // Get qubit indices that should be linked from input nodes

                for (size_t qubit_idx : ops[i].get_qubits())
                {
                    if (std::find(node_predecessors_qubits.begin(), node_predecessors_qubits.end(), qubit_idx) == node_predecessors_qubits.end())
                    {
                        dot_file << "    " << num_ops + qubit_idx << " -> " << i;
                        // Label the edge as q[qubit]
                        dot_file << " [label=\"q[" << qubit_idx << "]\"];\n";
                    }
                }

                for (OperationDependency op_dep : node_predecessors)
                {
                    dot_file << "    " << op_dep.node << " -> " << i;
                    // Label the edge as q[qubit]
                    dot_file << " [label=\"q[" << op_dep.qubit << "]\"];\n";
                }
            }

            dot_file << "}\n";
            return true;
        }

        /**
         * @brief Generate a PNG image from a DOT file using Graphviz
         *
         * @param dot_file Input DOT filename
         * @param output_file Output PNG filename
         * @return bool True if successful, false otherwise
         */
        bool generate_graph_image(const std::string &dot_file, const std::string &output_file) const
        {
            const std::string cmd = "dot -Tpng " + dot_file + " -o " + output_file;
            return system(cmd.c_str()) == 0;
        }
    };

} // namespace NWQEC
