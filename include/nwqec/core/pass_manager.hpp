#ifndef NWQEC_WITH_GRIDSYNTH_CPP
#define NWQEC_WITH_GRIDSYNTH_CPP 0
#endif

#pragma once

#include "nwqec/core/circuit.hpp"

#include "nwqec/passes/clifford_reduction_pass.hpp"
#include "nwqec/passes/pbc_pass.hpp"

#include "nwqec/passes/decompose_pass.hpp"
#include "nwqec/passes/remove_trivial_rz_pass.hpp"
#if NWQEC_WITH_GRIDSYNTH_CPP
#include "nwqec/passes/synthesize_rz_pass.hpp"
#endif
#include "nwqec/passes/gate_fusion_pass.hpp"
#include "nwqec/passes/tfuse_pass.hpp"
#include "nwqec/passes/remove_pauli_pass.hpp"

#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>

namespace NWQEC
{

    class PassManager
    {
    public:
        PassManager() = default;

        std::unique_ptr<Circuit> apply_passes(std::unique_ptr<Circuit> circuit,
                                              bool to_pbc = false,
                                              bool to_clifford_reduction = false,
                                              bool keep_cx = false,
                                              bool t_pauli_opt = false,
                                              bool remove_pauli = false,
                                              bool keep_ccx = false,
                                              bool silent = false,
                                              double epsilon_override = -1.0)
        {
            if (!silent)
            {
                std::cout << "\n=== Circuit Transpilation Summary ===\n";
            }

            if (!circuit->is_clifford_t())
            {
                if (!silent)
                {
                    std::cout << "\n--- Transpiling to Clifford+T ---\n";
                    print_table_header();
                }

                // Skip CCX decomposition if PBC or Clifford Reduction is enabled, or user explicitly requested to keep CCX gates
                bool keep_ccx_gates = to_pbc || to_clifford_reduction || keep_ccx;
                auto transpilation_passes = get_general_ct_transpilation(keep_ccx_gates, epsilon_override);

                for (const auto &pass_factory : transpilation_passes)
                {
                    auto pass = pass_factory();
                    bool modified = pass->run(*circuit);
                    if (!silent)
                    {
                        print_table_row(pass->get_name(), modified ? "YES" : "NO");
                    }
                }
                if (!silent)
                {
                    print_table_footer();
                    print_circuit_stats(*circuit, "After Clifford+T transpilation");
                }
            }

            // check that only one of to_pbc or to_clifford_reduction is true, all can be false
            int selected_passes = (to_pbc ? 1 : 0) + (to_clifford_reduction ? 1 : 0);
            if (selected_passes > 1)
            {
                throw std::runtime_error("Cannot transpile to multiple passes at the same time. Only one of PBC or Clifford Reduction can be enabled.");
            }

            if (to_pbc)
            {
                if (!silent)
                {
                    std::cout << "\n--- Transpiling to PBC ---\n";
                    print_table_header();
                }
                auto pbc_pass = std::make_unique<PbcPass>(keep_cx);
                bool modified = pbc_pass->run(*circuit);
                if (!silent)
                {
                    print_table_row(pbc_pass->get_name(), modified ? "YES" : "NO");
                }
                if (!silent)
                {
                    print_table_footer();
                    print_circuit_stats(*circuit, "After PBC transpilation");
                }
            }
            else if (to_clifford_reduction)
            {
                if (!silent)
                {
                    std::cout << "\n--- Applying Clifford Reduction ---\n";
                    print_table_header();
                }
                for (const auto &pass_factory : clifford_reduction_transpilation)
                {
                    auto pass = pass_factory();
                    bool modified = pass->run(*circuit);
                    if (!silent)
                    {
                        print_table_row(pass->get_name(), modified ? "YES" : "NO");
                    }
                }
                if (!silent)
                {
                    print_table_footer();
                    print_circuit_stats(*circuit, "After Clifford Reduction transpilation");
                }
            }

            // Apply T Pauli optimizer if enabled after PBC
            if (t_pauli_opt)
            {
                if (!silent)
                {
                    std::cout << "\n--- Applying T Pauli Optimizer ---\n";
                }

                // Count T_PAULI operations before optimization
                size_t t_pauli_count_before = count_t_pauli_operations(*circuit);

                if (!silent)
                {
                    print_table_header();
                }
                auto pass = std::make_unique<TfusePass>();
                bool modified = pass->run(*circuit);
                if (!silent)
                {
                    print_table_row(pass->get_name(), modified ? "YES" : "NO");
                    print_table_footer();

                    // Count T_PAULI operations after optimization and print statistics
                    size_t t_pauli_count_after = count_t_pauli_operations(*circuit);
                    print_t_pauli_reduction_stats(t_pauli_count_before, t_pauli_count_after);

                    print_circuit_stats(*circuit, "After T Pauli optimization");
                }
            }

            // Apply RemovePauliPass at the end if remove_pauli is enabled
            if (remove_pauli)
            {
                if (!silent)
                {
                    std::cout << "\n--- Applying Pauli Gate Removal ---\n";
                    print_table_header();
                }
                auto pauli_pass = std::make_unique<RemovePauliPass>();
                bool modified = pauli_pass->run(*circuit);
                if (!silent)
                {
                    print_table_row(pauli_pass->get_name(), modified ? "YES" : "NO");
                    print_table_footer();
                    print_circuit_stats(*circuit, "After Pauli gate removal");
                }
            }

            return circuit;
        }

    private:
        using PassFactory = std::function<std::unique_ptr<Pass>()>;

        // Helper functions for table formatting
        void print_table_header(bool show_modified = true) const
        {
            std::string separator = show_modified ? "┌─────────────────────────────────┬──────────────┐" : "┌─────────────────────────────────┐";
            std::cout << separator << "\n";

            std::string header = show_modified ? "│ Pass Name                       │ Modified     │" : "│ Pass Name                       │";
            std::cout << header << "\n";

            separator = show_modified ? "├─────────────────────────────────┼──────────────┤" : "├─────────────────────────────────┤";
            std::cout << separator << "\n";
        }

        void print_table_row(const std::string &pass_name, const std::string &status, bool show_modified = true) const
        {
            std::stringstream ss;
            if (show_modified)
            {
                ss << "│ " << std::left << std::setw(31) << pass_name
                   << " │ " << std::setw(12) << status << " │";
            }
            else
            {
                ss << "│ " << std::left << std::setw(31) << pass_name << " │";
            }
            std::cout << ss.str() << "\n";
        }

        void print_table_footer(bool show_modified = true) const
        {
            std::string separator = show_modified ? "└─────────────────────────────────┴──────────────┘" : "└─────────────────────────────────┘";
            std::cout << separator << "\n";
        }

        void print_circuit_stats(const Circuit &circuit, const std::string &stage) const
        {
            std::cout << "\n--- " << stage << " ---\n";
            circuit.print_stats(std::cout);
        }

        std::vector<PassFactory> get_general_ct_transpilation(bool keep_ccx = false, double epsilon_override = -1.0) const
        {
            std::vector<PassFactory> passes = {
                [keep_ccx]()
                { return std::make_unique<DecomposePass>(keep_ccx); },
                // []() { return std::make_unique<GateFusionPass>(); },
                []()
                { return std::make_unique<RemoveTrivialRzPass>(); },
            };
#if NWQEC_WITH_GRIDSYNTH_CPP
            passes.push_back([epsilon_override]()
                             { return std::make_unique<SynthesizeRzPass>(epsilon_override); });
#endif
            return passes;
        }

        // []()
        // { return std::make_unique<GateFusionPass>(); }

        std::vector<PassFactory> clifford_reduction_transpilation = {
            []()
            { return std::make_unique<CRPass>(); }};

        size_t count_t_pauli_operations(const Circuit &circuit) const
        {
            size_t count = 0;
            const auto &operations = circuit.get_operations();
            for (const auto &op : operations)
            {
                if (op.get_type() == Operation::Type::T_PAULI)
                {
                    count++;
                }
            }
            return count;
        }

        void print_t_pauli_reduction_stats(size_t before, size_t after) const
        {
            size_t reduction = (before >= after) ? (before - after) : 0;
            double percentage = (before > 0) ? (static_cast<double>(reduction) / before * 100.0) : 0.0;

            std::cout << "\n--- T-Pauli Reduction Statistics ---\n";
            std::cout << "Original T-Pauli count:  " << before << "\n";
            std::cout << "Optimized T-Pauli count: " << after << "\n";
            std::cout << "Reduction:               " << reduction << "\n";
            std::cout << "Reduction percentage:    " << std::fixed << std::setprecision(2) << percentage << "%\n";
        }
    };

}
