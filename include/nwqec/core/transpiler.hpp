#pragma once

#include "nwqec/core/circuit.hpp"
#include "nwqec/core/constants.hpp"
#include "nwqec/core/transpiler_passes.hpp"

#include "nwqec/passes/clifford_reduction_pass.hpp"
#include "nwqec/passes/pbc_pass.hpp"
#include "nwqec/passes/decompose_pass.hpp"
#include "nwqec/passes/remove_trivial_rz_pass.hpp"
#include "nwqec/passes/synthesize_rz_pass.hpp"
#include "nwqec/passes/gate_fusion_pass.hpp"
#include "nwqec/passes/tfuse_pass.hpp"
#include "nwqec/passes/remove_pauli_pass.hpp"

#include <string>
#include <vector>
#include <memory>
#include <optional>
#include <iostream>
#include <iomanip>

namespace NWQEC {

/**
 * @brief Configuration options for pass execution
 */
struct PassConfig {
    bool keep_ccx = false;          // Preserve CCX gates during decomposition
    bool keep_cx = false;           // Preserve CX gates in PBC format
    RzErrorPolicy rz_error_policy = DEFAULT_RZ_ERROR_POLICY;
    std::optional<double> rz_error_epsilon = std::nullopt;
    bool silent = false;            // Suppress output during pass execution
};

/**
 * @brief Core transpiler engine that executes sequences of passes
 */
class Transpiler {
public:
    Transpiler() = default;

    /**
     * @brief Execute a sequence of passes on a circuit
     * @param circuit Input circuit (will be modified)
     * @param passes Vector of passes to apply in order
     * @param config Configuration options for pass execution
     * @return Modified circuit
     */
    std::unique_ptr<Circuit> execute_passes(
        std::unique_ptr<Circuit> circuit,
        const std::vector<PassType>& passes,
        const PassConfig& config = {}
    );

    /**
     * @brief Execute a predefined pass sequence
     * @param circuit Input circuit
     * @param sequence Predefined sequence from PassSequences namespace
     * @param config Configuration options
     * @return Modified circuit
     */
    std::unique_ptr<Circuit> execute_sequence(
        std::unique_ptr<Circuit> circuit,
        const std::vector<PassType>& sequence,
        const PassConfig& config = {}
    ) {
        return execute_passes(std::move(circuit), sequence, config);
    }

private:
    /**
     * @brief Create and configure a pass instance
     */
    std::unique_ptr<Pass> create_pass(PassType type, const PassConfig& config);
    
    /**
     * @brief Print execution statistics for a pass
     */
    void print_pass_stats(const std::string& pass_name, size_t gates_before, const Circuit& after, bool modified, size_t depth_after);
    
    /**
     * @brief Print table header for pass execution log
     */
    void print_table_header();
};

// Inline implementations
inline std::unique_ptr<Circuit> Transpiler::execute_passes(
    std::unique_ptr<Circuit> circuit,
    const std::vector<PassType>& passes,
    const PassConfig& config
) {
    if (!config.silent) {
        std::cout << "\n=== Pass Execution Summary ===\n";
        print_table_header();
    }

    std::optional<size_t> cached_depth;

    for (PassType pass_type : passes) {
        auto pass = create_pass(pass_type, config);
        if (!pass) {
            if (!config.silent) {
                std::cerr << "Warning: Unknown pass type: " << pass_type_to_string(pass_type) << std::endl;
            }
            continue;
        }

        const size_t gates_before = circuit->get_operations().size();
        bool modified = pass->run(*circuit);
        
        if (!config.silent) {
            if (modified || !cached_depth.has_value()) {
                cached_depth = circuit->depth();
            }
            print_pass_stats(pass_type_to_string(pass_type), gates_before, *circuit, modified, *cached_depth);
        }
    }

    if (!config.silent) {
        std::cout << "\n=== Final Statistics ===\n";
        circuit->print_stats(std::cout);
    }

    return circuit;
}

inline std::unique_ptr<Pass> Transpiler::create_pass(PassType type, const PassConfig& config) {
    switch (type) {
        case PassType::DECOMPOSE:
            return std::make_unique<DecomposePass>(config.keep_ccx);
        
        case PassType::REMOVE_TRIVIAL_RZ:
            return std::make_unique<RemoveTrivialRzPass>();
        
        case PassType::GATE_FUSION:
            return std::make_unique<GateFusionPass>();
        
        case PassType::REMOVE_PAULI:
            return std::make_unique<RemovePauliPass>();
        
        case PassType::TO_PBC:
            return std::make_unique<PbcPass>(config.keep_cx);
        
        case PassType::CLIFFORD_REDUCTION:
            return std::make_unique<CRPass>();
        
        case PassType::SYNTHESIZE_RZ: {
            double epsilon = config.rz_error_epsilon.value_or(default_rz_error_epsilon(config.rz_error_policy));
            return std::make_unique<SynthesizeRzPass>(config.rz_error_policy, epsilon);
        }
        
        case PassType::TFUSE:
            return std::make_unique<TfusePass>();
        
        default:
            return nullptr;
    }
}

inline void Transpiler::print_table_header() {
    std::cout << std::left
              << std::setw(25) << "Pass"
              << std::setw(10) << "Modified"
              << std::setw(15) << "Gates Before"
              << std::setw(15) << "Gates After"
              << std::setw(10) << "Depth"
              << std::endl;
    std::cout << std::string(75, '-') << std::endl;
}

inline void Transpiler::print_pass_stats(const std::string& pass_name, size_t gates_before, const Circuit& after, bool modified, size_t depth_after) {
    const size_t gates_after = after.get_operations().size();

    std::cout << std::left
              << std::setw(25) << pass_name
              << std::setw(10) << (modified ? "Yes" : "No")
              << std::setw(15) << gates_before
              << std::setw(15) << gates_after
              << std::setw(10) << depth_after
              << std::endl;
}

} // namespace NWQEC
