#pragma once

#include <vector>
#include <string>

namespace NWQEC {

/**
 * @brief Enumeration of available transpilation passes
 */
enum class PassType {
    // Basic decomposition and cleanup
    DECOMPOSE,              // Decompose gates to basic gate set
    REMOVE_TRIVIAL_RZ,      // Remove RZ gates with zero or trivial angles
    GATE_FUSION,            // Fuse adjacent gates where possible
    REMOVE_PAULI,           // Remove Pauli gates from output
    
    // Circuit format conversions  
    TO_PBC,                 // Convert to Pauli-Based Circuit format
    CLIFFORD_REDUCTION,     // Apply Clifford reduction (TACO) pipeline
    
    // RZ synthesis (requires gridsynth backend)
    SYNTHESIZE_RZ,          // Synthesize RZ gates to Clifford+T (needs GMP/MPFR)
    
    // PBC-specific optimizations
    TFUSE                   // T-count optimization for PBC circuits
};

/**
 * @brief Convert PassType enum to string representation
 */
inline std::string pass_type_to_string(PassType type) {
    switch (type) {
        case PassType::DECOMPOSE:           return "DECOMPOSE";
        case PassType::REMOVE_TRIVIAL_RZ:   return "REMOVE_TRIVIAL_RZ";
        case PassType::GATE_FUSION:         return "GATE_FUSION";
        case PassType::REMOVE_PAULI:        return "REMOVE_PAULI";
        case PassType::TO_PBC:              return "TO_PBC";
        case PassType::CLIFFORD_REDUCTION:  return "CLIFFORD_REDUCTION";
        case PassType::SYNTHESIZE_RZ:       return "SYNTHESIZE_RZ";
        case PassType::TFUSE:               return "TFUSE";
        default:                            return "UNKNOWN";
    }
}

/**
 * @brief Predefined pass sequences for common workflows
 */
namespace PassSequences {
    // Convert to Clifford+T (requires gridsynth)
    inline const std::vector<PassType> TO_CLIFFORD_T = {
        PassType::DECOMPOSE,
        PassType::REMOVE_TRIVIAL_RZ,
        PassType::SYNTHESIZE_RZ,  // Requires gridsynth
        PassType::GATE_FUSION
    };
    
    // Convert to Clifford+T+RZ (no gridsynth needed - stops before RZ synthesis)
    inline const std::vector<PassType> TO_CLIFFORD_T_RZ = {
        PassType::DECOMPOSE,
        PassType::REMOVE_TRIVIAL_RZ
    };
    
    // Convert to PBC format
    inline const std::vector<PassType> TO_PBC_BASIC = {
        PassType::DECOMPOSE,
        PassType::REMOVE_TRIVIAL_RZ,
        PassType::SYNTHESIZE_RZ,  // Requires gridsynth
        PassType::TO_PBC
    };
    
    // PBC with T-count optimization  
    inline const std::vector<PassType> TO_PBC_OPTIMIZED = {
        PassType::DECOMPOSE,
        PassType::REMOVE_TRIVIAL_RZ,
        PassType::SYNTHESIZE_RZ,  // Requires gridsynth
        PassType::TO_PBC,
        PassType::TFUSE
    };
    
    // Clifford reduction (TACO) pipeline
    inline const std::vector<PassType> CLIFFORD_REDUCTION = {
        PassType::DECOMPOSE,
        PassType::REMOVE_TRIVIAL_RZ,
        PassType::SYNTHESIZE_RZ,  // Requires gridsynth
        PassType::CLIFFORD_REDUCTION
    };
    
    // Post-synthesis cleanup (after manual RZ synthesis)
    inline const std::vector<PassType> POST_SYNTHESIS_CLEANUP = {
        PassType::GATE_FUSION,
        PassType::REMOVE_TRIVIAL_RZ
    };
}

} // namespace NWQEC