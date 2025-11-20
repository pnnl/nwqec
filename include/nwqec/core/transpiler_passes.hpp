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
    CLIFFORD_REDUCTION,     // Apply Clifford reduction optimization
    
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
 * 
 * The transpilation workflow follows this pattern:
 * 1. Basic Processing: DECOMPOSE → REMOVE_TRIVIAL_RZ → SYNTHESIZE_RZ
 * 2. Choose transformation: Clifford+T, PBC, or Clifford Reduction
 * 3. Optional optimizations (T-fusion for PBC, cleanup passes)
 */
namespace PassSequences {
    // === BASIC PROCESSING ===
    
    // Standard preprocessing: decompose and clean up trivial RZ gates
    inline const std::vector<PassType> BASIC_PREPROCESSING = {
        PassType::DECOMPOSE,
        PassType::REMOVE_TRIVIAL_RZ
    };
    
    // Full preprocessing including RZ synthesis to Clifford+T
    inline const std::vector<PassType> FULL_PREPROCESSING = {
        PassType::DECOMPOSE,
        PassType::REMOVE_TRIVIAL_RZ,
        PassType::SYNTHESIZE_RZ
    };
    
    // === COMPLETE WORKFLOWS ===
    
    // Convert to Clifford+T with final gate fusion
    inline const std::vector<PassType> TO_CLIFFORD_T = {
        PassType::DECOMPOSE,
        PassType::REMOVE_TRIVIAL_RZ,
        PassType::SYNTHESIZE_RZ,
        PassType::GATE_FUSION
    };
    
    // Convert to PBC format (synthesize RZ first, then convert to PBC)
    inline const std::vector<PassType> TO_PBC = {
        PassType::DECOMPOSE,
        PassType::REMOVE_TRIVIAL_RZ,
        PassType::SYNTHESIZE_RZ,
        PassType::TO_PBC
    };
    
    // Convert to PBC format with T-count optimization
    inline const std::vector<PassType> TO_PBC_OPTIMIZED = {
        PassType::DECOMPOSE,
        PassType::REMOVE_TRIVIAL_RZ,
        PassType::SYNTHESIZE_RZ,
        PassType::TO_PBC,
        PassType::TFUSE
    };
    
    // Apply Clifford reduction optimization (from FTQC codesign paper)
    inline const std::vector<PassType> TO_CLIFFORD_REDUCTION = {
        PassType::DECOMPOSE,
        PassType::REMOVE_TRIVIAL_RZ,
        PassType::SYNTHESIZE_RZ,
        PassType::CLIFFORD_REDUCTION
    };
    
    // === INCREMENTAL PASSES (for composing workflows) ===
    
    // Just T-optimization (assumes input is already PBC)
    inline const std::vector<PassType> T_OPTIMIZATION_ONLY = {
        PassType::TFUSE
    };
    
    // Post-processing cleanup
    inline const std::vector<PassType> CLEANUP = {
        PassType::GATE_FUSION,
        PassType::REMOVE_TRIVIAL_RZ
    };
}

} // namespace NWQEC