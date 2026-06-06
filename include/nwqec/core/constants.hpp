// Common constants shared across the project
#pragma once

namespace NWQEC
{
    enum class RzErrorPolicy
    {
        PER_GATE,
        TOTAL,
        RELATIVE
    };

    inline constexpr RzErrorPolicy DEFAULT_RZ_ERROR_POLICY = RzErrorPolicy::PER_GATE;

    // Default fixed absolute synthesis tolerance per RZ gate.
    inline constexpr double DEFAULT_RZ_PER_GATE_EPSILON = 1e-10;

    // Default total synthesis error budget. Split evenly over RZ occurrences.
    inline constexpr double DEFAULT_RZ_TOTAL_EPSILON = 1e-2;

    // Default relative synthesis tolerance: epsilon(theta) = |theta| * factor.
    inline constexpr double DEFAULT_RZ_RELATIVE_EPSILON = 1e-2;

    // Retained for standalone gridsynth's angle-relative default.
    inline constexpr double DEFAULT_EPSILON_MULTIPLIER = DEFAULT_RZ_RELATIVE_EPSILON;

    inline constexpr double default_rz_error_epsilon(RzErrorPolicy policy)
    {
        switch (policy)
        {
        case RzErrorPolicy::PER_GATE:
            return DEFAULT_RZ_PER_GATE_EPSILON;
        case RzErrorPolicy::TOTAL:
            return DEFAULT_RZ_TOTAL_EPSILON;
        case RzErrorPolicy::RELATIVE:
            return DEFAULT_RZ_RELATIVE_EPSILON;
        default:
            return DEFAULT_RZ_PER_GATE_EPSILON;
        }
    }

    // Default numerical precision for mpmath (Python fallback)
    inline constexpr int DEFAULT_MPMATH_DPS = 128;

    // Default timeouts for gridsynth diophantine and factoring (milliseconds)
    inline constexpr int DEFAULT_DIOPHANTINE_TIMEOUT_MS = 200;
    inline constexpr int DEFAULT_FACTORING_TIMEOUT_MS = 50;

}
