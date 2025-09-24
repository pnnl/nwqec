// Common constants shared across the project
#pragma once

namespace NWQEC
{
    // Default epsilon selection: two orders of magnitude smaller than theta
    inline constexpr double DEFAULT_EPSILON_MULTIPLIER = 1e-2;

    // Default numerical precision for mpmath (Python fallback)
    inline constexpr int DEFAULT_MPMATH_DPS = 128;

    // Default timeouts for gridsynth diophantine and factoring (milliseconds)
    inline constexpr int DEFAULT_DIOPHANTINE_TIMEOUT_MS = 200;
    inline constexpr int DEFAULT_FACTORING_TIMEOUT_MS = 50;

}
