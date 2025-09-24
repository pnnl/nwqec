#pragma once
#include <iostream>
#include <string>
#include <type_traits>

#include "gmp_integer.hpp"
#include "gmp_float.hpp"

/**
 * types.hpp
 *
 * Defines the standard integer and floating-point types used throughout
 * the gridsynth library. This provides a central location to configure
 * the precision and range of numerical types.
 */

namespace gridsynth
{
    // using Integer = Int128;
    // using Integer = long long;
    using Integer = GMPInteger;

    using Float = GMPFloat;

} // namespace gridsynth
