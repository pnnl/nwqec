#pragma once

// Core components
#include "nwqec/core/circuit.hpp"

#include <vector>
#include <memory>
#include <string>
#include <set>


namespace NWQEC
{
    /**
     * @brief Base class for all circuit transformation passes
     *
     * A pass represents a transformation that can be applied to a quantum circuit.
     */
    class Pass
    {
    public:
        /**
         * @brief Constructor for Pass
         */
        Pass() = default;

        virtual ~Pass() = default;

        /**
         * @brief Run the pass on the given circuit
         *
         * @param circuit The circuit to transform
         * @return true if the circuit was modified, false otherwise
         */
        virtual bool run(Circuit &circuit) = 0;

        /**
         * @brief Get the name of the pass
         *
         * @return The pass name
         */
        virtual std::string get_name() const = 0;
    };

} // namespace NWQEC