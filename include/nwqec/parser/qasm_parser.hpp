#pragma once

#include "nwqec/parser/circuit_parser.hpp"

#include "nwqec/core/circuit.hpp"

#include <string>
#include <fstream>
#include <iostream>
#include <memory>

namespace NWQEC
{

    /**
     * Main API class for parsing QASM code
     */
    class QASMParser
    {
    private:
        std::unique_ptr<Circuit> circuit;
        std::string lastError;
        bool hasError = false;

    public:
        /**
         * Parse QASM code from a string
         *
         * @param source QASM code string to parse
         * @return true if parsing succeeded, false otherwise
         */
        bool parse_string(const std::string &source)
        {
            hasError = false;
            lastError = "";
            circuit.reset();

            try
            {
                // Build the flattened circuit directly from tokens.
                CircuitParser circuit_parser(source);
                circuit = std::make_unique<Circuit>(circuit_parser.parse());

                return true;
            }
            catch (const std::exception &e)
            {
                lastError = e.what();
                hasError = true;
                return false;
            }
        }

        /**
         * Parse QASM code from a file
         *
         * @param filename Path to the QASM file
         * @return true if parsing succeeded, false otherwise
         */
        bool parse_file(const std::string &filename)
        {
            std::ifstream file(filename, std::ios::binary | std::ios::ate);
            if (!file.is_open())
            {
                lastError = "Could not open file: " + filename;
                hasError = true;
                return false;
            }

            const std::streamsize size = file.tellg();
            if (size < 0)
            {
                lastError = "Could not determine file size: " + filename;
                hasError = true;
                return false;
            }

            std::string source(static_cast<size_t>(size), '\0');
            file.seekg(0, std::ios::beg);
            if (size > 0 && !file.read(source.data(), size))
            {
                lastError = "Could not read file: " + filename;
                hasError = true;
                return false;
            }

            return parse_string(source);
        }

        /**
         * Get a mutable copy of the circuit
         *
         * @return A unique_ptr to a copy of the circuit that can be modified, or nullptr if parsing failed
         */
        std::unique_ptr<Circuit> get_circuit()
        {
            if (hasError || !circuit)
                return nullptr;
            return std::move(circuit);
        }

        /**
         * Check if the last parsing operation encountered an error
         *
         * @return true if an error occurred, false otherwise
         */
        bool has_parse_error() const
        {
            return hasError;
        }

        /**
         * Get the error message from the last parsing operation
         *
         * @return error message, or empty string if no error occurred
         */
        const std::string &get_error_message() const
        {
            return lastError;
        }

        /**
         * Print the flattened circuit to the specified output stream
         */
        void print_circuit(std::ostream &os) const
        {
            if (!hasError && circuit)
            {
                circuit->print(os);
            }
        }
    };

} // namespace NWQEC
