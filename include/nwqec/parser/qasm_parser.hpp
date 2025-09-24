#pragma once

#include "nwqec/parser/token.hpp"
#include "nwqec/parser/lexer.hpp"
#include "nwqec/parser/ast.hpp"
#include "nwqec/parser/ast_generator.hpp"
#include "nwqec/parser/ast_converter.hpp"

#include "nwqec/core/circuit.hpp"

#include <string>
#include <fstream>
#include <sstream>
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
        std::unique_ptr<ASTProgram> program;
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

            try
            {
                // Lexical analysis
                Lexer lexer(source);
                std::vector<Token> tokens = lexer.tokenize();

                // Check for lexer errors
                for (const auto &token : tokens)
                {
                    if (token.type == TokenType::INVALID)
                    {
                        lastError = "Lexical error at line " + std::to_string(token.line) +
                                    ", column " + std::to_string(token.column) +
                                    ": Invalid token '" + token.lexeme + "'";
                        hasError = true;
                        return false;
                    }
                }

                // Parse tokens into AST
                ASTGenerator ast_gen(tokens);
                program = std::make_unique<ASTProgram>(ast_gen.parse());

                // Build the flattened circuit
                ASTCircuitConverter builder;
                circuit = std::make_unique<Circuit>(builder.build(program.get()));

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
            std::ifstream file(filename);
            if (!file.is_open())
            {
                lastError = "Could not open file: " + filename;
                hasError = true;
                return false;
            }

            std::stringstream buffer;
            buffer << file.rdbuf();
            return parse_string(buffer.str());
        }

        /**
         * Get the parsed program (AST representation)
         *
         * @return pointer to the parsed program, or nullptr if parsing failed
         */
        const ASTProgram *get_program() const
        {
            if (hasError || !program)
                return nullptr;
            return program.get();
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
         * Execute the parsed program
         *
         * @return true if execution succeeded, false otherwise
         */
        bool execute()
        {
            if (hasError || !program)
            {
                return false;
            }

            try
            {
                program->execute();
                return true;
            }
            catch (const std::exception &e)
            {
                lastError = "Execution error: " + std::string(e.what());
                hasError = true;
                return false;
            }
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