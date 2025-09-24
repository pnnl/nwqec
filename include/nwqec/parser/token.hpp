#pragma once

#include <string>
#include <variant>
#include <vector>
#include <unordered_map>
#include <iostream>

namespace NWQEC
{

    /**
     * TokenType defines all possible types of tokens in QASM code
     */
    enum class TokenType
    {
        // Keywords
        OPENQASM,
        INCLUDE,
        QREG,
        CREG,
        GATE,
        MEASURE,
        RESET,
        IF,
        BARRIER,

        // Identifiers and literals
        IDENTIFIER,
        REAL,
        INTEGER,
        STRING,

        // Operators
        PLUS,   // +
        MINUS,  // -
        TIMES,  // *
        DIVIDE, // /
        POWER,  // ^
        ASSIGN, // =
        EQUALS, // ==
        ARROW,  // ->

        // Delimiters
        SEMICOLON, // ;
        COMMA,     // ,
        LPAREN,    // (
        RPAREN,    // )
        LBRACE,    // {
        RBRACE,    // }
        LBRACKET,  // [
        RBRACKET,  // ]

        // Special
        PI, // pi
        EOF_TOKEN,
        INVALID
    };

    /**
     * Token represents a lexical token in QASM code
     */
    struct Token
    {
        TokenType type;
        std::string lexeme;
        int line;
        int column;

        Token(TokenType type, std::string lexeme, int line, int column)
            : type(type), lexeme(std::move(lexeme)), line(line), column(column) {}

        std::string to_string() const
        {
            // Print token information for debugging
            std::string typeStr;

            switch (type)
            {
            case TokenType::OPENQASM:
                typeStr = "OPENQASM";
                break;
            case TokenType::INCLUDE:
                typeStr = "INCLUDE";
                break;
            case TokenType::QREG:
                typeStr = "QREG";
                break;
            case TokenType::CREG:
                typeStr = "CREG";
                break;
            case TokenType::GATE:
                typeStr = "GATE";
                break;
            case TokenType::MEASURE:
                typeStr = "MEASURE";
                break;
            case TokenType::RESET:
                typeStr = "RESET";
                break;
            case TokenType::IF:
                typeStr = "IF";
                break;
            case TokenType::BARRIER:
                typeStr = "BARRIER";
                break;
            case TokenType::IDENTIFIER:
                typeStr = "IDENTIFIER";
                break;
            case TokenType::REAL:
                typeStr = "REAL";
                break;
            case TokenType::INTEGER:
                typeStr = "INTEGER";
                break;
            case TokenType::STRING:
                typeStr = "STRING";
                break;
            case TokenType::PLUS:
                typeStr = "PLUS";
                break;
            case TokenType::MINUS:
                typeStr = "MINUS";
                break;
            case TokenType::TIMES:
                typeStr = "TIMES";
                break;
            case TokenType::DIVIDE:
                typeStr = "DIVIDE";
                break;
            case TokenType::POWER:
                typeStr = "POWER";
                break;
            case TokenType::ASSIGN:
                typeStr = "ASSIGN";
                break;
            case TokenType::EQUALS:
                typeStr = "EQUALS";
                break;
            case TokenType::ARROW:
                typeStr = "ARROW";
                break;
            case TokenType::SEMICOLON:
                typeStr = "SEMICOLON";
                break;
            case TokenType::COMMA:
                typeStr = "COMMA";
                break;
            case TokenType::LPAREN:
                typeStr = "LPAREN";
                break;
            case TokenType::RPAREN:
                typeStr = "RPAREN";
                break;
            case TokenType::LBRACE:
                typeStr = "LBRACE";
                break;
            case TokenType::RBRACE:
                typeStr = "RBRACE";
                break;
            case TokenType::LBRACKET:
                typeStr = "LBRACKET";
                break;
            case TokenType::RBRACKET:
                typeStr = "RBRACKET";
                break;
            case TokenType::PI:
                typeStr = "PI";
                break;
            case TokenType::EOF_TOKEN:
                typeStr = "EOF";
                break;
            case TokenType::INVALID:
                typeStr = "INVALID";
                break;
            default:
                typeStr = "UNKNOWN";
                break;
            }

            return "Token(" + typeStr + ", '" + lexeme + "', line=" +
                   std::to_string(line) + ", col=" + std::to_string(column) + ")";
        }
    };

} // namespace NWQEC