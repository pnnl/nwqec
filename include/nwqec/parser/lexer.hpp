#pragma once

#include "token.hpp"
#include <string>
#include <vector>
#include <unordered_map>
#include <cctype>
#include <stdexcept>
#include <regex>

namespace NWQEC
{

    /**
     * Lexer class tokenizes a QASM program string into a sequence of tokens
     */
    class Lexer
    {
    private:
        std::string source;
        std::vector<Token> tokens;
        size_t start = 0;
        size_t current = 0;
        int line = 1;
        int column = 1;

        // Keywords map for fast lookup
        const std::unordered_map<std::string, TokenType> keywords = {
            {"OPENQASM", TokenType::OPENQASM},
            {"include", TokenType::INCLUDE},
            {"qreg", TokenType::QREG},
            {"creg", TokenType::CREG},
            {"gate", TokenType::GATE},
            {"measure", TokenType::MEASURE},
            {"reset", TokenType::RESET},
            {"if", TokenType::IF},
            {"barrier", TokenType::BARRIER},
            {"pi", TokenType::PI}};

    public:
        Lexer(const std::string &source) : source(source) {}

        std::vector<Token> tokenize()
        {
            tokens.clear();
            start = 0;
            current = 0;
            line = 1;
            column = 1;

            while (!isAtEnd())
            {
                start = current;
                scanToken();
            }

            // Add EOF token
            tokens.emplace_back(TokenType::EOF_TOKEN, "", line, column);
            return tokens;
        }

    private:
        bool isAtEnd() const
        {
            return current >= source.length();
        }

        char advance()
        {
            column++;
            return source[current++];
        }

        char peek() const
        {
            if (isAtEnd())
                return '\0';
            return source[current];
        }

        char peekNext() const
        {
            if (current + 1 >= source.length())
                return '\0';
            return source[current + 1];
        }

        bool match(char expected)
        {
            if (isAtEnd() || source[current] != expected)
                return false;

            current++;
            column++;
            return true;
        }

        void addToken(TokenType type)
        {
            std::string lexeme = source.substr(start, current - start);
            tokens.emplace_back(type, lexeme, line, column - (current - start));
        }

        void scanToken()
        {
            char c = advance();

            switch (c)
            {
            // Single character tokens
            case '(':
                addToken(TokenType::LPAREN);
                break;
            case ')':
                addToken(TokenType::RPAREN);
                break;
            case '{':
                addToken(TokenType::LBRACE);
                break;
            case '}':
                addToken(TokenType::RBRACE);
                break;
            case '[':
                addToken(TokenType::LBRACKET);
                break;
            case ']':
                addToken(TokenType::RBRACKET);
                break;
            case ',':
                addToken(TokenType::COMMA);
                break;
            case ';':
                addToken(TokenType::SEMICOLON);
                break;
            case '+':
                addToken(TokenType::PLUS);
                break;
            case '-':
                // Check if it's the arrow operator '->'
                if (match('>'))
                {
                    addToken(TokenType::ARROW);
                }
                else
                {
                    addToken(TokenType::MINUS);
                }
                break;
            case '*':
                addToken(TokenType::TIMES);
                break;
            case '/':
                if (match('/'))
                {
                    // Comment until end of line
                    while (peek() != '\n' && !isAtEnd())
                        advance();
                }
                else
                {
                    addToken(TokenType::DIVIDE);
                }
                break;
            case '^':
                addToken(TokenType::POWER);
                break;

            // Two-character tokens
            case '=':
                if (match('='))
                {
                    addToken(TokenType::EQUALS);
                }
                else
                {
                    addToken(TokenType::ASSIGN);
                }
                break;

            // String literals
            case '"':
                string();
                break;

            // Whitespace
            case ' ':
            case '\r':
            case '\t':
                // Ignore whitespace
                break;
            case '\n':
                line++;
                column = 1;
                break;

            default:
                if (std::isdigit(c))
                {
                    // Number literal
                    number();
                }
                else if (std::isalpha(c) || c == '_')
                {
                    // Identifier or keyword
                    identifier();
                }
                else
                {
                    // Invalid character
                    tokens.emplace_back(TokenType::INVALID, std::string(1, c), line, column - 1);
                }
                break;
            }
        }

        void identifier()
        {
            while (std::isalnum(peek()) || peek() == '_')
                advance();

            std::string text = source.substr(start, current - start);

            // Check if it's a keyword
            auto it = keywords.find(text);
            TokenType type = (it != keywords.end()) ? it->second : TokenType::IDENTIFIER;

            addToken(type);
        }

        void number()
        {
            bool isReal = false;

            // Integer part
            while (std::isdigit(peek()))
                advance();

            // Look for decimal point
            if (peek() == '.' && std::isdigit(peekNext()))
            {
                isReal = true;

                // Consume the '.'
                advance();

                // Fractional part
                while (std::isdigit(peek()))
                    advance();
            }

            // Look for scientific notation (e.g., 1e10, 1.5e-7)
            if ((peek() == 'e' || peek() == 'E'))
            {
                isReal = true;

                // Consume the 'e' or 'E'
                advance();

                // Check for +/- sign
                if (peek() == '+' || peek() == '-')
                    advance();

                // Must have at least one digit after 'e'
                if (!std::isdigit(peek()))
                {
                    // Error: invalid scientific notation
                    tokens.emplace_back(TokenType::INVALID, source.substr(start, current - start),
                                        line, column - (current - start));
                    return;
                }

                // Consume digits after 'e'
                while (std::isdigit(peek()))
                    advance();
            }

            addToken(isReal ? TokenType::REAL : TokenType::INTEGER);
        }

        void string()
        {
            // Consume characters until closing quote
            while (peek() != '"' && !isAtEnd())
            {
                if (peek() == '\n')
                {
                    line++;
                    column = 1;
                }
                advance();
            }

            // Unterminated string
            if (isAtEnd())
            {
                tokens.emplace_back(TokenType::INVALID, source.substr(start, current - start),
                                    line, column - (current - start));
                return;
            }

            // Consume the closing quote
            advance();

            // Extract the string value (excluding quotes)
            std::string value = source.substr(start + 1, current - start - 2);
            tokens.emplace_back(TokenType::STRING, value, line, column - (current - start));
        }
    };

} // namespace NWQEC