#pragma once

#include "token.hpp"
#include "ast.hpp"
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <limits>

namespace NWQEC
{

    class ParseError : public std::runtime_error
    {
    public:
        ParseError(const std::string &message, int line, int column)
            : std::runtime_error("Parse error at line " + std::to_string(line) +
                                 ", column " + std::to_string(column) + ": " + message) {}
    };

    /**
     * Parser class that converts tokens into an AST
     */
    class ASTGenerator
    {
    private:
        std::vector<Token> tokens;
        size_t current = 0;

        // Maps for fast checking if a given identifier is a valid gate or not
        std::unordered_map<std::string, bool> predefined_gates = {
            {"x", true},
            {"y", true},
            {"z", true},
            {"h", true},
            {"s", true},
            {"sdg", true},
            {"t", true},
            {"tdg", true},
            {"rx", true},
            {"ry", true},
            {"rz", true},
            {"sx", true},
            {"sxdg", true},
            {"rxp4", true},
            {"rxp4dg", true},
            {"p", true},
            {"u", true},
            {"u1", true},
            {"u2", true},
            {"u3", true},
            {"cx", true},
            {"cy", true},
            {"cz", true},
            {"ch", true},
            {"cs", true},
            {"csdg", true},
            {"ct", true},
            {"ctdg", true},
            {"ecr", true},
            {"crx", true},
            {"cry", true},
            {"crz", true},
            {"csx", true},
            {"cp", true},
            {"cu", true},
            {"cu1", true},
            {"cu3", true},
            {"rxx", true},
            {"ryy", true},
            {"rzz", true},
            {"id", true},
            {"swap", true},
            {"ccx", true},
            {"cswap", true},
            {"rccx", true},
            {"t_pauli", true},
            {"s_pauli", true},
            {"z_pauli", true},
            {"m_pauli", true}};

        // User-defined gates (would be populated as we parse gate declarations)
        std::unordered_map<std::string, bool> user_defined_gates;

    public:
        ASTGenerator(const std::vector<Token> &tokens) : tokens(tokens) {}

        ASTProgram parse()
        {
            ASTProgram program;

            while (!is_at_end())
            {
                try
                {
                    program.add_statement(declaration());
                }
                catch (const ParseError &e)
                {
                    // Report error and synchronize
                    std::cerr << e.what() << std::endl;
                    synchronize();
                }
            }

            return program;
        }

    private:
        //------------------------
        // Helper methods
        //------------------------

        bool is_at_end() const
        {
            return peek().type == TokenType::EOF_TOKEN;
        }

        const Token &peek() const
        {
            return tokens[current];
        }

        const Token &previous() const
        {
            return tokens[current - 1];
        }

        const Token &advance()
        {
            if (!is_at_end())
                current++;
            return previous();
        }

        bool check(TokenType type) const
        {
            if (is_at_end())
                return false;
            return peek().type == type;
        }

        bool match(std::initializer_list<TokenType> types)
        {
            for (TokenType type : types)
            {
                if (check(type))
                {
                    advance();
                    return true;
                }
            }
            return false;
        }

        const Token &consume(TokenType type, const std::string &message)
        {
            if (check(type))
                return advance();

            throw ParseError(message, peek().line, peek().column);
        }

        void synchronize()
        {
            advance();

            while (!is_at_end())
            {
                if (previous().type == TokenType::SEMICOLON)
                    return;

                switch (peek().type)
                {
                case TokenType::OPENQASM:
                case TokenType::INCLUDE:
                case TokenType::QREG:
                case TokenType::CREG:
                case TokenType::GATE:
                case TokenType::MEASURE:
                case TokenType::RESET:
                case TokenType::BARRIER:
                case TokenType::IF:
                    return;
                default:
                    break;
                }

                advance();
            }
        }

        //------------------------
        // Parsing methods
        //------------------------

        std::unique_ptr<Stmt> declaration()
        {
            if (match({TokenType::OPENQASM}))
                return version_declaration();
            if (match({TokenType::INCLUDE}))
                return include_statement();
            if (match({TokenType::QREG}))
                return qreg_declaration();
            if (match({TokenType::CREG}))
                return creg_declaration();
            if (match({TokenType::GATE}))
                return gate_declaration();

            return statement();
        }

        std::unique_ptr<Stmt> version_declaration()
        {
            Token version_token = consume(TokenType::REAL, "Expected version number after OPENQASM.");
            std::string version = version_token.lexeme;
            consume(TokenType::SEMICOLON, "Expected ';' after version number.");

            return std::make_unique<VersionDecl>(version);
        }

        std::unique_ptr<Stmt> include_statement()
        {
            Token filename_token = consume(TokenType::STRING, "Expected file name after include.");
            std::string filename = filename_token.lexeme;
            consume(TokenType::SEMICOLON, "Expected ';' after include statement.");

            return std::make_unique<IncludeStmt>(filename);
        }

        size_t parse_register_size(const Token &size_token)
        {
            try
            {
                unsigned long long raw = std::stoull(size_token.lexeme);
                if (raw > std::numeric_limits<size_t>::max())
                {
                    throw ParseError("Register size " + size_token.lexeme + " exceeds supported range.",
                                     size_token.line, size_token.column);
                }
                return static_cast<size_t>(raw);
            }
            catch (const std::invalid_argument &)
            {
                throw ParseError("Invalid register size: " + size_token.lexeme + ".",
                                 size_token.line, size_token.column);
            }
            catch (const std::out_of_range &)
            {
                throw ParseError("Register size " + size_token.lexeme + " is out of range.",
                                 size_token.line, size_token.column);
            }
        }

        long long parse_integer_literal(const Token &token)
        {
            try
            {
                return std::stoll(token.lexeme);
            }
            catch (const std::invalid_argument &)
            {
                throw ParseError("Invalid integer literal: " + token.lexeme + ".",
                                 token.line, token.column);
            }
            catch (const std::out_of_range &)
            {
                throw ParseError("Integer literal " + token.lexeme + " is out of range.",
                                 token.line, token.column);
            }
        }

        std::unique_ptr<Stmt> qreg_declaration()
        {
            Token name_token = consume(TokenType::IDENTIFIER, "Expected register name after qreg.");
            std::string name = name_token.lexeme;

            consume(TokenType::LBRACKET, "Expected '[' after register name.");
            Token size_token = consume(TokenType::INTEGER, "Expected size after '['.");
            size_t size = parse_register_size(size_token);
            consume(TokenType::RBRACKET, "Expected ']' after size.");
            consume(TokenType::SEMICOLON, "Expected ';' after register declaration.");

            return std::make_unique<QRegDecl>(name, size);
        }

        std::unique_ptr<Stmt> creg_declaration()
        {
            Token name_token = consume(TokenType::IDENTIFIER, "Expected register name after creg.");
            std::string name = name_token.lexeme;

            consume(TokenType::LBRACKET, "Expected '[' after register name.");
            Token size_token = consume(TokenType::INTEGER, "Expected size after '['.");
            size_t size = parse_register_size(size_token);
            consume(TokenType::RBRACKET, "Expected ']' after size.");
            consume(TokenType::SEMICOLON, "Expected ';' after register declaration.");

            return std::make_unique<CRegDecl>(name, size);
        }

        std::unique_ptr<Stmt> gate_declaration()
        {
            // Parse gate name
            Token name_token = consume(TokenType::IDENTIFIER, "Expected gate name.");
            std::string name = name_token.lexeme;

            // Register this gate name so we can use it later
            user_defined_gates[name] = true;

            std::vector<std::string> params;
            std::vector<std::string> qubits;

            // Parse parameters (optional)
            if (match({TokenType::LPAREN}))
            {
                if (!check(TokenType::RPAREN))
                {
                    do
                    {
                        Token param_token = consume(TokenType::IDENTIFIER, "Expected parameter name.");
                        params.push_back(param_token.lexeme);
                    } while (match({TokenType::COMMA}));
                }
                consume(TokenType::RPAREN, "Expected ')' after parameters.");
            }

            // Parse qubits
            if (!check(TokenType::LBRACE))
            {
                do
                {
                    Token qubit_token = consume(TokenType::IDENTIFIER, "Expected qubit name.");
                    qubits.push_back(qubit_token.lexeme);
                } while (match({TokenType::COMMA}));
            }

            // Parse gate body
            std::vector<std::unique_ptr<Stmt>> body = gate_body();

            return std::make_unique<GateDeclStmt>(name, params, qubits, std::move(body));
        }

        std::vector<std::unique_ptr<Stmt>> gate_body()
        {
            consume(TokenType::LBRACE, "Expected '{' before gate body.");

            std::vector<std::unique_ptr<Stmt>> statements;
            while (!check(TokenType::RBRACE) && !is_at_end())
            {
                statements.push_back(gate_operation());
            }

            consume(TokenType::RBRACE, "Expected '}' after gate body.");
            return statements;
        }

        std::unique_ptr<Stmt> gate_operation()
        {
            // Only certain statements are allowed in gate body
            Token token = peek();

            if (match({TokenType::IDENTIFIER}))
            {
                // This is a gate application inside a gate definition
                std::string gate_name = previous().lexeme;

                // Parse parameters if any
                std::vector<std::unique_ptr<Expr>> params;
                if (match({TokenType::LPAREN}))
                {
                    if (!check(TokenType::RPAREN))
                    {
                        params.push_back(expression());
                        while (match({TokenType::COMMA}))
                        {
                            params.push_back(expression());
                        }
                    }
                    consume(TokenType::RPAREN, "Expected ')' after gate parameters.");
                }

                // Parse qubits
                std::vector<std::unique_ptr<Expr>> qubits;
                qubits.push_back(std::make_unique<VariableExpr>(consume(
                                                                    TokenType::IDENTIFIER, "Expected qubit argument.")
                                                                    .lexeme));

                while (match({TokenType::COMMA}))
                {
                    qubits.push_back(std::make_unique<VariableExpr>(consume(
                                                                        TokenType::IDENTIFIER, "Expected qubit argument.")
                                                                        .lexeme));
                }

                consume(TokenType::SEMICOLON, "Expected ';' after gate operation.");
                return std::make_unique<GateStmt>(gate_name, std::move(params), std::move(qubits));
            }

            throw ParseError("Expected gate operation in gate body.", token.line, token.column);
        }

        std::unique_ptr<Stmt> statement()
        {
            if (match({TokenType::MEASURE}))
                return measure_statement();
            if (match({TokenType::RESET}))
                return reset_statement();
            if (match({TokenType::IF}))
                return if_statement();
            if (match({TokenType::LBRACE}))
                return block_statement();
            if (match({TokenType::BARRIER}))
                return barrier_statement();

            // Must be a gate application
            return gate_statement();
        }

        std::unique_ptr<Stmt> barrier_statement()
        {
            std::vector<std::unique_ptr<Expr>> qubits;
            do
            {
                qubits.push_back(primary_expr());
            } while (match({TokenType::COMMA}));

            consume(TokenType::SEMICOLON, "Expected ';' after barrier statement.");

            return std::make_unique<BarrierStmt>(std::move(qubits));
        }

        std::unique_ptr<Stmt> measure_statement()
        {
            std::unique_ptr<Expr> qubit = primary_expr();

            // Look for -> operator (now a single ARROW token)
            consume(TokenType::ARROW, "Expected '->' after qubit in measure statement.");

            std::unique_ptr<Expr> bit = primary_expr();
            consume(TokenType::SEMICOLON, "Expected ';' after measure statement.");

            return std::make_unique<MeasureStmt>(std::move(qubit), std::move(bit));
        }

        std::unique_ptr<Stmt> reset_statement()
        {
            std::unique_ptr<Expr> qubit = primary_expr();
            consume(TokenType::SEMICOLON, "Expected ';' after reset statement.");

            return std::make_unique<ResetStmt>(std::move(qubit));
        }

        std::unique_ptr<Stmt> if_statement()
        {
            consume(TokenType::LPAREN, "Expected '(' after if.");
            Token creg_token = consume(TokenType::IDENTIFIER, "Expected classical register name.");
            consume(TokenType::EQUALS, "Expected '==' after register name.");
            Token value_token = consume(TokenType::INTEGER, "Expected integer after '=='.");
            consume(TokenType::RPAREN, "Expected ')' after condition.");

            std::unique_ptr<Stmt> then_branch = statement();

            return std::make_unique<IfStmt>(creg_token.lexeme, parse_integer_literal(value_token), std::move(then_branch));
        }

        std::unique_ptr<Stmt> block_statement()
        {
            std::vector<std::unique_ptr<Stmt>> statements;

            while (!check(TokenType::RBRACE) && !is_at_end())
            {
                statements.push_back(declaration());
            }

            consume(TokenType::RBRACE, "Expected '}' after block.");

            return std::make_unique<BlockStmt>(std::move(statements));
        }

        std::unique_ptr<Stmt> gate_statement()
        {
            Token name_token = consume(TokenType::IDENTIFIER, "Expected gate name.");
            std::string name = name_token.lexeme;

            // Check if this is a valid gate
            bool is_valid_gate = predefined_gates.find(name) != predefined_gates.end() ||
                                 user_defined_gates.find(name) != user_defined_gates.end();

            if (!is_valid_gate)
            {
                throw ParseError("Unknown gate: " + name, name_token.line, name_token.column);
            }

            // Special handling for Pauli gates (t_pauli, m_pauli, s_pauli)
            if (name == "t_pauli" || name == "s_pauli" || name == "z_pauli" || name == "m_pauli")
            {
                // Expect: GATE_NAME (+|-)??PAULISTRING; (sign optional, defaults to '+')
                std::string pauli_string;
                if (match({TokenType::PLUS}))
                {
                    pauli_string.push_back('+');
                }
                else if (match({TokenType::MINUS}))
                {
                    pauli_string.push_back('-');
                }
                else
                {
                    pauli_string.push_back('+');
                }

                Token pauli_token = consume(TokenType::IDENTIFIER, "Expected Pauli string after " + name + ": e.g., +XYZI");
                pauli_string += pauli_token.lexeme;

                consume(TokenType::SEMICOLON, "Expected ';' after Pauli gate.");

                return std::make_unique<PauliStmt>(name, pauli_string);
            }

            // Parse parameters if any
            std::vector<std::unique_ptr<Expr>> params;
            if (match({TokenType::LPAREN}))
            {
                if (!check(TokenType::RPAREN))
                {
                    params.push_back(expression());
                    while (match({TokenType::COMMA}))
                    {
                        params.push_back(expression());
                    }
                }
                consume(TokenType::RPAREN, "Expected ')' after gate parameters.");
            }

            // Parse qubits
            std::vector<std::unique_ptr<Expr>> qubits;
            qubits.push_back(primary_expr());

            while (match({TokenType::COMMA}))
            {
                qubits.push_back(primary_expr());
            }

            consume(TokenType::SEMICOLON, "Expected ';' after gate application.");

            return std::make_unique<GateStmt>(name, std::move(params), std::move(qubits));
        }

        //------------------------
        // Expression parsing
        //------------------------

        std::unique_ptr<Expr> expression()
        {
            return additive_expr();
        }

        std::unique_ptr<Expr> additive_expr()
        {
            std::unique_ptr<Expr> expr = multiplicative_expr();

            while (match({TokenType::PLUS, TokenType::MINUS}))
            {
                Token op = previous();
                std::unique_ptr<Expr> right = multiplicative_expr();

                BinaryExpr::Op opType;
                if (op.type == TokenType::PLUS)
                {
                    opType = BinaryExpr::Op::PLUS;
                }
                else
                {
                    opType = BinaryExpr::Op::MINUS;
                }

                expr = std::make_unique<BinaryExpr>(std::move(expr), opType, std::move(right));
            }

            return expr;
        }

        std::unique_ptr<Expr> multiplicative_expr()
        {
            std::unique_ptr<Expr> expr = unary_expr();

            while (match({TokenType::TIMES, TokenType::DIVIDE}))
            {
                Token op = previous();
                std::unique_ptr<Expr> right = unary_expr();

                BinaryExpr::Op opType;
                if (op.type == TokenType::TIMES)
                {
                    opType = BinaryExpr::Op::MULTIPLY;
                }
                else
                {
                    opType = BinaryExpr::Op::DIVIDE;
                }

                expr = std::make_unique<BinaryExpr>(std::move(expr), opType, std::move(right));
            }

            return expr;
        }

        std::unique_ptr<Expr> unary_expr()
        {
            if (match({TokenType::MINUS}))
            {
                std::unique_ptr<Expr> right = unary_expr();
                // Represent -x as 0 - x
                auto zero = std::make_unique<NumberExpr>(0.0, true);
                return std::make_unique<BinaryExpr>(std::move(zero), BinaryExpr::Op::MINUS, std::move(right));
            }

            return power_expr();
        }

        std::unique_ptr<Expr> power_expr()
        {
            std::unique_ptr<Expr> expr = primary_expr();

            if (match({TokenType::POWER}))
            {
                std::unique_ptr<Expr> right = unary_expr();
                expr = std::make_unique<BinaryExpr>(std::move(expr), BinaryExpr::Op::POWER, std::move(right));
            }

            return expr;
        }

        std::unique_ptr<Expr> primary_expr()
        {
            if (match({TokenType::INTEGER}))
            {
                long long value = parse_integer_literal(previous());
                return std::make_unique<NumberExpr>(static_cast<double>(value), true);
            }

            if (match({TokenType::REAL}))
            {
                double value = std::stod(previous().lexeme);
                return std::make_unique<NumberExpr>(value, false);
            }

            if (match({TokenType::PI}))
            {
                return std::make_unique<PiExpr>();
            }

            if (match({TokenType::LPAREN}))
            {
                std::unique_ptr<Expr> expr = expression();
                consume(TokenType::RPAREN, "Expected ')' after expression.");
                return std::make_unique<ParenExpr>(std::move(expr));
            }

            if (match({TokenType::IDENTIFIER}))
            {
                std::string name = previous().lexeme;

                // Check if it's a register index
                if (match({TokenType::LBRACKET}))
                {
                    std::unique_ptr<Expr> index = expression();
                    consume(TokenType::RBRACKET, "Expected ']' after index.");
                    return std::make_unique<IndexExpr>(name, std::move(index));
                }

                // Otherwise it's a variable reference
                return std::make_unique<VariableExpr>(name);
            }

            throw ParseError("Expected expression.", peek().line, peek().column);
        }

        // Helper to peek one token ahead
        const Token &peekNext() const
        {
            if (current + 1 >= tokens.size())
                return tokens.back();
            return tokens[current + 1];
        }
    };

} // namespace NWQEC
