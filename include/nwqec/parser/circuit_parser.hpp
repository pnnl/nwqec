#pragma once

#include "nwqec/core/circuit.hpp"
#include "nwqec/parser/parse_error.hpp"
#include "nwqec/parser/token.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cctype>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace NWQEC
{
    class SourceLexer
    {
    public:
        struct State
        {
            size_t start = 0;
            size_t current = 0;
            int line = 1;
            int column = 1;
        };

    private:
        const std::string &source;
        State state;

    public:
        explicit SourceLexer(const std::string &source) : source(source) {}

        State get_state() const { return state; }
        void set_state(const State &new_state) { state = new_state; }

        Token next_token()
        {
            skip_ignored();
            state.start = state.current;

            if (is_at_end())
                return Token(TokenType::EOF_TOKEN, "", state.line, state.column);

            char c = advance();
            switch (c)
            {
            case '(':
                return make_token(TokenType::LPAREN);
            case ')':
                return make_token(TokenType::RPAREN);
            case '{':
                return make_token(TokenType::LBRACE);
            case '}':
                return make_token(TokenType::RBRACE);
            case '[':
                return make_token(TokenType::LBRACKET);
            case ']':
                return make_token(TokenType::RBRACKET);
            case ',':
                return make_token(TokenType::COMMA);
            case ';':
                return make_token(TokenType::SEMICOLON);
            case '+':
                return make_token(TokenType::PLUS);
            case '-':
                if (match('>'))
                    return make_token(TokenType::ARROW);
                return make_token(TokenType::MINUS);
            case '*':
                return make_token(TokenType::TIMES);
            case '/':
                return make_token(TokenType::DIVIDE);
            case '^':
                return make_token(TokenType::POWER);
            case '=':
                if (match('='))
                    return make_token(TokenType::EQUALS);
                return make_token(TokenType::ASSIGN);
            case '"':
                return string_token();
            default:
                if (std::isdigit(static_cast<unsigned char>(c)))
                    return number_token();
                if (std::isalpha(static_cast<unsigned char>(c)) || c == '_')
                    return identifier_token();
                return make_token(TokenType::INVALID);
            }
        }

    private:
        bool is_at_end() const { return state.current >= source.size(); }

        char advance()
        {
            ++state.column;
            return source[state.current++];
        }

        char peek() const
        {
            if (is_at_end())
                return '\0';
            return source[state.current];
        }

        char peek_next() const
        {
            if (state.current + 1 >= source.size())
                return '\0';
            return source[state.current + 1];
        }

        bool match(char expected)
        {
            if (is_at_end() || source[state.current] != expected)
                return false;
            ++state.current;
            ++state.column;
            return true;
        }

        void skip_ignored()
        {
            for (;;)
            {
                if (is_at_end())
                    return;

                char c = peek();
                if (c == ' ' || c == '\r' || c == '\t')
                {
                    advance();
                    continue;
                }
                if (c == '\n')
                {
                    advance();
                    ++state.line;
                    state.column = 1;
                    continue;
                }
                if (c == '/' && peek_next() == '/')
                {
                    while (peek() != '\n' && !is_at_end())
                        advance();
                    continue;
                }
                return;
            }
        }

        Token make_token(TokenType type) const
        {
            return Token(type, std::string_view(source.data() + state.start, state.current - state.start),
                         state.line, state.column - static_cast<int>(state.current - state.start));
        }

        Token identifier_token()
        {
            while (std::isalnum(static_cast<unsigned char>(peek())) || peek() == '_')
                advance();

            std::string_view text(source.data() + state.start, state.current - state.start);
            return Token(keyword_type(text), text, state.line,
                         state.column - static_cast<int>(state.current - state.start));
        }

        static TokenType keyword_type(std::string_view text)
        {
            switch (text.size())
            {
            case 2:
                if (text == "if")
                    return TokenType::IF;
                if (text == "pi")
                    return TokenType::PI;
                break;
            case 4:
                if (text == "qreg")
                    return TokenType::QREG;
                if (text == "creg")
                    return TokenType::CREG;
                if (text == "gate")
                    return TokenType::GATE;
                break;
            case 5:
                if (text == "reset")
                    return TokenType::RESET;
                break;
            case 7:
                if (text == "include")
                    return TokenType::INCLUDE;
                if (text == "measure")
                    return TokenType::MEASURE;
                if (text == "barrier")
                    return TokenType::BARRIER;
                break;
            case 8:
                if (text == "OPENQASM")
                    return TokenType::OPENQASM;
                break;
            default:
                break;
            }
            return TokenType::IDENTIFIER;
        }

        Token number_token()
        {
            bool is_real = false;
            while (std::isdigit(static_cast<unsigned char>(peek())))
                advance();

            if (peek() == '.' && std::isdigit(static_cast<unsigned char>(peek_next())))
            {
                is_real = true;
                advance();
                while (std::isdigit(static_cast<unsigned char>(peek())))
                    advance();
            }

            if (peek() == 'e' || peek() == 'E')
            {
                is_real = true;
                advance();
                if (peek() == '+' || peek() == '-')
                    advance();
                if (!std::isdigit(static_cast<unsigned char>(peek())))
                    return make_token(TokenType::INVALID);
                while (std::isdigit(static_cast<unsigned char>(peek())))
                    advance();
            }

            return make_token(is_real ? TokenType::REAL : TokenType::INTEGER);
        }

        Token string_token()
        {
            while (peek() != '"' && !is_at_end())
            {
                if (peek() == '\n')
                {
                    ++state.line;
                    state.column = 1;
                }
                advance();
            }

            if (is_at_end())
                return make_token(TokenType::INVALID);

            advance();
            return Token(TokenType::STRING,
                         std::string_view(source.data() + state.start + 1, state.current - state.start - 2),
                         state.line, state.column - static_cast<int>(state.current - state.start));
        }
    };

    class CircuitParser
    {
    private:
        struct Operand
        {
            std::string name;
            bool has_index = false;
            size_t index = 0;
        };

        struct GateApplication
        {
            std::string name;
            std::vector<double> parameters;
            std::vector<Operand> qubits;
        };

        struct GateDefinition
        {
            std::vector<std::string> qubits;
            std::vector<GateApplication> body;
        };

        const std::string &source;
        SourceLexer lexer;
        Token current_token{TokenType::EOF_TOKEN, "", 1, 1};
        Token previous_token{TokenType::EOF_TOKEN, "", 1, 1};
        Circuit circuit;
        std::unordered_map<std::string_view, size_t> qubit_register_starts;
        std::unordered_map<std::string_view, size_t> qubit_register_sizes;
        std::unordered_map<std::string_view, size_t> bit_register_starts;
        std::unordered_map<std::string_view, size_t> bit_register_sizes;
        std::unordered_map<std::string, GateDefinition> gate_definitions;
        std::vector<std::unordered_map<std::string, std::vector<size_t>>> qubit_binding_stack;

    public:
        explicit CircuitParser(const std::string &source) : source(source), lexer(source) {}

        Circuit parse()
        {
            circuit = Circuit();
            circuit.reserve_operations(count_statements());
            qubit_register_starts.clear();
            qubit_register_sizes.clear();
            bit_register_starts.clear();
            bit_register_sizes.clear();
            gate_definitions.clear();
            qubit_binding_stack.clear();
            current_token = lexer.next_token();

            while (!is_at_end())
            {
                declaration();
            }

            return std::move(circuit);
        }

    private:
        size_t count_statements() const
        {
            size_t count = 0;
            for (char c : source)
            {
                if (c == ';')
                    ++count;
            }
            return count;
        }

        struct ParserState
        {
            SourceLexer::State lexer_state;
            Token current_token;
            Token previous_token;
        };

        ParserState get_state() const
        {
            return {lexer.get_state(), current_token, previous_token};
        }

        void set_state(const ParserState &state)
        {
            lexer.set_state(state.lexer_state);
            current_token = state.current_token;
            previous_token = state.previous_token;
        }

        bool is_at_end() const { return current_token.type == TokenType::EOF_TOKEN; }

        const Token &peek() const { return current_token; }
        const Token &previous() const { return previous_token; }

        const Token &advance()
        {
            if (!is_at_end())
            {
                previous_token = current_token;
                current_token = lexer.next_token();
            }
            return previous_token;
        }

        bool check(TokenType type) const
        {
            return !is_at_end() && peek().type == type;
        }

        bool match(TokenType type)
        {
            if (!check(type))
                return false;
            advance();
            return true;
        }

        const Token &consume(TokenType type, const std::string &message)
        {
            if (check(type))
                return advance();
            throw ParseError(message, peek().line, peek().column);
        }

        static std::string to_string(std::string_view view)
        {
            return std::string(view.data(), view.size());
        }

        static std::string lowercase(std::string_view view)
        {
            std::string out(view.data(), view.size());
            std::transform(out.begin(), out.end(), out.begin(),
                           [](unsigned char c)
                           { return static_cast<char>(std::tolower(c)); });
            return out;
        }

        size_t parse_size(const Token &token, const std::string &kind)
        {
            unsigned long long raw = 0;
            const char *begin = token.lexeme.data();
            const char *end = begin + token.lexeme.size();
            auto result = std::from_chars(begin, end, raw);
            if (result.ec == std::errc::invalid_argument || result.ptr != end)
            {
                throw ParseError("Invalid " + kind + ": " + to_string(token.lexeme) + ".",
                                 token.line, token.column);
            }
            if (result.ec == std::errc::result_out_of_range || raw > std::numeric_limits<size_t>::max())
            {
                throw ParseError(kind + " " + to_string(token.lexeme) + " is out of range.",
                                 token.line, token.column);
            }
            return static_cast<size_t>(raw);
        }

        void declaration()
        {
            if (match(TokenType::OPENQASM))
            {
                if (!match(TokenType::REAL) && !match(TokenType::INTEGER))
                {
                    consume(TokenType::REAL, "Expected version number after OPENQASM.");
                }
                consume(TokenType::SEMICOLON, "Expected ';' after version number.");
                return;
            }
            if (match(TokenType::INCLUDE))
            {
                consume(TokenType::STRING, "Expected file name after include.");
                consume(TokenType::SEMICOLON, "Expected ';' after include statement.");
                return;
            }
            if (match(TokenType::QREG))
            {
                qreg_declaration();
                return;
            }
            if (match(TokenType::CREG))
            {
                creg_declaration();
                return;
            }
            if (match(TokenType::GATE))
            {
                gate_declaration();
                return;
            }
            statement();
        }

        void qreg_declaration()
        {
            std::string_view name_view = consume(TokenType::IDENTIFIER, "Expected register name after qreg.").lexeme;
            std::string name = to_string(name_view);
            consume(TokenType::LBRACKET, "Expected '[' after register name.");
            size_t size = parse_size(consume(TokenType::INTEGER, "Expected size after '['."), "register size");
            consume(TokenType::RBRACKET, "Expected ']' after size.");
            consume(TokenType::SEMICOLON, "Expected ';' after register declaration.");

            qubit_register_starts[name_view] = circuit.get_num_qubits();
            qubit_register_sizes[name_view] = size;
            circuit.add_qreg(name, size);
        }

        void creg_declaration()
        {
            std::string_view name_view = consume(TokenType::IDENTIFIER, "Expected register name after creg.").lexeme;
            std::string name = to_string(name_view);
            consume(TokenType::LBRACKET, "Expected '[' after register name.");
            size_t size = parse_size(consume(TokenType::INTEGER, "Expected size after '['."), "register size");
            consume(TokenType::RBRACKET, "Expected ']' after size.");
            consume(TokenType::SEMICOLON, "Expected ';' after register declaration.");

            bit_register_starts[name_view] = circuit.get_num_bits();
            bit_register_sizes[name_view] = size;
            circuit.add_creg(name, size);
        }

        void gate_declaration()
        {
            std::string name = lowercase(consume(TokenType::IDENTIFIER, "Expected gate name.").lexeme);
            GateDefinition def;

            if (match(TokenType::LPAREN))
            {
                if (!check(TokenType::RPAREN))
                {
                    do
                    {
                        consume(TokenType::IDENTIFIER, "Expected parameter name.");
                    } while (match(TokenType::COMMA));
                }
                consume(TokenType::RPAREN, "Expected ')' after parameters.");
            }

            if (!check(TokenType::LBRACE))
            {
                do
                {
                    def.qubits.push_back(to_string(consume(TokenType::IDENTIFIER, "Expected qubit name.").lexeme));
                } while (match(TokenType::COMMA));
            }

            consume(TokenType::LBRACE, "Expected '{' before gate body.");
            while (!check(TokenType::RBRACE) && !is_at_end())
            {
                def.body.push_back(parse_gate_application(/*allow_indexed_qubits=*/false, /*allow_pauli_statement=*/false));
            }
            consume(TokenType::RBRACE, "Expected '}' after gate body.");

            gate_definitions[name] = std::move(def);
        }

        void statement()
        {
            if (match(TokenType::MEASURE))
            {
                measure_statement();
                return;
            }
            if (match(TokenType::RESET))
            {
                reset_statement();
                return;
            }
            if (match(TokenType::IF))
            {
                throw ParseError("Conditional statements not yet supported in circuit flattening", previous().line, previous().column);
            }
            if (match(TokenType::LBRACE))
            {
                while (!check(TokenType::RBRACE) && !is_at_end())
                {
                    declaration();
                }
                consume(TokenType::RBRACE, "Expected '}' after block.");
                return;
            }
            if (match(TokenType::BARRIER))
            {
                barrier_statement();
                return;
            }

            gate_statement();
        }

        void measure_statement()
        {
            std::vector<size_t> qubits = resolve_quantum_operand(parse_operand(/*allow_index=*/true));
            consume(TokenType::ARROW, "Expected '->' after qubit in measure statement.");
            std::vector<size_t> bits = resolve_classical_operand(parse_operand(/*allow_index=*/true));
            consume(TokenType::SEMICOLON, "Expected ';' after measure statement.");

            size_t count = std::max(qubits.size(), bits.size());
            for (size_t i = 0; i < count; ++i)
            {
                size_t q = qubits.size() > 1 ? qubits[i] : qubits[0];
                size_t b = bits.size() > 1 ? bits[i] : bits[0];
                circuit.add_operation(Operation(Operation::Type::MEASURE, {q}, {}, {b}));
            }
        }

        void reset_statement()
        {
            std::vector<size_t> qubits = resolve_quantum_operand(parse_operand(/*allow_index=*/true));
            consume(TokenType::SEMICOLON, "Expected ';' after reset statement.");
            for (size_t q : qubits)
            {
                circuit.add_operation(Operation(Operation::Type::RESET, {q}));
            }
        }

        void barrier_statement()
        {
            std::vector<size_t> qubits;
            do
            {
                auto sub = resolve_quantum_operand(parse_operand(/*allow_index=*/true));
                qubits.insert(qubits.end(), sub.begin(), sub.end());
            } while (match(TokenType::COMMA));
            consume(TokenType::SEMICOLON, "Expected ';' after barrier statement.");
            circuit.add_operation(Operation(Operation::Type::BARRIER, std::move(qubits)));
        }

        void gate_statement()
        {
            const Token &name_token = consume(TokenType::IDENTIFIER, "Expected gate name.");
            std::string_view name_view = name_token.lexeme;

            if (name_view == "t_pauli" || name_view == "s_pauli" ||
                name_view == "z_pauli" || name_view == "m_pauli")
            {
                parse_pauli_statement(to_string(name_view));
                return;
            }

            Operation::Type op_type;
            std::string lowered_name;
            bool is_builtin = Operation::try_lowercase_name_to_type(name_view, op_type);
            if (!is_builtin)
            {
                lowered_name = lowercase(name_view);
                is_builtin = Operation::try_lowercase_name_to_type(lowered_name, op_type);
            }

            std::vector<double> parameters = parse_optional_parameters();

            if (is_builtin)
            {
                parse_and_apply_builtin_operands(op_type, parameters);
                return;
            }

            GateApplication app;
            app.name = lowered_name.empty() ? to_string(name_view) : std::move(lowered_name);
            app.parameters = std::move(parameters);
            app.qubits.push_back(parse_operand(/*allow_index=*/true));
            while (match(TokenType::COMMA))
            {
                app.qubits.push_back(parse_operand(/*allow_index=*/true));
            }
            consume(TokenType::SEMICOLON, "Expected ';' after gate application.");
            apply_gate_application(app);
        }

        std::vector<double> parse_optional_parameters()
        {
            std::vector<double> parameters;
            if (match(TokenType::LPAREN))
            {
                if (!check(TokenType::RPAREN))
                {
                    parameters.push_back(expression());
                    while (match(TokenType::COMMA))
                    {
                        parameters.push_back(expression());
                    }
                }
                consume(TokenType::RPAREN, "Expected ')' after gate parameters.");
            }
            return parameters;
        }

        void parse_and_apply_builtin_operands(Operation::Type op_type, const std::vector<double> &parameters)
        {
            const ParserState operands_start = get_state();
            if (qubit_binding_stack.empty())
            {
                std::vector<size_t> direct_qubits;
                direct_qubits.reserve(3);
                bool direct = true;
                do
                {
                    size_t index = 0;
                    if (try_parse_direct_quantum_index(index))
                    {
                        direct_qubits.push_back(index);
                    }
                    else
                    {
                        direct = false;
                        break;
                    }
                } while (match(TokenType::COMMA));

                if (direct && match(TokenType::SEMICOLON))
                {
                    if (parameters.empty())
                    {
                        circuit.add_operation(Operation(op_type, std::move(direct_qubits)));
                    }
                    else
                    {
                        circuit.add_operation(Operation(op_type, std::move(direct_qubits), parameters));
                    }
                    return;
                }

                set_state(operands_start);
            }

            std::vector<Operand> fallback_operands;
            fallback_operands.push_back(parse_operand(/*allow_index=*/true));
            while (match(TokenType::COMMA))
            {
                fallback_operands.push_back(parse_operand(/*allow_index=*/true));
            }
            consume(TokenType::SEMICOLON, "Expected ';' after gate application.");

            apply_builtin_operands(op_type, parameters, fallback_operands);
        }

        bool try_parse_direct_quantum_index(size_t &index)
        {
            if (!check(TokenType::IDENTIFIER))
                return false;

            std::string_view name_view = advance().lexeme;
            if (!match(TokenType::LBRACKET))
                return false;

            size_t raw_index = parse_register_index();
            if (!match(TokenType::RBRACKET))
                return false;

            auto start = qubit_register_starts.find(name_view);
            if (start == qubit_register_starts.end())
                return false;

            index = start->second + raw_index;
            return true;
        }

        GateApplication parse_gate_application(bool allow_indexed_qubits, bool allow_pauli_statement)
        {
            GateApplication app;
            const Token &name_token = consume(TokenType::IDENTIFIER, "Expected gate name.");
            app.name = lowercase(name_token.lexeme);

            if (allow_pauli_statement &&
                (app.name == "t_pauli" || app.name == "s_pauli" ||
                 app.name == "z_pauli" || app.name == "m_pauli"))
            {
                parse_pauli_statement(app.name);
                app.name.clear();
                return app;
            }

            if (match(TokenType::LPAREN))
            {
                if (!check(TokenType::RPAREN))
                {
                    app.parameters.push_back(expression());
                    while (match(TokenType::COMMA))
                    {
                        app.parameters.push_back(expression());
                    }
                }
                consume(TokenType::RPAREN, "Expected ')' after gate parameters.");
            }

            app.qubits.push_back(parse_operand(allow_indexed_qubits));
            while (match(TokenType::COMMA))
            {
                app.qubits.push_back(parse_operand(allow_indexed_qubits));
            }

            consume(TokenType::SEMICOLON, "Expected ';' after gate application.");
            return app;
        }

        void parse_pauli_statement(const std::string &name)
        {
            std::string pauli_string;
            if (match(TokenType::PLUS))
                pauli_string.push_back('+');
            else if (match(TokenType::MINUS))
                pauli_string.push_back('-');
            else
                pauli_string.push_back('+');

            pauli_string += to_string(consume(TokenType::IDENTIFIER, "Expected Pauli string after " + name + ": e.g., +XYZI").lexeme);
            consume(TokenType::SEMICOLON, "Expected ';' after Pauli gate.");

            Operation::Type type = Operation::Type::T_PAULI;
            if (name == "m_pauli")
                type = Operation::Type::M_PAULI;
            else if (name == "s_pauli")
                type = Operation::Type::S_PAULI;
            else if (name == "z_pauli")
                type = Operation::Type::Z_PAULI;

            PauliOp pauli_op(circuit.get_num_qubits());
            pauli_op.from_string(pauli_string);
            circuit.add_operation(Operation(type, {}, {}, {}, pauli_op));
        }

        Operand parse_operand(bool allow_index)
        {
            Operand operand;
            operand.name = to_string(consume(TokenType::IDENTIFIER, "Expected qubit argument.").lexeme);
            if (allow_index && match(TokenType::LBRACKET))
            {
                size_t raw_index = parse_register_index();
                operand.has_index = true;
                operand.index = raw_index;
                consume(TokenType::RBRACKET, "Expected ']' after index.");
            }
            return operand;
        }

        size_t parse_register_index()
        {
            if (check(TokenType::INTEGER))
            {
                return parse_size(advance(), "integer literal");
            }

            double raw_index = expression();
            if (raw_index < 0)
            {
                throw ParseError("Negative register index is not supported.", previous().line, previous().column);
            }
            return static_cast<size_t>(raw_index);
        }

        std::vector<size_t> resolve_quantum_operand(const Operand &operand) const
        {
            for (auto it = qubit_binding_stack.rbegin(); it != qubit_binding_stack.rend(); ++it)
            {
                auto binding = it->find(operand.name);
                if (binding != it->end())
                    return binding->second;
            }

            auto start = qubit_register_starts.find(std::string_view(operand.name));
            if (start == qubit_register_starts.end())
            {
                throw std::runtime_error("Unknown quantum register: " + operand.name);
            }

            if (operand.has_index)
            {
                return {start->second + operand.index};
            }

            auto size = qubit_register_sizes.find(std::string_view(operand.name));
            std::vector<size_t> indices;
            indices.reserve(size->second);
            for (size_t i = 0; i < size->second; ++i)
            {
                indices.push_back(start->second + i);
            }
            return indices;
        }

        std::vector<size_t> resolve_classical_operand(const Operand &operand) const
        {
            auto start = bit_register_starts.find(std::string_view(operand.name));
            if (start == bit_register_starts.end())
            {
                throw std::runtime_error("Unknown classical register: " + operand.name);
            }

            if (operand.has_index)
            {
                return {start->second + operand.index};
            }

            auto size = bit_register_sizes.find(std::string_view(operand.name));
            std::vector<size_t> indices;
            indices.reserve(size->second);
            for (size_t i = 0; i < size->second; ++i)
            {
                indices.push_back(start->second + i);
            }
            return indices;
        }

        void apply_gate_application(const GateApplication &app)
        {
            if (app.name.empty())
                return;

            Operation::Type op_type;
            if (Operation::try_lowercase_name_to_type(app.name, op_type))
            {
                apply_builtin_gate(app, op_type);
                return;
            }

            auto gate = gate_definitions.find(app.name);
            if (gate == gate_definitions.end())
            {
                throw std::runtime_error("Unknown gate: " + app.name);
            }

            std::unordered_map<std::string, std::vector<size_t>> bindings;
            bindings.reserve(gate->second.qubits.size());
            for (size_t i = 0; i < gate->second.qubits.size() && i < app.qubits.size(); ++i)
            {
                bindings[gate->second.qubits[i]] = resolve_quantum_operand(app.qubits[i]);
            }

            qubit_binding_stack.push_back(std::move(bindings));
            for (const auto &body_app : gate->second.body)
            {
                apply_gate_application(body_app);
            }
            qubit_binding_stack.pop_back();
        }

        void apply_builtin_gate(const GateApplication &app, Operation::Type op_type)
        {
            if (qubit_binding_stack.empty())
            {
                std::vector<size_t> direct_qubits;
                direct_qubits.reserve(app.qubits.size());
                bool direct = !app.qubits.empty();
                for (const auto &operand : app.qubits)
                {
                    if (!operand.has_index)
                    {
                        direct = false;
                        break;
                    }
                    auto start = qubit_register_starts.find(std::string_view(operand.name));
                    if (start == qubit_register_starts.end())
                    {
                        direct = false;
                        break;
                    }
                    direct_qubits.push_back(start->second + operand.index);
                }

                if (direct)
                {
                    if (app.parameters.empty())
                    {
                        circuit.add_operation(Operation(op_type, std::move(direct_qubits)));
                    }
                    else
                    {
                        circuit.add_operation(Operation(op_type, std::move(direct_qubits), app.parameters));
                    }
                    return;
                }
            }

            apply_builtin_operands(op_type, app.parameters, app.qubits);
        }

        void apply_builtin_operands(
            Operation::Type op_type,
            const std::vector<double> &parameters,
            const std::vector<Operand> &qubit_operands)
        {
            std::vector<std::vector<size_t>> operands;
            operands.reserve(qubit_operands.size());
            size_t max_size = 0;
            for (const auto &operand : qubit_operands)
            {
                operands.push_back(resolve_quantum_operand(operand));
                max_size = std::max(max_size, operands.back().size());
            }

            for (size_t i = 0; i < max_size; ++i)
            {
                std::vector<size_t> qubits;
                qubits.reserve(operands.size());
                for (const auto &operand : operands)
                {
                    qubits.push_back(operand.size() > 1 ? operand[i] : operand[0]);
                }

                if (parameters.empty())
                {
                    circuit.add_operation(Operation(op_type, std::move(qubits)));
                }
                else
                {
                    circuit.add_operation(Operation(op_type, std::move(qubits), parameters));
                }
            }
        }

        double expression() { return additive_expr(); }

        double additive_expr()
        {
            double value = multiplicative_expr();
            while (match(TokenType::PLUS) || match(TokenType::MINUS))
            {
                TokenType op = previous().type;
                double right = multiplicative_expr();
                value = (op == TokenType::PLUS) ? value + right : value - right;
            }
            return value;
        }

        double multiplicative_expr()
        {
            double value = unary_expr();
            while (match(TokenType::TIMES) || match(TokenType::DIVIDE))
            {
                TokenType op = previous().type;
                double right = unary_expr();
                value = (op == TokenType::TIMES) ? value * right : value / right;
            }
            return value;
        }

        double unary_expr()
        {
            if (match(TokenType::MINUS))
                return -unary_expr();
            return power_expr();
        }

        double power_expr()
        {
            double value = primary_expr();
            if (match(TokenType::POWER))
            {
                value = std::pow(value, unary_expr());
            }
            return value;
        }

        double primary_expr()
        {
            if (match(TokenType::INTEGER))
                return static_cast<double>(parse_size(previous(), "integer literal"));

            if (match(TokenType::REAL))
                return std::stod(to_string(previous().lexeme));

            if (match(TokenType::PI))
                return M_PI;

            if (match(TokenType::LPAREN))
            {
                double value = expression();
                consume(TokenType::RPAREN, "Expected ')' after expression.");
                return value;
            }

            if (match(TokenType::IDENTIFIER))
                return 0.0;

            throw ParseError("Expected expression.", peek().line, peek().column);
        }
    };
}
