#pragma once

#include "token.hpp"
#include <vector>
#include <memory>
#include <string>
#include <map>
#include <unordered_map>
#include <cmath>
#include <iostream>
#include <algorithm>

namespace NWQEC
{

    // Forward declarations
    class Expr;
    class Stmt;

    // Value types that can be used in QASM expressions
    // Using simple union-like approach instead of std::variant for C++11 compatibility
    class Value
    {
    private:
        enum class Type
        {
            DOUBLE,
            INT,
            STRING,
            BOOL
        };

        Type type;
        union
        {
            double double_value;
            int int_value;
            bool bool_value;
        };
        std::string string_value; // Outside union since it has non-trivial constructor

    public:
        Value(double value) : type(Type::DOUBLE), double_value(value) {}
        Value(int value) : type(Type::INT), int_value(value) {}
        Value(bool value) : type(Type::BOOL), bool_value(value) {}
        Value(const std::string &value) : type(Type::STRING), string_value(value) {}
        Value(const char *value) : type(Type::STRING), string_value(value) {}

        bool is_double() const { return type == Type::DOUBLE; }
        bool is_int() const { return type == Type::INT; }
        bool is_string() const { return type == Type::STRING; }
        bool is_bool() const { return type == Type::BOOL; }

        double as_double() const
        {
            switch (type)
            {
            case Type::DOUBLE:
                return double_value;
            case Type::INT:
                return static_cast<double>(int_value);
            case Type::BOOL:
                return bool_value ? 1.0 : 0.0;
            case Type::STRING:
                return 0.0; // Can't convert string to double
            }
            return 0.0;
        }

        int as_int() const
        {
            switch (type)
            {
            case Type::DOUBLE:
                return static_cast<int>(double_value);
            case Type::INT:
                return int_value;
            case Type::BOOL:
                return bool_value ? 1 : 0;
            case Type::STRING:
                return 0; // Can't convert string to int
            }
            return 0;
        }

        const std::string &as_string() const
        {
            return string_value;
        }

        bool as_bool() const
        {
            switch (type)
            {
            case Type::DOUBLE:
                return double_value != 0.0;
            case Type::INT:
                return int_value != 0;
            case Type::BOOL:
                return bool_value;
            case Type::STRING:
                return !string_value.empty();
            }
            return false;
        }
    };

    /**
     * Base class for all AST expressions
     */
    class Expr
    {
    public:
        virtual ~Expr() = default;
        virtual Value evaluate() const = 0;

        // Clone method for creating deep copies
        virtual Expr *clone() const = 0;
    };

    /**
     * Number literal expression (either integer or real)
     */
    class NumberExpr : public Expr
    {
    private:
        double value;
        bool is_integer;

    public:
        NumberExpr(double value, bool is_integer) : value(value), is_integer(is_integer) {}

        Value evaluate() const override
        {
            if (is_integer)
            {
                return static_cast<int>(value);
            }
            return value;
        }

        Expr *clone() const override
        {
            return new NumberExpr(value, is_integer);
        }
    };

    /**
     * Pi constant expression
     */
    class PiExpr : public Expr
    {
    public:
        PiExpr() = default;

        Value evaluate() const override
        {
            return M_PI;
        }

        Expr *clone() const override
        {
            return new PiExpr();
        }
    };

    /**
     * Variable reference expression
     */
    class VariableExpr : public Expr
    {
    private:
        std::string name;

    public:
        VariableExpr(std::string name) : name(std::move(name)) {}

        Value evaluate() const override
        {
            // In a real implementation, this would look up the variable value
            // For now, just return the name
            return name;
        }

        const std::string &get_name() const
        {
            return name;
        }

        Expr *clone() const override
        {
            return new VariableExpr(name);
        }
    };

    /**
     * Binary operation expression (e.g., a + b, a * b)
     */
    class BinaryExpr : public Expr
    {
    public:
        enum class Op
        {
            PLUS,
            MINUS,
            MULTIPLY,
            DIVIDE,
            POWER
        };

    private:
        std::unique_ptr<Expr> left;
        Op op;
        std::unique_ptr<Expr> right;

    public:
        BinaryExpr(std::unique_ptr<Expr> left, Op op, std::unique_ptr<Expr> right)
            : left(std::move(left)), op(op), right(std::move(right)) {}

        Value evaluate() const override
        {
            Value left_val = left->evaluate();
            Value right_val = right->evaluate();

            // For simplicity, convert both values to numeric
            double left_num = left_val.as_double();
            double right_num = right_val.as_double();

            switch (op)
            {
            case Op::PLUS:
                return left_num + right_num;
            case Op::MINUS:
                return left_num - right_num;
            case Op::MULTIPLY:
                return left_num * right_num;
            case Op::DIVIDE:
                return left_num / right_num;
            case Op::POWER:
                return std::pow(left_num, right_num);
            default:
                throw std::runtime_error("Unknown binary operator");
            }
        }

        Expr *clone() const override
        {
            return new BinaryExpr(
                std::unique_ptr<Expr>(left->clone()),
                op,
                std::unique_ptr<Expr>(right->clone()));
        }
    };

    /**
     * Parenthesized expression
     */
    class ParenExpr : public Expr
    {
    private:
        std::unique_ptr<Expr> expression;

    public:
        ParenExpr(std::unique_ptr<Expr> expression) : expression(std::move(expression)) {}

        Value evaluate() const override
        {
            return expression->evaluate();
        }

        Expr *clone() const override
        {
            return new ParenExpr(std::unique_ptr<Expr>(expression->clone()));
        }
    };

    /**
     * Register indexing expression (e.g., q[0], x[i+1])
     */
    class IndexExpr : public Expr
    {
    private:
        std::string name;
        std::unique_ptr<Expr> index;

    public:
        IndexExpr(std::string name, std::unique_ptr<Expr> index)
            : name(std::move(name)), index(std::move(index)) {}

        Value evaluate() const override
        {
            // In a real implementation, this would evaluate the index and look up the register value
            // For now, just return the name
            return name + "[index]";
        }

        const std::string &get_name() const
        {
            return name;
        }

        const Expr *get_index() const
        {
            return index.get();
        }

        Expr *clone() const override
        {
            return new IndexExpr(name, std::unique_ptr<Expr>(index->clone()));
        }
    };

    /**
     * Base class for all AST statements
     */
    class Stmt
    {
    public:
        virtual ~Stmt() = default;
        virtual void execute() const = 0;

        // Clone method for creating deep copies
        virtual Stmt *clone() const = 0;
    };

    /**
     * QASM version declaration (e.g., OPENQASM 2.0;)
     */
    class VersionDecl : public Stmt
    {
    private:
        std::string version;

    public:
        VersionDecl(std::string version) : version(std::move(version)) {}

        void execute() const override
        {
            // Set QASM version
        }

        const std::string &get_version() const
        {
            return version;
        }

        Stmt *clone() const override
        {
            return new VersionDecl(version);
        }
    };

    /**
     * Include statement (e.g., include "qelib1.inc";)
     */
    class IncludeStmt : public Stmt
    {
    private:
        std::string filename;

    public:
        IncludeStmt(std::string filename) : filename(std::move(filename)) {}

        void execute() const override
        {
            // Include file
        }

        const std::string &get_filename() const
        {
            return filename;
        }

        Stmt *clone() const override
        {
            return new IncludeStmt(filename);
        }
    };

    /**
     * Quantum register declaration (e.g., qreg q[5];)
     */
    class QRegDecl : public Stmt
    {
    private:
        std::string name;
        size_t size;

    public:
        QRegDecl(std::string name, size_t size) : name(std::move(name)), size(size) {}

        void execute() const override
        {
            // Create quantum register
        }

        const std::string &get_name() const
        {
            return name;
        }

        size_t get_size() const
        {
            return size;
        }

        Stmt *clone() const override
        {
            return new QRegDecl(name, size);
        }
    };

    /**
     * Classical register declaration (e.g., creg c[5];)
     */
    class CRegDecl : public Stmt
    {
    private:
        std::string name;
        size_t size;

    public:
        CRegDecl(std::string name, size_t size) : name(std::move(name)), size(size) {}

        void execute() const override
        {
            // Create classical register
        }

        const std::string &get_name() const
        {
            return name;
        }

        size_t get_size() const
        {
            return size;
        }

        Stmt *clone() const override
        {
            return new CRegDecl(name, size);
        }
    };

    /**
     * Quantum gate application (e.g., h q[0];, cx q[0],q[1];)
     */
    class GateStmt : public Stmt
    {
    private:
        std::string name;
        std::vector<std::unique_ptr<Expr>> parameters;
        std::vector<std::unique_ptr<Expr>> qubits;

    public:
        GateStmt(std::string name,
                 std::vector<std::unique_ptr<Expr>> parameters,
                 std::vector<std::unique_ptr<Expr>> qubits)
            : name(std::move(name)),
              parameters(std::move(parameters)),
              qubits(std::move(qubits)) {}

        void execute() const override
        {
            // Apply gate
        }

        const std::string &get_name() const
        {
            return name;
        }

        const std::vector<std::unique_ptr<Expr>> &get_parameters() const
        {
            return parameters;
        }

        const std::vector<std::unique_ptr<Expr>> &get_qubits() const
        {
            return qubits;
        }

        Stmt *clone() const override
        {
            std::vector<std::unique_ptr<Expr>> paramsCopy;
            for (const auto &param : parameters)
            {
                paramsCopy.push_back(std::unique_ptr<Expr>(dynamic_cast<Expr *>(param->clone())));
            }

            std::vector<std::unique_ptr<Expr>> qubitsCopy;
            for (const auto &qubit : qubits)
            {
                qubitsCopy.push_back(std::unique_ptr<Expr>(dynamic_cast<Expr *>(qubit->clone())));
            }

            return new GateStmt(name, std::move(paramsCopy), std::move(qubitsCopy));
        }
    };

    /**
     * Measurement statement (e.g., measure q[0] -> c[0];)
     */
    class MeasureStmt : public Stmt
    {
    private:
        std::unique_ptr<Expr> qubit;
        std::unique_ptr<Expr> bit;

    public:
        MeasureStmt(std::unique_ptr<Expr> qubit, std::unique_ptr<Expr> bit)
            : qubit(std::move(qubit)), bit(std::move(bit)) {}

        void execute() const override
        {
            // Perform measurement
        }

        const Expr *get_qubit() const
        {
            return qubit.get();
        }

        const Expr *get_bit() const
        {
            return bit.get();
        }

        Stmt *clone() const override
        {
            return new MeasureStmt(
                std::unique_ptr<Expr>(dynamic_cast<Expr *>(qubit->clone())),
                std::unique_ptr<Expr>(dynamic_cast<Expr *>(bit->clone())));
        }
    };

    /**
     * Reset statement (e.g., reset q[0];)
     */
    class ResetStmt : public Stmt
    {
    private:
        std::unique_ptr<Expr> qubit;

    public:
        ResetStmt(std::unique_ptr<Expr> qubit) : qubit(std::move(qubit)) {}

        void execute() const override
        {
            // Reset qubit
        }

        const Expr *get_qubit() const
        {
            return qubit.get();
        }

        Stmt *clone() const override
        {
            return new ResetStmt(std::unique_ptr<Expr>(dynamic_cast<Expr *>(qubit->clone())));
        }
    };

    /**
     * Custom gate declaration (e.g., gate h q { u(pi/2,0,pi) q; })
     */
    class GateDeclStmt : public Stmt
    {
    private:
        std::string name;
        std::vector<std::string> params;
        std::vector<std::string> qubits;
        std::vector<std::unique_ptr<Stmt>> body;

    public:
        GateDeclStmt(std::string name,
                     std::vector<std::string> params,
                     std::vector<std::string> qubits,
                     std::vector<std::unique_ptr<Stmt>> body)
            : name(std::move(name)),
              params(std::move(params)),
              qubits(std::move(qubits)),
              body(std::move(body)) {}

        void execute() const override
        {
            // Define custom gate
        }

        const std::string &get_name() const
        {
            return name;
        }

        const std::vector<std::string> &get_params() const
        {
            return params;
        }

        const std::vector<std::string> &get_qubits() const
        {
            return qubits;
        }

        const std::vector<std::unique_ptr<Stmt>> &get_body() const
        {
            return body;
        }

        Stmt *clone() const override
        {
            std::vector<std::string> paramsCopy = params;
            std::vector<std::string> qubitsCopy = qubits;

            std::vector<std::unique_ptr<Stmt>> bodyCopy;
            for (const auto &stmt : body)
            {
                bodyCopy.push_back(std::unique_ptr<Stmt>(stmt->clone()));
            }

            return new GateDeclStmt(name, std::move(paramsCopy), std::move(qubitsCopy), std::move(bodyCopy));
        }
    };

    /**
     * Conditional statement (e.g., if(c==1) x q[0];)
     */
    class IfStmt : public Stmt
    {
    private:
        std::string creg;
        int value;
        std::unique_ptr<Stmt> then_branch;

    public:
        IfStmt(std::string creg, int value, std::unique_ptr<Stmt> then_branch)
            : creg(std::move(creg)), value(value), then_branch(std::move(then_branch)) {}

        void execute() const override
        {
            // Execute conditional statement
        }

        const std::string &get_creg() const
        {
            return creg;
        }

        int get_value() const
        {
            return value;
        }

        const Stmt *get_then_branch() const
        {
            return then_branch.get();
        }

        Stmt *clone() const override
        {
            return new IfStmt(creg, value, std::unique_ptr<Stmt>(then_branch->clone()));
        }
    };

    /**
     * Block statement containing multiple statements (e.g., { x q[0]; h q[1]; })
     */
    class BlockStmt : public Stmt
    {
    private:
        std::vector<std::unique_ptr<Stmt>> statements;

    public:
        BlockStmt(std::vector<std::unique_ptr<Stmt>> statements)
            : statements(std::move(statements)) {}

        void execute() const override
        {
            // Execute all statements in the block
            for (const auto &stmt : statements)
            {
                stmt->execute();
            }
        }

        const std::vector<std::unique_ptr<Stmt>> &get_statements() const
        {
            return statements;
        }

        Stmt *clone() const override
        {
            std::vector<std::unique_ptr<Stmt>> statementsCopy;
            for (const auto &stmt : statements)
            {
                statementsCopy.push_back(std::unique_ptr<Stmt>(stmt->clone()));
            }

            return new BlockStmt(std::move(statementsCopy));
        }
    };

    /**
     * Barrier statement (e.g., barrier q[0],q[1];)
     */
    class BarrierStmt : public Stmt
    {
    private:
        std::vector<std::unique_ptr<Expr>> qubits;

    public:
        BarrierStmt(std::vector<std::unique_ptr<Expr>> qubits) : qubits(std::move(qubits)) {}

        void execute() const override
        {
            // Execute barrier
        }

        const std::vector<std::unique_ptr<Expr>> &get_qubits() const
        {
            return qubits;
        }

        Stmt *clone() const override
        {
            std::vector<std::unique_ptr<Expr>> qubitsCopy;
            for (const auto &qubit : qubits)
            {
                qubitsCopy.push_back(std::unique_ptr<Expr>(qubit->clone()));
            }
            return new BarrierStmt(std::move(qubitsCopy));
        }
    };

    /**
     * Pauli gate statement (e.g., t_pauli +IXYZ;, s_pauli +XYZI;, z_pauli -ZYXI;, m_pauli +ZZII;)
     */
    class PauliStmt : public Stmt
    {
    private:
        std::string gate_name;  // "t_pauli", "s_pauli", "z_pauli", or "m_pauli"
        std::string pauli_string;  // The Pauli string like "+IXYZ" or "-ZYXI"

    public:
        PauliStmt(const std::string &gate_name, const std::string &pauli_string)
            : gate_name(gate_name), pauli_string(pauli_string) {}

        void execute() const override
        {
            // Execute Pauli gate operation
        }

        const std::string &get_gate_name() const
        {
            return gate_name;
        }

        const std::string &get_pauli_string() const
        {
            return pauli_string;
        }

        Stmt *clone() const override
        {
            return new PauliStmt(gate_name, pauli_string);
        }
    };

    /**
     * QASM ASTprogram, which consists of a sequence of statements
     */
    class ASTProgram
    {
    private:
        std::vector<std::unique_ptr<Stmt>> statements;

    public:
        ASTProgram() = default;

        void add_statement(std::unique_ptr<Stmt> stmt)
        {
            statements.push_back(std::move(stmt));
        }

        void execute() const
        {
            // Execute all statements in the program
            for (const auto &stmt : statements)
            {
                stmt->execute();
            }
        }

        const std::vector<std::unique_ptr<Stmt>> &get_statements() const
        {
            return statements;
        }
    };

} // namespace NWQEC