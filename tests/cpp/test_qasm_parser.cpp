#include "nwqec/parser/qasm_parser.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace
{
    using NWQEC::Circuit;
    using NWQEC::Operation;
    using NWQEC::QASMParser;

    void require(bool condition, const std::string &message)
    {
        if (!condition)
        {
            throw std::runtime_error(message);
        }
    }

    std::unique_ptr<Circuit> parse_string(const std::string &qasm)
    {
        QASMParser parser;
        require(parser.parse_string(qasm), "parse failed: " + parser.get_error_message());
        auto circuit = parser.get_circuit();
        require(static_cast<bool>(circuit), "parser returned null circuit");
        return circuit;
    }

    void expect_op(
        const Operation &op,
        Operation::Type type,
        const std::vector<size_t> &qubits,
        const std::vector<size_t> &bits = {},
        const std::vector<double> &params = {})
    {
        require(op.get_type() == type, "unexpected operation type");
        require(op.get_qubits() == qubits, "unexpected operation qubits");
        require(op.get_bits() == bits, "unexpected operation bits");
        require(op.get_parameters().size() == params.size(), "unexpected parameter count");
        for (size_t i = 0; i < params.size(); ++i)
        {
            require(std::fabs(op.get_parameters()[i] - params[i]) < 1e-12, "unexpected parameter value");
        }
    }

    void test_core_qasm_features()
    {
        const std::string qasm = R"qasm(
OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg c[3];
h q[0];
cx q[0],q[1];
rz(pi/2) q[2];
x q;
barrier q[0],q[2];
reset q[1];
measure q -> c;
)qasm";

        auto circuit = parse_string(qasm);
        require(circuit->get_num_qubits() == 3, "wrong qubit count");
        require(circuit->get_num_bits() == 3, "wrong bit count");

        const auto &ops = circuit->get_operations();
        require(ops.size() == 11, "wrong operation count");
        expect_op(ops[0], Operation::Type::H, {0});
        expect_op(ops[1], Operation::Type::CX, {0, 1});
        expect_op(ops[2], Operation::Type::RZ, {2}, {}, {M_PI / 2.0});
        expect_op(ops[3], Operation::Type::X, {0});
        expect_op(ops[4], Operation::Type::X, {1});
        expect_op(ops[5], Operation::Type::X, {2});
        expect_op(ops[6], Operation::Type::BARRIER, {0, 2});
        expect_op(ops[7], Operation::Type::RESET, {1});
        expect_op(ops[8], Operation::Type::MEASURE, {0}, {0});
        expect_op(ops[9], Operation::Type::MEASURE, {1}, {1});
        expect_op(ops[10], Operation::Type::MEASURE, {2}, {2});
    }

    void test_user_defined_gate_and_pauli()
    {
        const std::string qasm = R"qasm(
OPENQASM 2.0;
qreg q[3];
gate pair a,b {
  h a;
  cx a,b;
}
pair q[0],q[2];
t_pauli +XYZ;
m_pauli -ZZI;
)qasm";

        auto circuit = parse_string(qasm);
        const auto &ops = circuit->get_operations();
        require(ops.size() == 4, "wrong user gate operation count");
        expect_op(ops[0], Operation::Type::H, {0});
        expect_op(ops[1], Operation::Type::CX, {0, 2});
        require(ops[2].get_type() == Operation::Type::T_PAULI, "expected t_pauli");
        require(ops[2].get_pauli_string() == "+XYZ", "wrong t_pauli string");
        require(ops[3].get_type() == Operation::Type::M_PAULI, "expected m_pauli");
        require(ops[3].get_pauli_string() == "-ZZI", "wrong m_pauli string");
    }

    void test_parse_file_and_errors()
    {
        const std::string path = "parser_regression_tmp.qasm";
        {
            std::ofstream out(path);
            out << "OPENQASM 2.0;\nqreg q[1];\nh q[0];\n";
        }

        QASMParser file_parser;
        require(file_parser.parse_file(path), "parse_file failed: " + file_parser.get_error_message());
        auto circuit = file_parser.get_circuit();
        require(circuit && circuit->get_operations().size() == 1, "parse_file operation mismatch");
        std::remove(path.c_str());

        QASMParser bad_parser;
        require(!bad_parser.parse_string("OPENQASM 2.0;\ninvalid_statement;"),
                "invalid QASM unexpectedly parsed");
        require(!bad_parser.get_error_message().empty(), "missing parse error message");
    }
}

int main()
{
    try
    {
        test_core_qasm_features();
        test_user_defined_gate_and_pauli();
        test_parse_file_and_errors();
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}
