#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <set>
#include "pauli_op.hpp"

namespace NWQEC
{

    /**
     * Represents a quantum operation (gate, measurement, etc.) in a flattened circuit
     */
    class Operation
    {
    public:
        enum class Type
        {
            // Single-qubit gates
            X,
            Y,
            Z,
            H,
            S,
            SDG,
            T,
            TDG,
            SX,
            SXDG,
            P4,
            P8,
            P16,
            ID,
            // Parameterized single-qubit gates
            RX,
            RY,
            RZ,
            P,
            U,
            U1,
            U2,
            U3,
            // Two-qubit gates
            CX,
            CY,
            CZ,
            CH,
            CS,
            CSDG,
            CT,
            CTDG,
            CSX,
            SWAP,
            ECR,
            // Parameterized two-qubit gates
            CRX,
            CRY,
            CRZ,
            CP,
            CU,
            CU1,
            CU3,
            RXX,
            RYY,
            RZZ,
            // Three-qubit gates
            CCX,
            CSWAP,
            RCCX,
            // Measurement
            MEASURE,
            // Reset
            RESET,
            // Barrier
            BARRIER,
            T_PAULI,
            M_PAULI,
            S_PAULI,
            Z_PAULI,
            SWAP_BASIS,
        };

    private:
        Type type;
        std::vector<size_t> qubits;     // Global qubit indices
        std::vector<double> parameters; // Parameters for parameterized gates
        std::vector<size_t> bits;       // Classical bit indices (for measurement)
        PauliOp pauli_op;               // For T_PAULI, M_PAULI, S_PAULI, and Z_PAULI operations
        bool dagger;                    // Whether this is the dagger (conjugate transpose) of the operation
        bool x_rotation;                // Whether this operation includes x rotation

        std::vector<size_t> active_qubits(const PauliOp &pauli_op) const
        {
            std::vector<size_t> involved_qubits;
            
            // Get X and Z indices from PauliOp and merge them
            auto x_indices = pauli_op.get_x_indices();
            auto z_indices = pauli_op.get_z_indices();
            
            // Combine and sort unique indices
            std::set<size_t> unique_indices;
            unique_indices.insert(x_indices.begin(), x_indices.end());
            unique_indices.insert(z_indices.begin(), z_indices.end());
            
            involved_qubits.assign(unique_indices.begin(), unique_indices.end());
            return involved_qubits;
        }

    public:
        Operation(Type type,
                  std::vector<size_t> qubits,
                  std::vector<double> parameters = {},
                  std::vector<size_t> bits = {},
                  PauliOp pauli_op = PauliOp(),
                  bool dagger = false,
                  bool x_rotation = false)
            : type(type),
              qubits(qubits.empty() ? active_qubits(pauli_op) : std::move(qubits)),
              parameters(std::move(parameters)),
              bits(std::move(bits)),
              pauli_op(std::move(pauli_op)),
              dagger(dagger),
              x_rotation(x_rotation) {}

        Type get_type() const { return type; }
        const std::vector<size_t> &get_qubits() const { return qubits; }
        const std::vector<size_t> &get_bits() const { return bits; }
        const std::vector<double> &get_parameters() const { return parameters; }
        const PauliOp &get_pauli_op() const { return pauli_op; }
        std::string get_pauli_string() const { return pauli_op.to_string(); }
        bool get_dagger() const { return dagger; }
        bool get_x_rotation() const { return x_rotation; }

        // Get string representation of the operation type
        std::string get_type_name() const
        {
            // Handle Pn types with flags
            if (type == Type::P4 || type == Type::P8 || type == Type::P16)
            {
                int denominator = (type == Type::P4) ? 4 : (type == Type::P8) ? 8
                                                                              : 16;
                std::string sign = dagger ? "-" : "";
                std::string axis = x_rotation ? "rx" : "rz";

                // Special cases for P4 without x_rotation
                if (type == Type::P4 && !x_rotation)
                {
                    return dagger ? "tdg" : "t";
                }

                return axis + "(" + sign + "pi/" + std::to_string(denominator) + ")";
            }
            return get_type_name(type);
        }

        static std::string get_type_name(Operation::Type gate_type)
        {
            switch (gate_type)
            {
                // Single-qubit gates
            case Type::X:
                return "x";
            case Type::Y:
                return "y";
            case Type::Z:
                return "z";
            case Type::H:
                return "h";
            case Type::S:
                return "s";
            case Type::SDG:
                return "sdg";
            case Type::T:
                return "t";
            case Type::TDG:
                return "tdg";
            case Type::SX:
                return "sx";
            case Type::SXDG:
                return "sxdg";
            case Type::P4:
                return "p4";
            case Type::P8:
                return "p8";
            case Type::P16:
                return "p16";
            case Type::ID:
                return "id";
            // Parameterized single-qubit gates
            case Type::RX:
                return "rx";
            case Type::RY:
                return "ry";
            case Type::RZ:
                return "rz";
            case Type::P:
                return "p";
            case Type::U:
                return "u";
            case Type::U1:
                return "u1";
            case Type::U2:
                return "u2";
            case Type::U3:
                return "u3";
            // Two-qubit gates
            case Type::CX:
                return "cx";
            case Type::CY:
                return "cy";
            case Type::CZ:
                return "cz";
            case Type::CH:
                return "ch";
            case Type::CS:
                return "cs";
            case Type::CSDG:
                return "csdg";
            case Type::CT:
                return "ct";
            case Type::CTDG:
                return "ctdg";
            case Type::CSX:
                return "csx";
            case Type::SWAP:
                return "swap";
            case Type::ECR:
                return "ecr";
            // Parameterized two-qubit gates
            case Type::CRX:
                return "crx";
            case Type::CRY:
                return "cry";
            case Type::CRZ:
                return "crz";
            case Type::CP:
                return "cp";
            case Type::CU:
                return "cu";
            case Type::CU1:
                return "cu1";
            case Type::CU3:
                return "cu3";
            case Type::RXX:
                return "rxx";
            case Type::RYY:
                return "ryy";
            case Type::RZZ:
                return "rzz";
            // Three-qubit gates
            case Type::CCX:
                return "ccx";
            case Type::CSWAP:
                return "cswap";
            case Type::RCCX:
                return "rccx";
            // Measurement
            case Type::MEASURE:
                return "measure";
            case Type::M_PAULI:
                return "m_pauli";
            // Reset
            case Type::RESET:
                return "reset";
            // Barrier
            case Type::BARRIER:
                return "barrier";
            case Type::T_PAULI:
                return "t_pauli";
            case Type::S_PAULI:
                return "s_pauli";
            case Type::Z_PAULI:
                return "z_pauli";
            case Type::SWAP_BASIS:
                return "swap_basis";
            default:
                return "unknown";
            }
        }

        // Convert string gate name to operation type
        static Type name_to_type(const std::string &name)
        {
            std::string lowercase_name = name;
            std::transform(lowercase_name.begin(), lowercase_name.end(), lowercase_name.begin(), ::tolower);

            if (lowercase_name == "x")
                return Type::X;
            if (lowercase_name == "y")
                return Type::Y;
            if (lowercase_name == "z")
                return Type::Z;
            if (lowercase_name == "h")
                return Type::H;
            if (lowercase_name == "s")
                return Type::S;
            if (lowercase_name == "sdg")
                return Type::SDG;
            if (lowercase_name == "t")
                return Type::T;
            if (lowercase_name == "tdg")
                return Type::TDG;
            if (lowercase_name == "sx")
                return Type::SX;
            if (lowercase_name == "sxdg")
                return Type::SXDG;
            if (lowercase_name == "id")
                return Type::ID;

            // Parameterized single-qubit gates
            if (lowercase_name == "rx")
                return Type::RX;
            if (lowercase_name == "ry")
                return Type::RY;
            if (lowercase_name == "rz")
                return Type::RZ;
            if (lowercase_name == "p")
                return Type::P;
            if (lowercase_name == "u")
                return Type::U;
            if (lowercase_name == "u1")
                return Type::U1;
            if (lowercase_name == "u2")
                return Type::U2;
            if (lowercase_name == "u3")
                return Type::U3;

            // Two-qubit gates
            if (lowercase_name == "cx")
                return Type::CX;
            if (lowercase_name == "cy")
                return Type::CY;
            if (lowercase_name == "cz")
                return Type::CZ;
            if (lowercase_name == "ch")
                return Type::CH;
            if (lowercase_name == "cs")
                return Type::CS;
            if (lowercase_name == "csdg")
                return Type::CSDG;
            if (lowercase_name == "ct")
                return Type::CT;
            if (lowercase_name == "ctdg")
                return Type::CTDG;
            if (lowercase_name == "csx")
                return Type::CSX;
            if (lowercase_name == "swap")
                return Type::SWAP;
            if (lowercase_name == "ecr")
                return Type::ECR;

            // Parameterized two-qubit gates
            if (lowercase_name == "crx")
                return Type::CRX;
            if (lowercase_name == "cry")
                return Type::CRY;
            if (lowercase_name == "crz")
                return Type::CRZ;
            if (lowercase_name == "cp")
                return Type::CP;
            if (lowercase_name == "cu")
                return Type::CU;
            if (lowercase_name == "cu1")
                return Type::CU1;
            if (lowercase_name == "cu3")
                return Type::CU3;
            if (lowercase_name == "rxx")
                return Type::RXX;
            if (lowercase_name == "ryy")
                return Type::RYY;
            if (lowercase_name == "rzz")
                return Type::RZZ;

            // Three-qubit gates
            if (lowercase_name == "ccx")
                return Type::CCX;
            if (lowercase_name == "cswap")
                return Type::CSWAP;
            if (lowercase_name == "rccx")
                return Type::RCCX;
            // Measurement
            if (lowercase_name == "measure")
                return Type::MEASURE;
            // Reset
            if (lowercase_name == "reset")
                return Type::RESET;
            // Barrier
            if (lowercase_name == "barrier")
                return Type::BARRIER;

            if (lowercase_name == "t_pauli")
                return Type::T_PAULI;
            if (lowercase_name == "m_pauli")
                return Type::M_PAULI;
            if (lowercase_name == "s_pauli")
                return Type::S_PAULI;
            if (lowercase_name == "z_pauli")
                return Type::Z_PAULI;
            if (lowercase_name == "swap_basis")
                return Type::SWAP_BASIS;

            throw std::runtime_error("Unknown gate: " + name);
        }

        static bool is_builtin_gate(const std::string &name)
        {
            try
            {
                name_to_type(name);
                return true;
            }
            catch (const std::runtime_error &)
            {
                return false;
            }
        }

        std::string get_parameter_string(double param_value, int precision = 10, double eps = 1e-10) const
        {
            const double pi = M_PI;

            double multiplier = param_value / pi;

            // Attempt to approximate multiplier as a rational number numerator/denominator
            const int max_denominator = 100; // Reasonable limit to keep fractions simple

            int best_numerator = 0;
            int best_denominator = 1;
            double min_error = std::abs(multiplier); // Start with the worst case

            for (int denom = 1; denom <= max_denominator; ++denom)
            {
                int num = static_cast<int>(std::round(multiplier * denom));
                double error = std::fabs(multiplier - static_cast<double>(num) / denom);

                if (error < min_error - eps)
                {
                    min_error = error;
                    best_numerator = num;
                    best_denominator = denom;

                    if (min_error < eps)
                        break; // Good enough approximation found
                }
            }

            // If the approximation is within acceptable error, return symbolic representation
            if (min_error < eps)
            {
                // Simplify fraction
                int gcd = std::gcd(std::abs(best_numerator), best_denominator);
                best_numerator /= gcd;
                best_denominator /= gcd;

                // Handle special cases
                if (best_numerator == 0)
                    return "0";
                if (best_denominator == 1)
                {
                    if (best_numerator == 1)
                        return "pi";
                    if (best_numerator == -1)
                        return "-pi";
                    return std::to_string(best_numerator) + "pi";
                }

                std::ostringstream oss;
                if (best_numerator == 1)
                    oss << "pi/" << best_denominator;
                else if (best_numerator == -1)
                    oss << "-pi/" << best_denominator;
                else
                    oss << best_numerator << "*pi/" << best_denominator;

                return oss.str();
            }

            // Fallback to decimal representation with precision
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(precision) << param_value;
            std::string result = oss.str();

            // Trim trailing zeros and possible trailing decimal point
            result.erase(result.find_last_not_of('0') + 1);
            if (!result.empty() && result.back() == '.')
            {
                result.pop_back();
            }

            return result;
        }

        // Print the operation in QASM format
        void print(std::ostream &os) const
        {
            os << get_type_name();

            if (type == Type::T_PAULI || type == Type::M_PAULI || type == Type::S_PAULI || type == Type::Z_PAULI)
            {
                // Special case for T_PAULI, M_PAULI, S_PAULI, and Z_PAULI
                os << " " << pauli_op.to_string() << ";";
                return;
            }
            
            if (type == Type::SWAP_BASIS)
            {
                // Special case for SWAP_BASIS - single qubit operation
                if (!qubits.empty())
                {
                    os << " q[" << qubits[0] << "]";
                }
                os << ";";
                return;
            }

            // Print parameters if any
            if (!parameters.empty())
            {
                os << "(";
                for (size_t i = 0; i < parameters.size(); ++i)
                {
                    os << get_parameter_string(parameters[i]);
                    if (i < parameters.size() - 1)
                        os << ",";
                }
                os << ")";
            }

            // Print qubits
            os << " ";
            for (size_t i = 0; i < qubits.size(); ++i)
            {
                os << "q[" << qubits[i] << "]";
                if (i < qubits.size() - 1)
                    os << ",";
            }

            // Print bits for measurement
            if (type == Type::MEASURE && !bits.empty())
            {
                os << " -> ";
                for (size_t i = 0; i < bits.size(); ++i)
                {
                    os << "c[" << bits[i] << "]";
                    if (i < bits.size() - 1)
                        os << ",";
                }
            }

            os << ";";
        }
    };

} // namespace NWQEC