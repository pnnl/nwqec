// Python bindings for NWQEC using pybind11

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "nwqec/parser/qasm_parser.hpp"
#include "nwqec/core/operation.hpp"
#include "nwqec/core/pauli_op.hpp"
#include "nwqec/core/constants.hpp"

#include "nwqec/core/transpiler.hpp"

namespace py = pybind11;

namespace
{
    // Helper to render stats to a string
    std::string circuit_stats(const NWQEC::Circuit &c)
    {
        std::ostringstream oss;
        c.print_stats(oss);
        return oss.str();
    }
    std::string circuit_to_qasm(const NWQEC::Circuit &c)
    {
        std::ostringstream oss;
        c.print(oss);
        return oss.str();
    }
    void circuit_save_qasm(const NWQEC::Circuit &c, const std::string &filename)
    {
        std::ofstream ofs(filename);
        if (!ofs)
        {
            throw std::runtime_error("Failed to open file for writing: " + filename);
        }
        c.print(ofs);
    }

    // Helpers to enforce PBC vs standard gate exclusivity in Python API
    inline bool is_pauli_op(NWQEC::Operation::Type t)
    {
        using T = NWQEC::Operation::Type;
        return t == T::T_PAULI || t == T::S_PAULI || t == T::Z_PAULI || t == T::M_PAULI;
    }

    inline bool is_barrier(NWQEC::Operation::Type t)
    {
        return t == NWQEC::Operation::Type::BARRIER;
    }

    bool circuit_has_pauli_ops(const NWQEC::Circuit &c)
    {
        for (const auto &op : c.get_operations())
        {
            if (is_pauli_op(op.get_type()))
                return true;
        }
        return false;
    }

    bool circuit_has_non_pauli_ops(const NWQEC::Circuit &c)
    {
        for (const auto &op : c.get_operations())
        {
            if (!is_pauli_op(op.get_type()) && !is_barrier(op.get_type()))
                return true;
        }
        return false;
    }

    // Helper to run transforms using the Transpiler
    std::unique_ptr<NWQEC::Circuit> apply_transforms(const NWQEC::Circuit &circuit,
                                                     bool to_pbc,
                                                     bool to_clifford_reduction,
                                                     bool keep_cx,
                                                     bool t_pauli_opt,
                                                     bool remove_pauli,
                                                     bool keep_ccx,
                                                     bool silent,
                                                     double epsilon_override = -1.0)
    {
        NWQEC::Transpiler transpiler;
        NWQEC::PassConfig config;
        config.keep_ccx = keep_ccx;
        config.keep_cx = keep_cx;
        config.epsilon_override = epsilon_override;
        config.silent = silent;
        
        auto circuit_copy = std::make_unique<NWQEC::Circuit>(circuit);
        
        // Choose the appropriate pass sequence
        std::vector<NWQEC::PassType> passes;
        
        if (t_pauli_opt) {
            // Just T-optimization for PBC circuits
            passes = {NWQEC::PassType::TFUSE};
        } else if (to_pbc) {
            passes = NWQEC::PassSequences::TO_PBC_BASIC;
        } else if (to_clifford_reduction) {
            passes = NWQEC::PassSequences::CLIFFORD_REDUCTION;
        } else {
            // Default: Clifford+T (always available)
            passes = NWQEC::PassSequences::TO_CLIFFORD_T;
        }
        
        // Add cleanup passes if requested
        if (remove_pauli) {
            passes.push_back(NWQEC::PassType::REMOVE_PAULI);
        }
        
        return transpiler.execute_passes(std::move(circuit_copy), passes, config);
    }
}

PYBIND11_MODULE(_core, m)
{
    m.doc() = "NWQEC Python bindings";

#ifdef NWQEC_WITH_GRIDSYNTH_CPP
    m.attr("WITH_GRIDSYNTH_CPP") = py::bool_(NWQEC_WITH_GRIDSYNTH_CPP);
#else
    m.attr("WITH_GRIDSYNTH_CPP") = py::bool_(false);
#endif

    // Circuit class (owned by Python via unique_ptr)
    py::class_<NWQEC::Circuit, std::unique_ptr<NWQEC::Circuit>>(m, "Circuit")
        // Circuit constructor
        .def(py::init([](size_t num_qubits)
                      {
                 auto c = std::make_unique<NWQEC::Circuit>();
                 if (num_qubits > 0) c->add_qreg("q", num_qubits);
                 return c; }),
             py::arg("num_qubits"))
        .def("x", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             {
                 if (circuit_has_pauli_ops(c))
                     throw std::runtime_error("Cannot mix Pauli-based operations with standard gates in one circuit (PBC-only).");
                 c.add_operation({NWQEC::Operation::Type::X, {q}});
                 return c; }, py::arg("q"), "Apply Pauli-X to qubit q.")
        .def("y", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             {
                 if (circuit_has_pauli_ops(c))
                     throw std::runtime_error("Cannot mix Pauli-based operations with standard gates in one circuit (PBC-only).");
                 c.add_operation({NWQEC::Operation::Type::Y, {q}});
                 return c; }, py::arg("q"), "Apply Pauli-Y to qubit q.")
        .def("z", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             {
                 if (circuit_has_pauli_ops(c))
                     throw std::runtime_error("Cannot mix Pauli-based operations with standard gates in one circuit (PBC-only).");
                 c.add_operation({NWQEC::Operation::Type::Z, {q}});
                 return c; }, py::arg("q"), "Apply Pauli-Z to qubit q.")
        .def("h", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::H, {q}}); return c; }, py::arg("q"), "Apply Hadamard to qubit q.")
        .def("s", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::S, {q}}); return c; }, py::arg("q"), "Apply phase S (π/2 about Z) to qubit q.")
        .def("sdg", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::SDG, {q}}); return c; }, py::arg("q"), "Apply S† to qubit q.")
        .def("t", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::T, {q}}); return c; }, py::arg("q"), "Apply T (π/4 about Z) to qubit q.")
        .def("tdg", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::TDG, {q}}); return c; }, py::arg("q"), "Apply T† to qubit q.")
        .def("sx", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::SX, {q}}); return c; }, py::arg("q"), "Apply √X to qubit q.")
        .def("sxdg", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::SXDG, {q}}); return c; }, py::arg("q"), "Apply (√X)† to qubit q.")
        .def("cx", [](NWQEC::Circuit &c, size_t q0, size_t q1) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::CX, {q0, q1}}); return c; }, py::arg("q0"), py::arg("q1"), "Apply CX(control=q0, target=q1).")
        .def("ccx", [](NWQEC::Circuit &c, size_t q0, size_t q1, size_t q2) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::CCX, {q0, q1, q2}}); return c; }, py::arg("q0"), py::arg("q1"), py::arg("q2"), "Apply CCX(control=q0,q1; target=q2).")
        .def("cz", [](NWQEC::Circuit &c, size_t q0, size_t q1) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::CZ, {q0, q1}}); return c; }, py::arg("q0"), py::arg("q1"), "Apply CZ between q0 and q1.")
        .def("swap", [](NWQEC::Circuit &c, size_t q0, size_t q1) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::SWAP, {q0, q1}}); return c; }, py::arg("q0"), py::arg("q1"), "Swap states of q0 and q1.")
        .def("rx", [](NWQEC::Circuit &c, size_t q, double th) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::RX, {q}, {th}}); return c; }, py::arg("q"), py::arg("theta"))
        .def("rxp", [](NWQEC::Circuit &c, size_t q, double x) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::RX, {q}, {x * M_PI}}); return c; }, py::arg("q"), py::arg("x_pi"))
        .def("ry", [](NWQEC::Circuit &c, size_t q, double th) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::RY, {q}, {th}}); return c; }, py::arg("q"), py::arg("theta"))
        .def("ryp", [](NWQEC::Circuit &c, size_t q, double x) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::RY, {q}, {x * M_PI}}); return c; }, py::arg("q"), py::arg("x_pi"))
        .def("rz", [](NWQEC::Circuit &c, size_t q, double th) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::RZ, {q}, {th}}); return c; }, py::arg("q"), py::arg("theta"))
        .def("rzp", [](NWQEC::Circuit &c, size_t q, double x) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::RZ, {q}, {x * M_PI}}); return c; }, py::arg("q"), py::arg("x_pi"))
        .def("measure", [](NWQEC::Circuit &c, size_t q, size_t b) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::MEASURE, {q}, {}, {b}}); return c; }, py::arg("q"), py::arg("cbit"))
        .def("reset", [](NWQEC::Circuit &c, size_t q) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::RESET, {q}}); return c; }, py::arg("q"))
        .def("barrier", [](NWQEC::Circuit &c, const std::vector<size_t> &qs) -> NWQEC::Circuit &
             { c.add_operation({NWQEC::Operation::Type::BARRIER, qs}); return c; }, py::arg("qubits"))
        // Clean Pauli helpers: accept only a string
        .def("t_pauli", [](NWQEC::Circuit &c, const std::string &p) -> NWQEC::Circuit &
             {
                if (circuit_has_non_pauli_ops(c))
                    throw std::runtime_error("Pauli-based operations are valid only in PBC circuits; do not mix with standard gates.");
                NWQEC::PauliOp pop(c.get_num_qubits()); pop.from_string(p);
                c.add_operation(NWQEC::Operation(NWQEC::Operation::Type::T_PAULI, {}, {}, {}, pop));
                return c; }, py::arg("pauli"), "Apply rotation by π/4 about the given Pauli string (e.g., '+XIZ').")
        .def("m_pauli", [](NWQEC::Circuit &c, const std::string &p) -> NWQEC::Circuit &
             {
                if (circuit_has_non_pauli_ops(c))
                    throw std::runtime_error("Pauli-based operations are valid only in PBC circuits; do not mix with standard gates.");
                NWQEC::PauliOp pop(c.get_num_qubits()); pop.from_string(p);
                c.add_operation(NWQEC::Operation(NWQEC::Operation::Type::M_PAULI, {}, {}, {}, pop));
                return c; }, py::arg("pauli"), "Measure the given multi‑qubit Pauli string (projective measurement).")
        .def("s_pauli", [](NWQEC::Circuit &c, const std::string &p) -> NWQEC::Circuit &
             {
                if (circuit_has_non_pauli_ops(c))
                    throw std::runtime_error("Pauli-based operations are valid only in PBC circuits; do not mix with standard gates.");
                NWQEC::PauliOp pop(c.get_num_qubits()); pop.from_string(p);
                c.add_operation(NWQEC::Operation(NWQEC::Operation::Type::S_PAULI, {}, {}, {}, pop));
                return c; }, py::arg("pauli"), "Apply rotation by π/2 about the given Pauli string.")
        .def("z_pauli", [](NWQEC::Circuit &c, const std::string &p) -> NWQEC::Circuit &
             {
                if (circuit_has_non_pauli_ops(c))
                    throw std::runtime_error("Pauli-based operations are valid only in PBC circuits; do not mix with standard gates.");
                NWQEC::PauliOp pop(c.get_num_qubits()); pop.from_string(p);
                c.add_operation(NWQEC::Operation(NWQEC::Operation::Type::Z_PAULI, {}, {}, {}, pop));
                return c; }, py::arg("pauli"), "Apply rotation by π about the given Pauli string.")
        .def("num_qubits", &NWQEC::Circuit::get_num_qubits)
    .def("count_ops", &NWQEC::Circuit::count_ops)
    .def("is_clifford_t", &NWQEC::Circuit::is_clifford_t)
        .def("stats", &circuit_stats)
        .def("duration", &NWQEC::Circuit::duration, py::arg("code_distance"))
        .def("depth", &NWQEC::Circuit::depth)
        .def("to_qasm", &circuit_to_qasm)
        .def("to_qasm_str", &circuit_to_qasm)
        .def("save_qasm", &circuit_save_qasm, py::arg("path"))
        .def("to_qasm_file", &circuit_save_qasm, py::arg("filename"));

    // Module-level transforms: clean entrypoints
    m.def(
        "to_clifford_t",
        [](const NWQEC::Circuit &circuit, bool keep_ccx, py::object epsilon)
        {
            double eps_override = epsilon.is_none() ? -1.0 : epsilon.cast<double>();
            return apply_transforms(circuit,
                                    /*to_pbc=*/false,
                                    /*to_clifford_reduction=*/false,
                                    /*keep_cx=*/false,
                                    /*t_pauli_opt=*/false,
                                    /*remove_pauli=*/false,
                                    /*keep_ccx=*/keep_ccx,
                                    /*silent=*/true,
                                    /*epsilon_override=*/eps_override);
        },
        py::arg("circuit"),
        py::arg("keep_ccx") = false,
        py::arg("epsilon") = py::none(),
        "Convert the input circuit to a Clifford+T-only circuit and return a new Circuit.\n"
        "- keep_ccx: preserve CCX gates during decomposition\n"
        "- epsilon: optional absolute tolerance for RZ synthesis (applied to all angles)");

    m.def(
        "to_pbc",
        [](const NWQEC::Circuit &circuit, bool keep_cx, py::object epsilon)
        {
            double eps_override = epsilon.is_none() ? -1.0 : epsilon.cast<double>();
            return apply_transforms(circuit,
                                    /*to_pbc=*/true,
                                    /*to_clifford_reduction=*/false,
                                    /*keep_cx=*/keep_cx,
                                    /*t_pauli_opt=*/false,
                                    /*remove_pauli=*/false,
                                    /*keep_ccx=*/false,
                                    /*silent=*/true,
                                    /*epsilon_override=*/eps_override);
        },
        py::arg("circuit"),
        py::arg("keep_cx") = false,
        py::arg("epsilon") = py::none(),
        "Transpile the input circuit to a Pauli-Based Circuit (PBC) form and return a new Circuit.\n"
        "- keep_cx: preserve CX gates where possible in the PBC form\n"
        "- epsilon: optional absolute tolerance for RZ synthesis (applied to all angles)");

    m.def(
        "to_taco",
        [](const NWQEC::Circuit &circuit, py::object epsilon)
        {
            double eps_override = epsilon.is_none() ? -1.0 : epsilon.cast<double>();
            return apply_transforms(circuit,
                                    /*to_pbc=*/false,
                                    /*to_clifford_reduction=*/true,
                                    /*keep_cx=*/false,
                                    /*t_pauli_opt=*/false,
                                    /*remove_pauli=*/false,
                                    /*keep_ccx=*/false,
                                    /*silent=*/true,
                                    /*epsilon_override=*/eps_override);
        },
        py::arg("circuit"),
        py::arg("epsilon") = py::none(),
        "Apply the Clifford reduction (TACO) optimisation pipeline and return a new Circuit.\n"
        "- epsilon: optional absolute tolerance for RZ synthesis (applied to all angles)");

    // fuse_t: apply only the T-Pauli fusion stage within the PBC pipeline
    m.def(
        "fuse_t",
        [](const NWQEC::Circuit &circuit, py::object epsilon)
        {
            double eps_override = epsilon.is_none() ? -1.0 : epsilon.cast<double>();
            return apply_transforms(circuit,
                                    /*to_pbc=*/false,
                                    /*to_clifford_reduction=*/false,
                                    /*keep_cx=*/false,
                                    /*t_pauli_opt=*/true,
                                    /*remove_pauli=*/false,
                                    /*keep_ccx=*/false,
                                    /*silent=*/true,
                                    /*epsilon_override=*/eps_override);
        },
        py::arg("circuit"),
        py::arg("epsilon") = py::none(),
        "Optimize the number of T rotations within a Pauli-Based Circuit (PBC) and return a new Circuit.\n"
        "- epsilon: optional absolute tolerance for any RZ synthesis still required");

    m.def("load_qasm", [](const std::string &filename)
          {
        NWQEC::QASMParser p;
        if (!p.parse_file(filename))
        {
            throw std::runtime_error("Failed to parse QASM: " + p.get_error_message());
        }
        return p.get_circuit(); }, py::arg("filename"));
}
