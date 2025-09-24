// Python bindings for NWQEC using pybind11

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "nwqec/parser/qasm_parser.hpp"
#include "nwqec/core/pass_manager.hpp"
#include "nwqec/core/operation.hpp"
#include "nwqec/core/pauli_op.hpp"
#include "nwqec/core/constants.hpp"

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

    // Internal helper to run transforms with optional Python RZ synthesis fallback
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
        NWQEC::PassManager pm;
        auto up = std::make_unique<NWQEC::Circuit>(circuit); // copy for ownership
        auto out = pm.apply_passes(std::move(up), to_pbc, to_clifford_reduction, keep_cx, t_pauli_opt, remove_pauli, keep_ccx, silent, epsilon_override);

#if !NWQEC_WITH_GRIDSYNTH_CPP
        // Attempt transparent Python fallback for RZ synthesis if needed
        bool has_rz = false;
        for (const auto &op : out->get_operations())
        {
            if (op.get_type() == NWQEC::Operation::Type::RZ)
            {
                has_rz = true;
                break;
            }
        }
        if (has_rz)
        {
            try
            {
                // Defaults for fallback
                const int dps = NWQEC::DEFAULT_MPMATH_DPS;
                const char *module_name = "pygridsynth.gridsynth";

                py::module_ mp = py::module_::import("mpmath");
                mp.attr("mp").attr("dps") = dps; // precision
                py::object mpmathify = mp.attr("mpmathify");
                py::module_ mod = py::module_::import(module_name);
                if (!py::hasattr(mod, "gridsynth_gates"))
                {
                    throw std::runtime_error("pygridsynth module missing 'gridsynth_gates'");
                }
                py::object gridsynth_gates = mod.attr("gridsynth_gates");

                auto new_circuit = std::make_unique<NWQEC::Circuit>();
                new_circuit->add_qreg("q", out->get_num_qubits());
                new_circuit->add_creg("c", out->get_num_bits());

                const auto &ops = out->get_operations();
                for (const auto &op : ops)
                {
                    if (op.get_type() != NWQEC::Operation::Type::RZ)
                    {
                        new_circuit->add_operation(op);
                        continue;
                    }
                    const auto &params = op.get_parameters();
                    if (params.empty())
                        continue;
                    std::string theta_str = std::to_string(params[0]);
                    // Determine epsilon per angle: default to 2 orders smaller than theta, or use override if provided
                    double eps_val = (epsilon_override >= 0.0) ? epsilon_override : std::abs(params[0]) * NWQEC::DEFAULT_EPSILON_MULTIPLIER;
                    std::ostringstream eps_ss;
                    eps_ss.setf(std::ios::scientific);
                    eps_ss << std::setprecision(16) << eps_val;
                    std::string eps_str = eps_ss.str();
                    py::object theta = mpmathify(theta_str);
                    py::object epsilon = mpmathify(eps_str);
                    std::string gates = py::cast<std::string>(gridsynth_gates(theta, epsilon));
                    const auto &qs = op.get_qubits();
                    for (const char &g : gates)
                    {
                        switch (g)
                        {
                        case 'X':
                            new_circuit->add_operation(NWQEC::Operation(NWQEC::Operation::Type::X, qs));
                            break;
                        case 'Y':
                            new_circuit->add_operation(NWQEC::Operation(NWQEC::Operation::Type::Y, qs));
                            break;
                        case 'Z':
                            new_circuit->add_operation(NWQEC::Operation(NWQEC::Operation::Type::Z, qs));
                            break;
                        case 'H':
                            new_circuit->add_operation(NWQEC::Operation(NWQEC::Operation::Type::H, qs));
                            break;
                        case 'S':
                            new_circuit->add_operation(NWQEC::Operation(NWQEC::Operation::Type::S, qs));
                            break;
                        case 'T':
                            new_circuit->add_operation(NWQEC::Operation(NWQEC::Operation::Type::T, qs));
                            break;
                        case 'W':
                            break;
                        default:
                            throw std::runtime_error(std::string("Unknown gate from pygridsynth: ") + g);
                        }
                    }
                }
                out = std::move(new_circuit);
            }
            catch (const std::exception &)
            {
                // Neither C++ gridsynth nor pygridsynth available
                std::cerr << "RZ synthesis not available. Install GMP+MPFR and reinstall the module, or `pip install pygridsynth mpmath`.\n";
            }
        }
#endif

        return out;
    }
}

PYBIND11_MODULE(_core, m)
{
    m.doc() = "NWQEC Python bindings";

#if NWQEC_WITH_GRIDSYNTH_CPP
    m.attr("WITH_GRIDSYNTH_CPP") = py::bool_(true);
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
