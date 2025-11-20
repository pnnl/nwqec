#include "nwqec/parser/qasm_parser.hpp"
#include "nwqec/core/transpiler.hpp"

#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <chrono>
#include <iomanip>
#include <filesystem>
#include <cmath>
#include <fstream>

// PROJECT_ROOT_DIR is defined by CMake during compilation
// This fallback is for IDE IntelliSense support only
#ifndef PROJECT_ROOT_DIR
#define PROJECT_ROOT_DIR "."
#endif

std::unique_ptr<NWQEC::Circuit> generate_qft_circuit(int n_qubits);
std::unique_ptr<NWQEC::Circuit> generate_shor_circuit(int n_bits);

int main(int argc, char *argv[])
{
    // Parse the QASM code
    NWQEC::QASMParser parser;
    bool success;
    std::unique_ptr<NWQEC::Circuit> circuit;
    std::string qasm_file;
    bool generate_qft = false;
    int qft_qubits = 0;
    bool generate_shor = false;
    int shor_bits = 0;
    bool save_to_file = true;
    std::string output_filename = "";
    bool to_pbc = false;
    bool to_clifford_reduction = false;
    bool t_pauli_opt = false;
    bool remove_pauli = false;
    bool keep_ccx = false;
    bool keep_cx = false;

    // Helper function to print usage
    auto print_usage = [&](bool detailed = false)
    {
        std::cout << "NWQEC Quantum Circuit Transpiler" << std::endl;
        std::cout << "" << std::endl;
        std::cout << "Usage: " << argv[0] << " [OPTIONS] <INPUT>" << std::endl;
        std::cout << "" << std::endl;

        if (detailed)
        {
            std::cout << "DESCRIPTION:" << std::endl;
            std::cout << "  Transpiles quantum circuits to Clifford+T representation with various" << std::endl;
            std::cout << "  optimization passes. Supports QASM input files or generates test circuits." << std::endl;
            std::cout << "" << std::endl;
        }

        std::cout << "INPUT OPTIONS:" << std::endl;
        std::cout << "  <qasm_file>           Path to QASM file to transpile" << std::endl;
        std::cout << "  --qft <n_qubits>      Generate QFT circuit with n_qubits qubits" << std::endl;
        std::cout << "  --shor <n_bits>       Generate Shor test circuit for n_bits-bit number" << std::endl;
        std::cout << "" << std::endl;

        std::cout << "TRANSPILATION OPTIONS:" << std::endl;
        std::cout << "  --pbc                 Convert to Pauli-Based Circuit (PBC) format" << std::endl;
        std::cout << "  --cr                  Apply Clifford Reduction optimization" << std::endl;
        std::cout << "  --t-opt               Apply T-count optimization (works on PBC circuits)" << std::endl;
        std::cout << "  --keep-ccx            Preserve CCX gates during Clifford+T conversion" << std::endl;
        std::cout << "  --keep-cx             Preserve CX gates during PBC conversion" << std::endl;
        std::cout << "" << std::endl;

        std::cout << "ANALYSIS OPTIONS:" << std::endl;
        std::cout << "  --remove-pauli        Remove all Pauli gates (X, Y, Z) from final circuit" << std::endl;
        std::cout << "" << std::endl;

        std::cout << "OUTPUT OPTIONS:" << std::endl;
        std::cout << "  --no-save             Don't save transpiled circuit to file" << std::endl;
        std::cout << "  -o <file>             Specify output filename for transpiled circuit" << std::endl;
        std::cout << "" << std::endl;

        std::cout << "OTHER OPTIONS:" << std::endl;
        std::cout << "  -h                    Show this help message" << std::endl;

        if (detailed)
        {
            std::cout << "" << std::endl;
            std::cout << "EXAMPLES:" << std::endl;
            std::cout << "  " << argv[0] << " circuit.qasm" << std::endl;
            std::cout << "    Convert to Clifford+T format" << std::endl;
            std::cout << "" << std::endl;
            std::cout << "  " << argv[0] << " circuit.qasm --keep-ccx" << std::endl;
            std::cout << "    Convert to Clifford+T while preserving CCX gates" << std::endl;
            std::cout << "" << std::endl;
            std::cout << "  " << argv[0] << " circuit.qasm --pbc" << std::endl;
            std::cout << "    Convert to Pauli-Based Circuit format" << std::endl;
            std::cout << "" << std::endl;
            std::cout << "  " << argv[0] << " circuit.qasm --pbc --keep-cx" << std::endl;
            std::cout << "    Convert to PBC while preserving CX gates" << std::endl;
            std::cout << "" << std::endl;
            std::cout << "  " << argv[0] << " circuit.qasm --pbc --t-opt" << std::endl;
            std::cout << "    Convert to PBC and apply T-count optimization" << std::endl;
            std::cout << "" << std::endl;
            std::cout << "  " << argv[0] << " pbc_circuit.qasm --t-opt" << std::endl;
            std::cout << "    Apply T-count optimization to existing PBC circuit" << std::endl;
            std::cout << "" << std::endl;
            std::cout << "  " << argv[0] << " circuit.qasm --cr" << std::endl;
            std::cout << "    Apply Clifford Reduction optimization" << std::endl;
            std::cout << "" << std::endl;
            std::cout << "  " << argv[0] << " circuit.qasm -o my_output.qasm" << std::endl;
            std::cout << "    Transpile and save to custom filename" << std::endl;
            std::cout << "" << std::endl;
            std::cout << "WORKFLOW:" << std::endl;
            std::cout << "  1. Basic processing: decompose → remove trivial RZ → synthesize RZ" << std::endl;
            std::cout << "  2. Choose format: Clifford+T (default), PBC, or Clifford Reduction" << std::endl;
            std::cout << "  3. Optional: T-optimization (for PBC), cleanup passes" << std::endl;
            std::cout << "" << std::endl;
            std::cout << "NOTES:" << std::endl;
            std::cout << "  - Format options (--pbc, --cr) are mutually exclusive" << std::endl;
            std::cout << "  - T-optimization (--t-opt) works on PBC circuits (combine with --pbc or use alone)" << std::endl;
            std::cout << "  - CCX preservation (--keep-ccx) applies to Clifford+T conversion" << std::endl;
            std::cout << "  - CX preservation (--keep-cx) applies to PBC conversion" << std::endl;
            std::cout << "  - Clifford Reduction (--cr) optimizes circuit depth while preserving parallelism" << std::endl;
            std::cout << "  - Output files use '_transpiled.qasm' suffix by default" << std::endl;
        }
    };

    // Use example string or file from command line
    auto start_parse = std::chrono::high_resolution_clock::now();
    if (argc < 2)
    {
        print_usage();
        return 1;
    }

    // Parse command line arguments - process all arguments to handle flags in any position
    for (int arg_index = 1; arg_index < argc; arg_index++)
    {
        std::string arg = argv[arg_index];

        // Handle help first
        if (arg == "-h")
        {
            print_usage(true);
            return 0;
        }
        else if (arg == "--qft")
        {
            if (arg_index + 1 >= argc)
            {
                std::cout << "Error: --qft requires number of qubits" << std::endl;
                std::cout << "Usage: --qft <n_qubits>" << std::endl;
                return 1;
            }
            if (generate_shor)
            {
                std::cout << "Error: Cannot specify both --qft and --shor" << std::endl;
                return 1;
            }
            generate_qft = true;
            arg_index++;
            try
            {
                qft_qubits = std::stoi(argv[arg_index]);
                if (qft_qubits <= 0)
                {
                    std::cout << "Error: number of qubits must be positive, got: " << qft_qubits << std::endl;
                    return 1;
                }
                if (qft_qubits > 20)
                {
                    std::cout << "Warning: Large QFT circuit (" << qft_qubits << " qubits) may take significant time and memory" << std::endl;
                }
            }
            catch (const std::exception &e)
            {
                std::cout << "Error: invalid number of qubits: '" << argv[arg_index] << "' (must be a positive integer)" << std::endl;
                return 1;
            }
        }
        else if (arg == "--shor")
        {
            if (arg_index + 1 >= argc)
            {
                std::cout << "Error: --shor requires number of bits" << std::endl;
                std::cout << "Usage: --shor <n_bits>" << std::endl;
                return 1;
            }
            if (generate_qft)
            {
                std::cout << "Error: Cannot specify both --qft and --shor" << std::endl;
                return 1;
            }
            generate_shor = true;
            arg_index++;
            try
            {
                shor_bits = std::stoi(argv[arg_index]);
                if (shor_bits <= 0)
                {
                    std::cout << "Error: number of bits must be positive, got: " << shor_bits << std::endl;
                    return 1;
                }
                if (shor_bits > 15)
                {
                    std::cout << "Warning: Large Shor circuit (" << shor_bits << " bits) may take significant time and memory" << std::endl;
                }
            }
            catch (const std::exception &e)
            {
                std::cout << "Error: invalid number of bits: '" << argv[arg_index] << "' (must be a positive integer)" << std::endl;
                return 1;
            }
        }
        else if (arg == "--no-save")
        {
            save_to_file = false;
            std::cout << "File saving disabled" << std::endl;
        }
        else if (arg == "-o")
        {
            if (arg_index + 1 >= argc)
            {
                std::cout << "Error: -o requires output filename" << std::endl;
                std::cout << "Usage: -o <filename>" << std::endl;
                return 1;
            }
            arg_index++;
            output_filename = argv[arg_index];
            std::cout << "Output filename set to: " << output_filename << std::endl;
        }
        else if (arg == "--pbc")
        {
            if (to_clifford_reduction)
            {
                std::cout << "Error: Cannot specify --pbc with other transpilation passes (mutually exclusive)" << std::endl;
                return 1;
            }
            to_pbc = true;
            std::cout << "PBC pass enabled" << std::endl;
        }
        else if (arg == "--cr")
        {
            if (to_pbc)
            {
                std::cout << "Error: Cannot specify --cr with other transpilation passes (mutually exclusive)" << std::endl;
                return 1;
            }
            to_clifford_reduction = true;
            std::cout << "Clifford Reduction pass enabled" << std::endl;
        }
        else if (arg == "--t-opt")
        {
            t_pauli_opt = true;
            std::cout << "T Pauli optimizer enabled" << std::endl;
        }
        else if (arg == "--remove-pauli")
        {
            remove_pauli = true;
            std::cout << "Pauli gate removal enabled" << std::endl;
        }
        else if (arg == "--keep-ccx")
        {
            keep_ccx = true;
            std::cout << "CCX gate preservation enabled" << std::endl;
        }
        else if (arg == "--keep-cx")
        {
            keep_cx = true;
            std::cout << "CX gate preservation enabled" << std::endl;
        }
        else if (arg[0] == '-')
        {
            std::cout << "Error: unknown option '" << arg << "'" << std::endl;
            std::cout << "" << std::endl;
            print_usage(true);
            return 1;
        }
        else
        {
            // This is a non-flag argument (QASM file)
            if (!generate_qft && !generate_shor && qasm_file.empty())
            {
                // Check if file exists (for better error messages)
                std::ifstream file_check(arg);
                if (!file_check.good() && !arg.empty())
                {
                    std::cout << "Warning: File '" << arg << "' does not exist or is not readable" << std::endl;
                }
                qasm_file = arg;
            }
            else if (!generate_qft && !generate_shor)
            {
                std::cout << "Error: multiple input files specified" << std::endl;
                std::cout << "Current file: '" << qasm_file << "', additional file: '" << arg << "'" << std::endl;
                return 1;
            }
            else
            {
                std::cout << "Error: cannot specify both generated circuit and input file" << std::endl;
                return 1;
            }
        }
    }

    // Validate arguments after parsing
    if (!generate_qft && !generate_shor && qasm_file.empty())
    {
        std::cout << "Error: No input specified" << std::endl;
        std::cout << "Please provide a QASM file, or use --qft or --shor to generate a test circuit" << std::endl;
        std::cout << "Use -h for more information" << std::endl;
        return 1;
    }

    // T-optimization can be used alone (for existing PBC circuits) or with --pbc

    // Validate that CX preservation is only used with PBC
    if (keep_cx && !to_pbc)
    {
        std::cout << "Error: CX gate preservation (--keep-cx) requires PBC pass (--pbc)" << std::endl;
        std::cout << "Please add --pbc flag when using --keep-cx" << std::endl;
        return 1;
    }

    // Show configuration summary
    {
        std::cout << "\n==== Configuration Summary ====" << std::endl;
        if (generate_qft)
        {
            std::cout << "Input: QFT circuit (" << qft_qubits << " qubits)" << std::endl;
        }
        else if (generate_shor)
        {
            std::cout << "Input: Shor circuit (" << shor_bits << " bits)" << std::endl;
        }
        else
        {
            std::cout << "Input: " << qasm_file << std::endl;
        }

        std::cout << "Passes: ";
        std::vector<std::string> passes;
        if (to_pbc)
            passes.push_back("PBC");
        if (to_clifford_reduction)
            passes.push_back("Clifford Reduction");
        if (t_pauli_opt)
            passes.push_back("T Optimization");
        if (passes.empty())
            passes.push_back("Standard Clifford+T");

        for (size_t i = 0; i < passes.size(); ++i)
        {
            std::cout << passes[i];
            if (i < passes.size() - 1)
                std::cout << ", ";
        }
        std::cout << std::endl;

        std::vector<std::string> options;
        if (remove_pauli)
            options.push_back("Remove Pauli gates from final circuit");
        if (keep_ccx)
            options.push_back("CCX gate preservation enabled");
        if (keep_cx)
            options.push_back("CX gate preservation enabled");
        
        if (!options.empty()) {
            std::cout << "Options: ";
            for (size_t i = 0; i < options.size(); ++i) {
                std::cout << options[i];
                if (i < options.size() - 1)
                    std::cout << ", ";
            }
            std::cout << std::endl;
        }
        if (!save_to_file)
            std::cout << "Output: No file will be saved" << std::endl;
        else if (!output_filename.empty())
            std::cout << "Output: Custom filename: " << output_filename << std::endl;
        std::cout << "==============================\n"
                  << std::endl;
    }

    if (generate_qft)
    {
        // Generate QFT circuit
        std::cout << "Generating QFT circuit with " << qft_qubits << " qubits..." << std::endl;
        circuit = generate_qft_circuit(qft_qubits);
        success = true; // QFT generation always succeeds
    }
    else if (generate_shor)
    {
        // Generate Shor test circuit
        std::cout << "Generating Shor test circuit for " << shor_bits << "-bit number..." << std::endl;
        circuit = generate_shor_circuit(shor_bits);
        success = true; // Shor generation always succeeds
    }
    else
    {
        // Parse QASM file
        success = parser.parse_file(qasm_file);
        if (success)
        {
            circuit = parser.get_circuit();
        }
    }

    auto end_parse = std::chrono::high_resolution_clock::now();
    auto parse_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_parse - start_parse).count();

    if (!success)
    {
        if (generate_qft)
        {
            std::cerr << "Error: Failed to generate QFT circuit with " << qft_qubits << " qubits" << std::endl;
            std::cerr << "This may be due to insufficient memory or invalid parameters" << std::endl;
        }
        else if (generate_shor)
        {
            std::cerr << "Error: Failed to generate Shor test circuit for " << shor_bits << " bits" << std::endl;
            std::cerr << "This may be due to insufficient memory or invalid parameters" << std::endl;
        }
        else
        {
            std::cerr << "Error: Failed to parse QASM file '" << qasm_file << "'" << std::endl;
            std::cerr << "Parser error: " << parser.get_error_message() << std::endl;
            std::cerr << "Please check that the file exists and contains valid QASM code" << std::endl;
        }
        return 1;
    }

    auto start_circuit = std::chrono::high_resolution_clock::now();
    // circuit is already obtained above
    auto end_circuit = std::chrono::high_resolution_clock::now();
    auto circuit_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_circuit - start_circuit).count();

    std::cout << "Circuit contains " << circuit->get_num_qubits() << " qubits, "
              << circuit->get_num_bits() << " classical bits, and "
              << circuit->get_operations().size() << " operations." << std::endl;

    {
        circuit->print_stats(std::cout);
    }

    // Apply transpilation passes
    auto start_transpile = std::chrono::high_resolution_clock::now();
    try
    {
        NWQEC::Transpiler transpiler;
        NWQEC::PassConfig config;
        config.keep_ccx = keep_ccx;
        config.keep_cx = keep_cx;
        config.silent = false;  // CLI always shows output
        
        // Choose the appropriate pass sequence based on the logical workflow
        std::vector<NWQEC::PassType> passes;
        
        if (to_pbc && t_pauli_opt) {
            // PBC with T-optimization: full pipeline
            passes = NWQEC::PassSequences::TO_PBC_OPTIMIZED;
        } else if (to_pbc) {
            // PBC only
            passes = NWQEC::PassSequences::TO_PBC;
        } else if (to_clifford_reduction) {
            // Clifford reduction optimization pipeline
            passes = NWQEC::PassSequences::TO_CLIFFORD_REDUCTION;
        } else if (t_pauli_opt) {
            // T-optimization on existing PBC circuit
            passes = NWQEC::PassSequences::T_OPTIMIZATION_ONLY;
        } else {
            // Default: Clifford+T conversion
            passes = NWQEC::PassSequences::TO_CLIFFORD_T;
        }
        
        // Add cleanup if requested
        if (remove_pauli) {
            passes.push_back(NWQEC::PassType::REMOVE_PAULI);
        }
        
        circuit = transpiler.execute_passes(std::move(circuit), passes, config);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error during transpilation: " << e.what() << std::endl;
        std::cerr << "This may be due to:" << std::endl;
        
        std::cerr << "  - Invalid circuit structure" << std::endl;
        std::cerr << "  - Insufficient memory for large circuits" << std::endl;
        return 1;
    }
    auto end_transpile = std::chrono::high_resolution_clock::now();
    auto transpile_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_transpile - start_transpile).count();

    circuit->print_stats(std::cout);

    // save the circuit to a file with the name of the input file + _transpiled.qasm
    auto start_save = std::chrono::high_resolution_clock::now();
    std::string filename;
    if (save_to_file)
    {
        if (!output_filename.empty())
        {
            filename = output_filename;
        }
        else if (generate_qft)
        {
            filename = "qft_n" + std::to_string(qft_qubits) + "_transpiled.qasm";
        }
        else if (generate_shor)
        {
            filename = "shor_n" + std::to_string(shor_bits) + "_transpiled.qasm";
        }
        else
        {
            filename = qasm_file;
            std::string post_fix = "_transpiled.qasm";
            size_t dot_pos = filename.find_last_of('.');
            if (dot_pos != std::string::npos)
            {
                filename = filename.substr(0, dot_pos) + post_fix;
            }
            else
            {
                filename = filename + post_fix;
            }
        }
        {
            std::ofstream outfile(filename);
            circuit->print(outfile);
            outfile.close();
            std::cout << "Saved transpiled circuit to: " << filename << std::endl;
        }
    }

    auto end_save = std::chrono::high_resolution_clock::now();
    auto save_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_save - start_save).count();

    // Print timing results
    // Print timing results
    std::cout << "\n---- Performance Metrics ----" << std::endl;
    std::cout << std::left << std::setw(20) << "Parsing time:"
              << std::right << std::setw(5) << parse_time << " ms" << std::endl;
    std::cout << std::left << std::setw(20) << "Conversion time:"
              << std::right << std::setw(5) << circuit_time << " ms" << std::endl;
    std::cout << std::left << std::setw(20) << "Transpilation time:"
              << std::right << std::setw(5) << transpile_time << " ms" << std::endl;
    std::cout << std::left << std::setw(20) << "Write to file:"
              << std::right << std::setw(5) << save_time << " ms" << std::endl;
    std::cout << std::left << std::setw(20) << "Total time:"
              << std::right << std::setw(5) << (parse_time + circuit_time + transpile_time + save_time) << " ms" << std::endl;

    return 0;
}

// Function to generate a QFT circuit for n qubits
std::unique_ptr<NWQEC::Circuit> generate_qft_circuit(int n_qubits)
{
    auto circuit = std::make_unique<NWQEC::Circuit>();

    // Add quantum register
    circuit->add_qreg("q", n_qubits);

    // Generate QFT gates
    for (int i = 0; i < n_qubits; i++)
    {
        // Apply Hadamard gate
        circuit->add_operation(NWQEC::Operation(NWQEC::Operation::Type::H, {static_cast<size_t>(i)}));

        // Apply controlled phase rotations
        for (int j = i + 1; j < n_qubits; j++)
        {
            double angle = M_PI / std::pow(2, j - i);
            circuit->add_operation(NWQEC::Operation(NWQEC::Operation::Type::CP,
                                                    {static_cast<size_t>(j), static_cast<size_t>(i)}, {angle}));
        }
    }

    // Swap qubits to reverse the order
    for (int i = 0; i < n_qubits / 2; i++)
    {
        circuit->add_operation(NWQEC::Operation(NWQEC::Operation::Type::SWAP,
                                                {static_cast<size_t>(i), static_cast<size_t>(n_qubits - 1 - i)}));
    }

    return circuit;
}

// Function to generate a Shor test circuit for n-bit numbers
std::unique_ptr<NWQEC::Circuit> generate_shor_circuit(int n_bits)
{
    auto circuit = std::make_unique<NWQEC::Circuit>();

    // Calculate qubits and gates based on Shor's algorithm scaling
    int num_qubits = static_cast<int>(3 * n_bits + 0.002 * n_bits * std::log2(n_bits));
    int num_toffolis = static_cast<int>(0.3 * std::pow(n_bits, 3) + 0.0005 * std::pow(n_bits, 3) * std::log2(n_bits));

    // Add quantum register
    circuit->add_qreg("q", num_qubits);

    // Add some initialization gates (random H gates)
    int init_gates = std::min(num_qubits / 4, 10);
    std::srand(42); // Fixed seed for reproducible results
    for (int i = 0; i < init_gates; i++)
    {
        int qubit = std::rand() % num_qubits;
        circuit->add_operation(NWQEC::Operation(NWQEC::Operation::Type::H, {static_cast<size_t>(qubit)}));
    }

    // Generate random Toffoli gates
    for (int i = 0; i < num_toffolis; i++)
    {
        // Generate three unique random qubits
        std::vector<int> qubits;
        while (qubits.size() < 3)
        {
            int qubit = std::rand() % num_qubits;
            if (std::find(qubits.begin(), qubits.end(), qubit) == qubits.end())
            {
                qubits.push_back(qubit);
            }
        }

        // Add CCX gate
        circuit->add_operation(NWQEC::Operation(NWQEC::Operation::Type::CCX,
                                                {static_cast<size_t>(qubits[0]),
                                                 static_cast<size_t>(qubits[1]),
                                                 static_cast<size_t>(qubits[2])}));
    }

    // Add inverse QFT before measurement
    for (int i = 0; i < num_qubits / 2; i++)
    {
        circuit->add_operation(NWQEC::Operation(NWQEC::Operation::Type::SWAP,
                                                {static_cast<size_t>(i), static_cast<size_t>(num_qubits - 1 - i)}));
    }

    // Apply inverse QFT gates (reverse order of forward QFT)
    for (int i = num_qubits - 1; i >= 0; i--)
    {
        // Apply inverse controlled phase rotations
        for (int j = num_qubits - 1; j > i; j--)
        {
            double angle = -M_PI / std::pow(2, j - i);
            circuit->add_operation(NWQEC::Operation(NWQEC::Operation::Type::CP,
                                                    {static_cast<size_t>(j), static_cast<size_t>(i)}, {angle}));
        }

        // Apply Hadamard gate
        circuit->add_operation(NWQEC::Operation(NWQEC::Operation::Type::H, {static_cast<size_t>(i)}));
    }

    return circuit;
}
