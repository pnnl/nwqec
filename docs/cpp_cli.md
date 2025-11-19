C++ CLI Guide
=============

Overview
--------
The repository provides two C++ command-line tools:
- `nwqec-cli`: parse OpenQASM, transpile to Clifford+T or PBC, optionally optimize T rotations, and export QASM/statistics.
- `gridsynth`: synthesize a single RZ angle into a Clifford+T sequence.

**Platform Support:**
NWQEC is supported on macOS and Linux. Both tools are available with automatic prebuilt GMP/MPFR download.

Build Requirements
------------------
- CMake ≥ 3.16 and a C++17 compiler
- GMP/MPFR automatically downloaded from prebuilt binaries

Building
--------
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Optional CMake flags:
- `-DNWQEC_ENABLE_LTO=ON|OFF` (default ON) - Enable link-time optimization
- `-DNWQEC_ENABLE_NATIVE=ON|OFF` (default OFF) - Enable -march=native optimization
- `-DNWQEC_BUILD_PYTHON=ON|OFF` (default OFF) - Build Python bindings

CLI Usage
---------
Basic syntax: `nwqec-cli [OPTIONS] <INPUT>`

Get help: `nwqec-cli --help` or `nwqec-cli -h`

### Input Sources
```bash
# Parse QASM file
nwqec-cli circuit.qasm

# Generate test circuits
nwqec-cli --qft 4        # QFT circuit with 4 qubits  
nwqec-cli --shor 3       # Shor test circuit for 3-bit numbers
```

##### Transpilation Passes
```bash
# Default: Clifford+T conversion
nwqec-cli circuit.qasm

# Pauli-Based Circuit (PBC) 
nwqec-cli circuit.qasm --pbc

# Clifford Reduction (TACO)
nwqec-cli circuit.qasm --cr  

# Restricted PBC (preserves CCX gates)
nwqec-cli circuit.qasm --red-pbc

# PBC with T-count optimization
nwqec-cli circuit.qasm --pbc --t-opt
```

**Note**: PBC (`--pbc`), Clifford Reduction (`--cr`), and Restricted PBC (`--red-pbc`) are mutually exclusive.

### Output Options
```bash
# Default: saves to <input>_transpiled.qasm
nwqec-cli circuit.qasm

# Custom output filename
nwqec-cli circuit.qasm -o my_output.qasm
nwqec-cli circuit.qasm --output my_output.qasm

# Don't save file (display stats only)
nwqec-cli circuit.qasm --no-save
```

### Analysis Options
```bash
# Remove Pauli gates from output
nwqec-cli circuit.qasm --remove-pauli

# Preserve CCX gates during decomposition
nwqec-cli circuit.qasm --keep-ccx
```

### Complete Examples
```bash
# Basic transpilation
nwqec-cli qft_n4.qasm

# PBC with T optimization, custom output
nwqec-cli circuit.qasm --pbc --t-opt -o optimized.qasm

# Generate QFT, apply Clifford reduction, don't save
nwqec-cli --qft 8 --cr --no-save

# Shor circuit with restricted PBC
nwqec-cli --shor 4 --red-pbc

# Advanced: PBC with all options
nwqec-cli large_circuit.qasm --pbc --t-opt --remove-pauli --keep-ccx
```

Gridsynth Usage
---------------
**Available on macOS/Linux only**

Syntax: `gridsynth <angle> <precision_bits>`

```bash
# Synthesize π/8 rotation with 12 bits of precision
gridsynth pi/8 12

# Synthesize π/4 rotation with 10 bits of precision  
gridsynth pi/4 10

# Synthesize arbitrary angle
gridsynth 0.785398 15  # approximately π/4
```


Performance Notes
-----------------
- **Large circuits**: QFT >20 qubits or Shor >15 bits may require significant time/memory
- **T optimization**: `--t-opt` can substantially reduce T-count but increases computation time
- **Timing metrics**: The tool reports parsing, transpilation, and file I/O times

Output Format
-------------
- **QASM files**: Standard OpenQASM 2.0 format
- **Statistics**: Gate counts, circuit depth, T-count, and performance metrics
- **Default naming**: Input file with `_transpiled.qasm` suffix

Installation
------------
After building:
```bash
cmake --build build --target install
```

This installs CLI binaries to `${CMAKE_INSTALL_BINDIR}` (typically `/usr/local/bin`) and exports CMake targets under `NWQEC::` namespace for downstream C++ projects.