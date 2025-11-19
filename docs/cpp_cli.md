C++ CLI Guide
=============

Overview
--------
The repository ships two C++ command-line tools:
- `transpiler`: parse OpenQASM, transpile to Clifford+T or PBC, optionally optimize T rotations, and export QASM/statistics.
- `gridsynth`: synthesize a single RZ angle into a Clifford+T sequence (requires GMP/MPFR).

**Platform Support:**
- **macOS/Linux**: Both tools available with automatic prebuilt GMP/MPFR download
- **Windows**: Only `transpiler` available; users should install `pygridsynth` for RZ synthesis

Build Requirements
------------------
- CMake ≥ 3.16 and a C++17 compiler
- **macOS/Linux**: GMP/MPFR automatically downloaded from prebuilt binaries
- **Windows**: No additional dependencies (C++ gridsynth disabled, use pygridsynth instead)

Building
--------
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Optional CMake flags:
- `-DNWQEC_ENABLE_LTO=ON|OFF` (default ON) - Enable link-time optimization
- `-DNWQEC_ENABLE_NATIVE=ON|OFF` (default OFF) - Enable -march=native optimization
- `-DNWQEC_BUILD_PYTHON=ON|OFF` (default ON) - Build Python bindings

Transpiler Usage
----------------
Basic syntax: `transpiler [OPTIONS] <INPUT>`

Get help: `transpiler --help` or `transpiler -h`

### Input Sources
```bash
# Parse QASM file
transpiler circuit.qasm

# Generate test circuits
transpiler --qft 4        # QFT circuit with 4 qubits  
transpiler --shor 3       # Shor test circuit for 3-bit numbers
```

### Transpilation Passes
```bash
# Default: Clifford+T conversion
transpiler circuit.qasm

# Pauli-Based Circuit (PBC) 
transpiler circuit.qasm --pbc

# Clifford Reduction (TACO)
transpiler circuit.qasm --cr  

# Restricted PBC (preserves CCX gates)
transpiler circuit.qasm --red-pbc

# PBC with T-count optimization
transpiler circuit.qasm --pbc --t-opt
```

**Note**: PBC (`--pbc`), Clifford Reduction (`--cr`), and Restricted PBC (`--red-pbc`) are mutually exclusive.

### Output Options
```bash
# Default: saves to <input>_transpiled.qasm
transpiler circuit.qasm

# Custom output filename
transpiler circuit.qasm -o my_output.qasm
transpiler circuit.qasm --output my_output.qasm

# Don't save file (display stats only)
transpiler circuit.qasm --no-save
```

### Analysis Options
```bash
# Remove Pauli gates from output
transpiler circuit.qasm --remove-pauli

# Preserve CCX gates during decomposition
transpiler circuit.qasm --keep-ccx
```

### Complete Examples
```bash
# Basic transpilation
transpiler qft_n4.qasm

# PBC with T optimization, custom output
transpiler circuit.qasm --pbc --t-opt -o optimized.qasm

# Generate QFT, apply Clifford reduction, don't save
transpiler --qft 8 --cr --no-save

# Shor circuit with restricted PBC
transpiler --shor 4 --red-pbc

# Advanced: PBC with all options
transpiler large_circuit.qasm --pbc --t-opt --remove-pauli --keep-ccx
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

**Windows Users**: Install pygridsynth instead:
```bash
pip install pygridsynth
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