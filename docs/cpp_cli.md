C++ CLI Guide
=============

Overview
--------
The repository ships two C++ command-line tools build by default:
- `transpiler`: parse OpenQASM, transpile to Clifford+T or PBC, optionally optimize T rotations, and export QASM/statistics.
- `gridsynth`: synthesize a single RZ angle into a Clifford+T sequence (requires GMP/MPFR).

Build Requirements
------------------
- CMake ≥ 3.16 and a C++17 compiler
- GMP and MPFR

Building
--------
```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```
Optional flags:
- `-DNWQEC_ENABLE_LTO=ON|OFF` (default ON)
- `-DNWQEC_ENABLE_NATIVE=ON|OFF` (default OFF)
- `-DNWQEC_BUILD_PYTHON=ON|OFF`

Usage
-----
`transpiler --help` prints the full command list. Frequently used combinations:
```
# Transpile a QASM file to Clifford+T and save the result
transpiler --input qft_n4.qasm --clifford-t --output qft_n4_clifford_t.qasm

# Pauli-Based Circuit path with Tfuse
transpiler --input circuit.qasm --pbc --opt-t --no-save

# Clifford-reduction (TACO) path
transpiler --input circuit.qasm --cr --no-save

# Generate a benchmark circuit (QFT) without a source file
transpiler --qft 4 --clifford-t --no-save
```
Mapping between CLI flags and Python helpers:
- `--clifford-t` ↔ `nwqec.to_clifford_t`
- `--pbc` ↔ `nwqec.to_pbc`
- `--cr` ↔ `nwqec.to_taco`
- `--opt-t` ↔ `nwqec.fuse_t`

`gridsynth` requires the native backend:
```
gridsynth pi/8 12  # target angle, precision bits
```

Installation Paths
------------------
`cmake --build build --target install` installs the CLI binaries into `${CMAKE_INSTALL_BINDIR}` (typically `/usr/local/bin`) and exports CMake targets under `NWQEC::` for downstream C++ projects.
