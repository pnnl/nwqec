Python API Reference
====================

The `nwqec` package exposes the core transpilation functionality to Python users. This document lists the available functions, classes, signatures, arguments, and return values.

Module Constants
----------------
- `WITH_GRIDSYNTH_CPP: bool`
  - `True` if the package was built with the C++ gridsynth backend using prebuilt GMP/MPFR libraries.
  - Always `True` on supported platforms (macOS/Linux) as prebuilt binaries are automatically downloaded.
- `__version__: str`
  - Package version string.

Top-Level Functions
-------------------
Keyword-only options must be supplied by name (see bullet lists). Each function returns a new `Circuit` instance and raises `RuntimeError` on failure.

- **`load_qasm(path: str) -> Circuit`**
  - `path`: filesystem path to an OpenQASM 2.0 file.
  - Parses the file into a circuit.

- **`to_clifford_t(circuit: Circuit, keep_ccx: bool = False, epsilon: float | None = None) -> Circuit`**
  - `circuit`: source circuit.
  - `keep_ccx`: preserve CCX gates when `True`.
  - `epsilon`: absolute error tolerance for RZ synthesis; defaults to `abs(theta) * DEFAULT_EPSILON_MULTIPLIER` per angle.
  - Produces a Clifford+T-only circuit.

- **`to_pbc(circuit: Circuit, keep_cx: bool = False, optimize_t_count: bool = False, epsilon: float | None = None) -> Circuit`**
  - `circuit`: source circuit.
  - `keep_cx`: preserve CX gates where possible in the PBC form.
  - `optimize_t_count`: apply T-count optimization after PBC conversion.
  - `epsilon`: absolute error tolerance for RZ synthesis.
  - Transpiles the circuit to a Pauli-Based Circuit (PBC).

- **`to_clifford_reduction(circuit: Circuit, epsilon: float | None = None) -> Circuit`**
  - `circuit`: source circuit.
  - `epsilon`: absolute error tolerance for RZ synthesis.
  - Applies the Clifford reduction optimization (preserves parallelism while reducing non-T overhead).
  - Based on the technique from Wang et al. "Optimizing FTQC Programs through QEC Transpiler and Architecture Codesign" (2024).

- **`fuse_t(circuit: Circuit, epsilon: float | None = None) -> Circuit`**
  - `circuit`: source circuit, consisting exclusively of Pauli-based operations.
  - `epsilon`: absolute error tolerance for any remaining RZ synthesis.
  - Applies the Tfuse optimisation to reduce T rotations.

Circuit Class
-------------
Create with `Circuit(num_qubits: int)`.

### Inspection
- `num_qubits() -> int`
- `count_ops() -> dict[str, int]`
- `depth() -> int`
- `stats() -> str`
- `duration(code_distance: float) -> float`
- `to_qasm_str() -> str`
- `save_qasm(path: str) -> None`
- `to_qasm_file(filename: str) -> None`  # alias for save_qasm
- `is_clifford_t() -> bool`

### Single-Qubit Gates
- `h(q: int)`
- `x(q: int)`
- `y(q: int)`
- `z(q: int)`
- `s(q: int)`
- `sdg(q: int)`
- `t(q: int)`
- `tdg(q: int)`
- `sx(q: int)`
- `sxdg(q: int)`
- `rx(q: int, theta: float)`
- `ry(q: int, theta: float)`
- `rz(q: int, theta: float)`
- `rxp(q: int, multiplier: float)`  # rotation by `multiplier * π`
- `ryp(q: int, multiplier: float)`
- `rzp(q: int, multiplier: float)`

### Multi-Qubit Gates
- `cx(control: int, target: int)`
- `cz(control: int, target: int)`
- `swap(q0: int, q1: int)`
- `ccx(c0: int, c1: int, target: int)`

### Classical/Utility
- `measure(q: int, cbit: int)`
- `reset(q: int)`
- `barrier(qubits: Sequence[int])`

### Pauli-Based Operations
Valid only for PBC circuits (mixing with standard gates raises `RuntimeError`).
- `t_pauli(pauli: str)` — rotation by π/4 about the given signed Pauli string.
- `s_pauli(pauli: str)` — rotation by π/2 about the given signed Pauli string.
- `z_pauli(pauli: str)` — rotation by π about the given signed Pauli string.
- `m_pauli(pauli: str)` — projective measurement of the given signed Pauli string.

`pauli` strings must start with `+` or `-`, followed by one character per qubit chosen from `{X, Y, Z, I}`.


Examples
--------
```python
import nwqec

# Load and inspect a circuit
c = nwqec.load_qasm("example_circuits/qft_n18.qasm")
print(c.stats())

# Clifford+T conversion (default)
ct = nwqec.to_clifford_t(c, keep_ccx=False, epsilon=1e-10)
print("Clifford+T gate counts:", ct.count_ops())

# PBC conversion
pbc = nwqec.to_pbc(c, keep_cx=False)
print("PBC gate counts:", pbc.count_ops())

# PBC with T-count optimization (single call)
pbc_opt = nwqec.to_pbc(c, optimize_t_count=True)
print("Optimized PBC gate counts:", pbc_opt.count_ops())

# Alternative: PBC then T-optimization (two-step process)
pbc = nwqec.to_pbc(c)
pbc_opt = nwqec.fuse_t(pbc)
print("T count after separate Tfuse:", pbc_opt.count_ops().get("t_pauli", 0))

# Clifford Reduction optimization
clifford_reduced = nwqec.to_clifford_reduction(c)
print("Clifford Reduction gate counts:", clifford_reduced.count_ops())

# Save results
pbc_opt.save_qasm("example_circuits/qft_n18_transpiled.qasm")
```
