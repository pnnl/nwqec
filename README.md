NWQEC
=====

## Overview
NWQEC is a C++/Python toolkit for fault-tolerant quantum circuit transpilation. Key capabilities include:
- Parsing OpenQASM 2.0 circuits into an internal circuit representation.
- Converting arbitrary circuits to Clifford+T using the gridsynth algorithm [1].
- Producing Pauli-Based Circuits (PBC) with Tfuse optimisation for T-count reduction [2,3].
- Applying the TACO (Clifford-reduction) pipeline [4], which preserves circuit parallelism while reducing non-T overhead.
- Leveraging a tableau-based IR for high-performance PBC passes [2].

## Requirements
- Python 3.8+
- C++17 toolchain and CMake ≥ 3.16
- Python build dependencies: `scikit-build-core`, `pybind11`
- Recommended: GMP and MPFR for the C++ gridsynth backend
- Fallback mode (optional): `pygridsynth` and `mpmath` when installing without GMP/MPFR

## Installing GMP/MPFR
- macOS (Homebrew): `brew install gmp mpfr`
- Ubuntu/Debian: `sudo apt-get install -y libgmp-dev libmpfr-dev`
- Fedora: `sudo dnf install gmp-devel mpfr-devel`
- Windows (MSYS2): `pacman -S mingw-w64-x86_64-gmp mingw-w64-x86_64-mpfr`

Note: With GMP and MPFR, the native C++ gridsynth backend can significantly accelerate Clifford+T synthesis (20× speed-up on an 18-qubit QFT benchmark).

## Installing NWQEC
```bash
pip install -U pip
pip install scikit-build-core pybind11
pip install .
```
If GMP/MPFR are not available the build will stop with a highlighted message. Either install the libraries (see above) or opt into the Python fallback:
```bash
pip install . --config-settings=cmake.define.NWQEC_ALLOW_NO_GMP=ON
pip install pygridsynth mpmath  # required for fallback RZ synthesis
```
Install GMP/MPFR later and reinstall without the flag to restore the C++ gridsynth backend.


## Quick Start (Python)
```python
import nwqec

circuit = nwqec.load_qasm("example_circuits/qft_n18.qasm")
print(circuit.stats())

clifford_t = nwqec.to_clifford_t(circuit, keep_ccx=False, epsilon=1e-10)
print("WITH_GRIDSYNTH_CPP:", nwqec.WITH_GRIDSYNTH_CPP)
print("Clifford+T gate counts:", clifford_t.count_ops())

# Pauli-Based Circuit + Tfuse
pbc = nwqec.to_pbc(circuit)
pbc_opt = nwqec.fuse_t(pbc)
print("T count before Tfuse:", pbc.count_ops().get("t_pauli", 0))
print("T count after Tfuse:", pbc_opt.count_ops().get("t_pauli", 0))

# Clifford reduction (TACO) pipeline
taco = nwqec.to_taco(circuit)
print(f"Depth reduction of TACO over Clifford+T: {clifford_t.depth() / taco.depth():.2f}x")
print(f"Depth reduction of TACO over PBC: {pbc.depth() / taco.depth():.2f}x")
```

## Further Documentation
- Python API reference: [docs/python_api.md](docs/python_api.md)
- C++ CLI guide: [docs/cpp_cli.md](docs/cpp_cli.md)

## Repository Layout
- `include/nwqec/` — public headers (core, parser, passes, gridsynth, tableau)
- `python/nwqec/` — Python package and pybind11 bindings
- `tools/` — C++ command-line utilities (`transpiler`, `gridsynth`)
- `docs/` — additional documentation
- `tests/` — Python tests

## References
1. Neil J. Ross, and Peter Selinger. "Optimal ancilla-free Clifford+ T approximation of z-rotations." arXiv preprint arXiv:1403.2975 (2014).
2. Meng Wang, Chenxu Liu, Sean Garner, Samuel Stein, Yufei Ding, Prashant J. Nair, and Ang Li. "Tableau-Based Framework for Efficient Logical Quantum Compilation." arXiv preprint arXiv:2509.02721 (2025).
3. Sean Garner, Chenxu Liu, Meng Wang, Samuel Stein, and Ang Li. "STABSim: A Parallelized Clifford Simulator with Features Beyond Direct Simulation." arXiv preprint arXiv:2507.03092 (2025).
4. Meng Wang, Chenxu Liu, Samuel Stein, Yufei Ding, Poulami Das, Prashant J. Nair, and Ang Li. "Optimizing FTQC Programs through QEC Transpiler and Architecture Codesign." arXiv preprint arXiv:2412.15434 (2024).

## Citation format

Please cite our papers:
 - Meng Wang, Chenxu Liu, Samuel Stein, Yufei Ding, Poulami Das, Prashant J. Nair, and Ang Li. "Optimizing FTQC Programs through QEC Transpiler and Architecture Codesign." arXiv preprint arXiv:2412.15434 (2024).
 - Meng Wang, Chenxu Liu, Sean Garner, Samuel Stein, Yufei Ding, Prashant J. Nair, and Ang Li. "Tableau-Based Framework for Efficient Logical Quantum Compilation." arXiv preprint arXiv:2509.02721 (2025).

Bibtex:
```text
@article{wang2024optimizing,
  title={Optimizing FTQC Programs through QEC Transpiler and Architecture Codesign},
  author={Wang, Meng and Liu, Chenxu and Stein, Samuel and Ding, Yufei and Das, Poulami and Nair, Prashant J and Li, Ang},
  journal={arXiv preprint arXiv:2412.15434},
  year={2024}
}
@article{wang2025tableau,
  title={Tableau-Based Framework for Efficient Logical Quantum Compilation},
  author={Wang, Meng and Liu, Chenxu and Garner, Sean and Stein, Samuel and Ding, Yufei and Nair, Prashant J and Li, Ang},
  journal={arXiv preprint arXiv:2509.02721},
  year={2025}
}
``` 


## Acknowledgments

**PNNL-IPID: 33474-E, iEdison No. 0685901-25-0171**

The development of this software was supported by the U.S. Department of Energy, Office of Science, National Quantum Information Science Research Centers, Co-design Center for Quantum Advantage (C2QA) under contract number DE-SC0012704, (Basic Energy Sciences, PNNL FWP 76274). It was supported by the U.S. Department of Energy, Office of Science, National Quantum Information Science Research Centers, Quantum Science Center (QSC). The development used resources of the National Energy Research Scientific Computing Center (NERSC), a U.S. Department of Energy Office of Science User Facility located at Lawrence Berkeley National Laboratory, operated under Contract No. DE-AC02-05CH11231. The development used resources of the Oak Ridge Leadership Computing Facility, which is a DOE Office of Science User Facility supported under Contract DE-AC05-00OR22725. This software was also supported by the Quantum Algorithms and Architecture for Domain Science Initiative (QuAADS), under the Laboratory Directed Research and Development (LDRD) Program at Pacific Northwest National Laboratory (PNNL). PNNL is a multi-program national laboratory operated for the U.S. Department of Energy (DOE) by Battelle Memorial Institute under Contract No. DE-AC05-76RL01830.
