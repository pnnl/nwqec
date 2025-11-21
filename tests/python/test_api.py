import importlib.util
from pathlib import Path
import pytest


def _require_nwqec():
    spec = importlib.util.find_spec("nwqec")
    if spec is None:
        pytest.skip("nwqec module not installed")
    return spec


def test_import_nwqec():
    _require_nwqec()


def test_load_and_stats(tmp_path):
    _require_nwqec()
    import nwqec

    repo_root = Path(__file__).resolve().parents[2]
    qasm_path = repo_root / "tests" / "python" / "fixtures" / "fixture_circuit.qasm"
    assert qasm_path.exists()

    circuit = nwqec.load_qasm(str(qasm_path))
    assert circuit.num_qubits() > 0
    stats = circuit.stats()
    assert "Circuit Statistics" in stats
    counts = circuit.count_ops()
    assert isinstance(counts, dict)

    clifford = nwqec.to_clifford_t(circuit, epsilon=1e-6)
    assert clifford.is_clifford_t()

    pbc = nwqec.to_pbc(circuit, epsilon=1e-6)
    fused = nwqec.fuse_t(pbc)
    assert fused.count_ops().get("t_pauli", 0) <= pbc.count_ops().get("t_pauli", 0)

    out_file = tmp_path / "out.qasm"
    clifford.to_qasm_file(str(out_file))
    assert out_file.exists()
