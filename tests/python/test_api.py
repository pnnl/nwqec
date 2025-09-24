import os
import importlib.util
import pytest


def test_import_nwqec():
    assert importlib.util.find_spec("nwqec") is not None


@pytest.mark.skipif(importlib.util.find_spec("nwqec") is None, reason="nwqec module not installed")
def test_load_and_stats():
    import nwqec
    # Use a small benchmark file
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    qasm_path = os.path.join(repo_root, "example_circuits", "qft_n18.qasm")
    assert os.path.exists(qasm_path)

    c = nwqec.load_qasm(qasm_path)
    # Basic sanity checks
    assert c.num_qubits() > 0
    s = c.stats()
    assert "Circuit Statistics" in s
    counts = c.count_ops()
    assert isinstance(counts, dict)
