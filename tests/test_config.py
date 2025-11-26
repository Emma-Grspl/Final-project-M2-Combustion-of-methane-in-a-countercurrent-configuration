from counterflow.config import default_params

def test_default_params():
    p = default_params()
    assert "Nx" in p
    assert "Ny" in p
    assert "dt" in p
    assert p["Nx"] > 0
    assert p["dt"] > 0
