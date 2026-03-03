"""Tests for virtual gel plotting."""

import os
import tempfile

import pytest

try:
    from polymerase_tm.gel import plot_virtual_gel, _HAS_VIZ
except ImportError:
    _HAS_VIZ = False


@pytest.mark.skipif(not _HAS_VIZ, reason="Seaborn/Matplotlib not installed")
def test_plot_virtual_gel_1kb_plus():
    """Test generating a virtual gel with the 1kb Plus ladder."""
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = os.path.join(tmpdir, "test_gel_1kb.png")
        
        # Test drawing a 2500bp amplicon
        success = plot_virtual_gel(2500, ladder_name="1kb_plus", output_path=output_path)
        
        assert success is True
        assert os.path.exists(output_path)
        assert os.path.getsize(output_path) > 1000  # Should be a real image file


@pytest.mark.skipif(not _HAS_VIZ, reason="Seaborn/Matplotlib not installed")
def test_plot_virtual_gel_100bp():
    """Test generating a virtual gel with the 100bp ladder."""
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = os.path.join(tmpdir, "test_gel_100bp.png")
        
        # Test drawing a 300bp amplicon
        success = plot_virtual_gel(300, ladder_name="100bp", output_path=output_path)
        
        assert success is True
        assert os.path.exists(output_path)
        assert os.path.getsize(output_path) > 1000
