"""
Unit and regression test for the ffparaim package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import ffparaim


def test_ffparaim_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "ffparaim" in sys.modules
