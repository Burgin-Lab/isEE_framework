"""
Unit and regression test for the isee_framework package.
"""

# Import package, test suite, and other packages as needed
import isee_framework
import pytest
import sys

def test_isee_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "isee_framework" in sys.modules
