"""
Unit and regression test for the cinoas package.
"""

# Import package, test suite, and other packages as needed
import cinoas
import sys
import pytest
import numpy as np
import psi4



def test_cinoas_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "cinoas" in sys.modules

def test_find_active_space():
    wfn = psi4.core.Wavefunction.from_file('wfn')  
    # get natural orbital occupations
    noocs_occ, noocs_vir, _ = cinoas.block_diag_opd(wfn)

    options = {
        'print_level': 1  ,
        'threshold_occ': 0.98,
        'threshold_vir': 0.98, 
    }

    results = cinoas.find_active_space(noocs_occ, noocs_vir, options)
    
    ref_nact = 6
    ref_active = [2, 0, 0, 2, 0, 0, 2, 0]
    
    assert results['nact'] == ref_nact
    assert np.array_equal(ref_active, results['active'])