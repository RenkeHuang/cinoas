"""
cinoas.py
Contains functions:
 - block_diag_opd
 - find_active_space
"""
import psi4
import numpy as np


def block_diag_opd(wfn):
    """ Block-diagonalize the one-particle reduced density matrix (1-RDM) of state-averaged configuration interaction singles.
        
    Parameters
    ----------
    wfn: instance of psi4.core.Wavefunction class
    
    Returns
    -------
    noocs_occ: instance of psi4.core.Vector class
        natural orbital occupations of the occupied block
    noocs_vir: instance of psi4.core.Vector class
        natural orbital occupations of the virtual block
    CaN: instance of psi4.core.Matrix
        coefficient matrix that transforms atomic-orbital basis to natural-orbital basis
    """
    
    docc = wfn.doccpi()
    da_mo = wfn.Da_subset('MO').to_array()  ## `numpy.ndarray` in C1, 'tuple' otherwise
    
    if wfn.nirrep() == 1:
        da_mo = [da_mo, ]
        print('Run in C1 symmetry')
    
    ### Get o-o/v-v blocks of 1-RDM in each irrep
    da_occ = [da_mo[h][:occ, :occ] for h, occ in enumerate(docc)]
    da_vir = [da_mo[h][occ:, occ:] for h, occ in enumerate(docc)]          

    da_occ = psi4.core.Matrix.from_array(da_occ)
    da_vir = psi4.core.Matrix.from_array(da_vir)

    ### Diagonalize Da_occ, Da_vir separately, U^\dagger*D*U = \Omega
    nocc = da_occ.coldim()   # same as docc
    u_o = psi4.core.Matrix("U_occ", nocc, nocc)
    noocs_occ = psi4.core.Vector("eval_occ", nocc)
    da_occ.diagonalize(u_o, noocs_occ, psi4.core.DiagonalizeOrder.Descending)

    nvir = da_vir.coldim()
    u_v = psi4.core.Matrix("U_vir", nvir, nvir)
    noocs_vir = psi4.core.Vector("eval_vir", nvir)
    da_vir.diagonalize(u_v, noocs_vir, psi4.core.DiagonalizeOrder.Descending)    

    ### Transform Ca from MO basis to CIS-NO basis
    if wfn.nirrep() == 1:
        ## C1 symmetry
        occ = docc.to_tuple()[0]
        Ca, CaNO = wfn.Ca().to_array(), wfn.Ca().to_array() 
        # occupied block 
        CaNO[:, :occ] = psi4.core.Matrix.doublet(psi4.core.Matrix.from_array(Ca[:,:occ]), u_o, False, False)  # core.doublet(A, B, transA, transB)
        # virtual block
        CaNO[:, occ:] = psi4.core.Matrix.doublet(psi4.core.Matrix.from_array(Ca[:, occ:]), u_v, False, False)
    
    else:
        CaNO = []
        # loop over each irrep
        for h, occ in enumerate(docc):
            Ca_h = wfn.Ca().to_array()[h]
            u_o_h = psi4.core.Matrix.from_array(u_o.to_array()[h])
            u_v_h = psi4.core.Matrix.from_array(u_v.to_array()[h])

            ### transform occupied block of irrep_h
            CaNO_h_occ = psi4.core.Matrix.doublet(psi4.core.Matrix.from_array(Ca_h[:,:occ]), u_o_h, False, False)
            Ca_h[:,:occ] = CaNO_h_occ
            ### transform virtual block of irrep_h
            CaNO_h_vir= psi4.core.Matrix.doublet(psi4.core.Matrix.from_array(Ca_h[:, occ:]), u_v_h, False, False)
            Ca_h[:, occ:] = CaNO_h_vir
            CaNO.append(Ca_h)

    CaNO = psi4.core.Matrix.from_array(CaNO)
    return noocs_occ, noocs_vir, CaNO



def find_active_space(noocs_occ, noocs_vir, options):
    """ Block-diagonalize the one-particle reduced density matrix (1-RDM) of state-averaged configuration interaction singles.
        
    Parameters
    ----------
    noocs_occ: instance of psi4.core.Vector class
        natural orbital occupations of the occupied block
    
    noocs_vir: instance of psi4.core.Vector class
        natural orbital occupations of the virtual block
    
    options: dict
        provide different metrices for selecting active space and printing preference.
        Available keys:
        - print_level: 1
        - threshold_occ:
        - threshold_vir:
        - num_act_occ:
        - num_act_vir

    Returns
    -------
    results: dict
        Available keys:
        - active: list, number of active orbitals per irep
        - nact: int
        - nact_occ: int
        - nact_vir: int        
        - sigma_o: float
        - sigma_v: float
        - occ_denominator: float
        - vir_denominator: float
    """
    
    print_level = options['print_level'] if 'print_level' in options.keys() else 0
    results = {}
    active = [0 for _ in range(noocs_occ.nirrep())]
    
    noocs_occupied = []
    noocs_virtual = []
    occ_denominator = 0
    vir_denominator = 0
    
    if noocs_occ.nirrep() == 1:
        # for C1 symmetry
        noocs_occupied += [(0, occus, index_in_irrep) for index_in_irrep, occus in enumerate(noocs_occ.to_array())]
        noocs_virtual  += [(0, occus, index_in_irrep) for index_in_irrep, occus in enumerate(noocs_vir.to_array())]

        occ_denominator += noocs_occ.to_array().size - sum(noocs_occ.to_array())
        vir_denominator += sum(noocs_vir.to_array())        
        if print_level > 1:
            print(f'occ_denominator = {occ_denominator}\nvir_denominator = {vir_denominator}')  
    
    else: 
        sym = 0
        for noocs_occupied_h, noocs_virtual_h in zip(noocs_occ.to_array(), noocs_vir.to_array()):
            n_occ_h = noocs_occupied_h.size
            noocs_occupied += [(sym, occus, index_in_irrep) for index_in_irrep, occus in enumerate(noocs_occupied_h)]
            noocs_virtual += [(sym, occus, index_in_irrep) for index_in_irrep, occus in enumerate(noocs_virtual_h)]

            occ_denominator += n_occ_h - sum(noocs_occupied_h)
            vir_denominator += sum(noocs_virtual_h)
            sym += 1
            if print_level > 1:
                print(f'occ_denominator = {occ_denominator}\nvir_denominator = {vir_denominator}')  
    
    noocs_occupied.sort(key=lambda h_occus_tuple: h_occus_tuple[1])
    noocs_virtual.sort(reverse=True, key=lambda h_occus_tuple: h_occus_tuple[1])
    
    occ_numerator = 0
    vir_numerator = 0
    ## Select active space based on CINO occupations
    if 'threshold_occ' in options.keys():
        nact_occ = 0
        for irrep, occus, i in noocs_occupied:
            occ_numerator += 1-occus
            if (occ_numerator/occ_denominator > options['threshold_occ']):
                occ_numerator -= 1-occus
                break
            active[irrep] += 1
            nact_occ += 1
            if print_level > 0:
                print(f'add occupied orb {i} in irrep {irrep}, CIS-NO occupation: {occus}')
        results['sigma_o'] = occ_numerator/occ_denominator
        results['nact_occ'] = nact_occ

        
    if 'threshold_vir' in options.keys():
        nact_vir = 0
        for irrep, occus, i in noocs_virtual:
            vir_numerator += occus
            if (vir_numerator/vir_denominator > options['threshold_vir']):
                vir_numerator -= occus
                break
            active[irrep] += 1  
            nact_vir += 1
            if print_level > 0:
                print(f'add  virtual orb {i} in irrep {irrep}, CIS-NO occupation: {occus}')
        results['sigma_v'] = vir_numerator/vir_denominator
        results['nact_vir'] = nact_vir
    
    
    if 'num_act_occ' in options.keys():
        nact_occ = options['num_act_occ']
        for k in range(nact_occ):
            irrep, occus, i = noocs_occupied[k]
            occ_numerator += 1-occus
            active[irrep] += 1
            if print_level > 0:
                print(f'add occupied orb {i} in irrep {irrep}, CIS-NO occupation: {occus}')
        results['sigma_o'] = occ_numerator/occ_denominator
            
    if 'num_act_vir' in options.keys():
        nact_vir = options['num_act_vir']
        for k in range(nact_vir):
            irrep, occus, i = noocs_virtual[k]
            vir_numerator += occus
            active[irrep] += 1
            if print_level > 0:
                print(f'add  virtual orb {i} in irrep {irrep}, CIS-NO occupation: {occus}')
        results['sigma_v'] = vir_numerator/vir_denominator
        
    results['active'] = active
    results['nact'] = sum(active)
    results['occ_denominator'] = occ_denominator
    results['vir_denominator'] = vir_denominator
    
    
    if print_level > 0: 
        print(f'\n  ==> Summary <==')
        print(f'num. of active orbitals = {sum(active)}')
        print(f'active   {active}')
        print('--------------------')

    return results



def main():
    wfn = psi4.core.Wavefunction.from_file('wfn.npy')  
    # get natural orbital occupations
    noocs_occ, noocs_vir, CaNO = block_diag_opd(wfn)

    wfn.Ca().copy(CaNO)
    wfn.Cb().copy(CaNO)

    sigma = 0.98
    options = {
        'print_level': 1  ,
        'threshold_occ': sigma,
        'threshold_vir': sigma, 
        # 'num_act_occ': 2,
        # 'num_act_vir': 2,
    }

    results = find_active_space(noocs_occ, noocs_vir, options)
    print(f'Results:\n{results}')
    
    
if __name__ == "__main__":
    main()

