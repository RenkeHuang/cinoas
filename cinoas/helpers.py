"""
helpers.py
Contains functions to generate `wfn.npy` from a Psi4 CIS computation
"""

import subprocess
import json

def psi_cis_input(**kwargs):
    with open('example.json', 'r') as file:
        external = json.loads(file.read())  
    ethene_data = external['ethene']
    
    kwargs.setdefault('geom', ethene_data['geom'])
    kwargs.setdefault('num_roots', ethene_data['num_roots'])
    kwargs.setdefault('avg_states', ethene_data['avg_states'])

    kwargs.setdefault('basis', 'cc-pVDZ')

    return f"""
molecule {{
{kwargs['geom']}
}}

set {{
    basis                   {kwargs['basis']}
    scf_type                df  
    reference               rhf
    e_convergence           10   
    ex_level                1     # run CIS
    opdm                    true  # return one-particle density matrix
    num_roots               {kwargs['num_roots']}    
    avg_states              {kwargs['avg_states']}
}}

e, wfn = psi4.energy("detci", return_wfn=True)
wfn.to_file('wfn.npy')
    """

def create_psi_cis_wfn():
    with open('input.dat', 'w+') as file:
        file.write(psi_cis_input())
        print('Generate psi4 input.')
    
    subprocess.call(['psi4', '-v'])
    print('Generate wfn.npy')


if __name__ == "__main__":
    create_psi_cis_wfn()