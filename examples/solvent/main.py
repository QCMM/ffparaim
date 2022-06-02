#!/usr/bin/python3

from ffparaim import FFparAIM

pdb_file = 'IPH_solv.pdb'
smiles = 'c1ccc(cc1)O'


# Script for non-bonded parameters derivation.
def main():
    nb_params = FFparAIM(pdb_file,
                         smiles,
                         qm_charge=0,
                         ligand_selection=':1',
                         n_updates=4,
                         sampling_time=0.1,
                         total_qm_calculations=4,
                         method='PBE',
                         basis='def2-TZVP',
                         forcefield='openff_unconstrained-2.0.0.offxml')
    nb_params.run(pickle=True, charges=True, lj=False)


if __name__ == '__main__':
    main()
