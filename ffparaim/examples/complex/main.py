#!/usr/bin/python3

from ffparaim import FFparAIM

pdb_file = 'holo_complex_solv.pdb'
smiles = 'c1ccc(cc1)O'
restraint = {'boresch': {'receptor': 'protein', 'ligand': 'resname IPH'}}


# Script for non-bonded parameters derivation.
def main():
    nb_params = FFparAIM(pdb_file,
                         smiles,
                         qm_charge=0,
                         ligand_selection=':IPH',
                         n_updates=3,
                         sampling_time=1,
                         total_qm_calculations=3,
                         method='PBE',
                         basis='def2-TZVP',
                         forcefield='openff_unconstrained-2.0.0.offxml')
    nb_params.run(restraint_dict=restraint,
                  pickle=True,
                  charges=True,
                  lj=True)


if __name__ == '__main__':
    main()
