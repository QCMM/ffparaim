#!/usr/bin/python3

from ffparaim import FFparAIM

pdb_file = 'holo_complex_solv.pdb'
smiles = 'c1ccc(cc1)O'
restraint = {'flatbottom': {'receptor': 'protein', 'ligand': 'resname IPH'}}


# Script for non-bonded parameters derivation.
def main():
    nb_params = FFparAIM(pdb_file,
                         smiles,
                         qm_charge=0,
                         ligand_selection=':163',
                         n_updates=1,
                         sampling_time=0.1,
                         total_qm_calculations=1,
                         method='B3LYP',
                         basis='def2-TZVP',
                         forcefield='openff_unconstrained-2.0.0.offxml')
    nb_params.run(restraint_dict=None,
                  pickle=True,
                  charges=True,
                  lj=True)


if __name__ == '__main__':
    main()
