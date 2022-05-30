#!/usr/bin/python3

from ffparaim import FFparAIM

pdb_file = 'solvent.pdb'
smiles = 'CC(=O)OCc1ccccc1'


# Script for non-bonded parameters derivation.
def main():
    nb_params = FFparAIM(pdb_file,
                         smiles,
                         qm_charge=0,
                         ligand_selection=':1',
                         receptor_selection=None,
                         n_updates=5,
                         sampling_time=25,
                         total_qm_calculations=100,
                         method='B3LYP',
                         basis='def2-TZVP',
                         forcefield='openff_unconstrained-2.0.0.offxml')
    nb_params.run(pickle=True, charges=True, lj=True)


if __name__ == '__main__':
    main()
