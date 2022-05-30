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
                         n_updates=1,
                         sampling_time=0.1,
                         total_qm_calculations=1,
                         method='PBE',
                         forcefield='openff_unconstrained-1.3.0.offxml')
    nb_params.run(pickle=True, charges=True, lj=True)


if __name__ == '__main__':
    main()
