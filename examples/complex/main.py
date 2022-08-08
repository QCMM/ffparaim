#!/usr/bin/python3

from ffparaim import FFparAIM

pdb_file = 'complex.pdb'
smiles = 'c1ccc(cc1)O'
restraint = {'boresch': {'receptor': [1829, 1737, 1607], 'ligand': [2, 3, 5]}}


# Script for non-bonded parameters derivation.
def main():

    nb_params = FFparAIM(qm_charge=0,
                         ligand_selection=':MOL',
                         n_updates=3,
                         sampling_time=0.5,
                         total_qm_calculations=3,
                         method='PBE')
    molecule, system_structure, system = nb_params.prepare(smiles, pdb_file)
    nb_params.run(molecule,
                  system_structure,
                  system,
                  restraint_dict=restraint,
                  off=True,
                  charges=True,
                  lj=True)


if __name__ == '__main__':
    main()
