#!/usr/bin/python3

from ffparaim import FFparAIM

pdb_file = 'IPH_solv.pdb'
smiles = 'c1ccc(cc1)O'

validation_dict = {'sampling_time': [1, 2, 3],
                   'n_updates': [1, 2, 3],
                   'total_qm_calculations': [1, 2, 3],
                   'method': ['BLYP', 'PBE'],
                   'basis': ['SVP', 'def2-TZVP']
                   }


# Script for non-bonded parameters derivation.
def main():
    nb_params = FFparAIM(pdb_file,
                         smiles,
                         qm_charge=0,
                         ligand_selection=':IPH',
                         n_updates=5,
                         sampling_time=25,
                         total_qm_calculations=100,
                         method='B3LYP',
                         basis='def2-TZVP',
                         forcefield='openff_unconstrained-2.0.0.offxml')
    nb_params.validation(sampling_time=validation_dict['sampling_time'],
                         n_updates=validation_dict['n_updates'],
                         total_qm_calculations=validation_dict['total_qm_calculations'],
                         method=validation_dict['method'],
                         basis=validation_dict['basis'])


if __name__ == '__main__':
    main()
