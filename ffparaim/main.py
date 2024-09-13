#!/usr/bin/python3
import sys

from ffparaim import FFparAIM

smiles = sys.argv[1]
pdb_file = sys.argv[2]


# Script for non-bonded parameters derivation.
def main():

    nb_params = FFparAIM(qm_charge=0,
                         ligand_selection='resname LIG',
                         n_updates=5,
                         sampling_time=0.5,
                         total_qm_calculations=5,
                         method='PBE')
    molecule, lig_structure, system_structure, system = nb_params.prepare(smiles, pdb_file)
    nb_params.run(molecule,
                  system_structure,
                  system,
                  off=True,
                  charges=True,
                  lj=True)


if __name__ == '__main__':
    main()
