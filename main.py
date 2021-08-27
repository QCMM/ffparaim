#!/usr/bin/python3
import sys

from ffparaim import FFparAIM

pdb_file = sys.argv[1]
smiles = sys.argv[2]


def main():
    charges_param = FFparAIM(pdb_file, smiles, qm_charge=0, ligand_selection=':1',
                             n_charge_updates=5, sampling_time=0.5, total_qm_calculations=5,
                             method='PBE', forcefield='openff_unconstrained-1.3.0.offxml')
    charges_param.run(json=True, charges=True, lj=True)


if __name__ == '__main__':
    main()
