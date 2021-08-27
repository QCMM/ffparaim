#!/usr/bin/python3

from iodata import load_one
from denspart.adapters.horton3 import prepare_input
from denspart.mbis import partition
from denspart.properties import compute_rcubed
from ffparaim import units
from ffparaim import utils
from rdkit.Chem.rdmolfiles import CanonicalRankAtoms

import numpy as np

# catch warnings
np.seterr(all='warn')


class ForceFieldDerivation(object):
    """docstring for ForceFieldDerivation."""

    def __init__(self):
        self.data = None  # data from the orca mkl file
        self.grid = None
        self.rho = None
        self.pro_model = None
        self.localgrids = None

    def load_data(self, infile):
        self.data = load_one(infile)  # data from the orca mkl file

    # set molgrid
    def set_molgrid(self, nrad, nang):
        self.grid, self.rho = prepare_input(self.data, nrad, nang, 10000)  # molgrid
        return

    # do partitioning depending on method
    def do_partitioning(self, method='mbis'):
        if method == 'hi':
            print('Pending')
            return
        elif method == 'mbis':
            # do MBIS partitioning
            print('MBIS partitioning ...')
            self.pro_model, self.localgrids = partition(
                self.data.atnums, self.data.atcoords, self.grid, self.rho)
        else:
            print('Invalid method')
            return
        return

    def get_charges(self):
        return self.pro_model.charges.tolist()

    def get_epol(self):
        orcalog = load_one('orca_pol_corr.out', fmt="orcalog")
        epol = (
            orcalog.extra["scf_energies"][0] - orcalog.extra["scf_energies"][-1]) / units.kcalmol
        return epol

    def get_rcubed(self):
        rc = compute_rcubed(self.pro_model, self.grid, self.rho, self.localgrids)
        return rc.tolist()


def symmetrize(molecule, params):
    # Create a RDKit Molecule object.
    mol = molecule.to_rdkit()
    # Get symmetry classes.
    symm_class = list(CanonicalRankAtoms(mol, breakTies=False))
    # Generate a symmetry dict.
    symm_dict = dict()
    for i, symm in enumerate(symm_class):
        if symm in symm_dict.keys():
            symm_dict[symm].append(params[i])
        else:
            symm_dict[symm] = [params[i]]
    # Create a symmetrized parameters list.
    symm_params = [np.array(symm_dict[symm]).mean() for symm in symm_class]
    return symm_params


def get_lj_params(molecule, rcubed, rcubed_table):
    sigma = []
    epsilon = []
    for i, atom in enumerate(molecule.atoms):
        vol_aim = rcubed[i]
        vol_isolated = rcubed_table[atom.atomic_number]
        scaling = vol_aim / vol_isolated
        alpha = utils.alpha_table[atom.atomic_number] * scaling
        c6 = utils.c6_table[atom.atomic_number] * scaling ** 2
        radius = 2.54 * alpha ** (1.0 / 7.0)
        rmin = 2 * radius
        sig = (rmin / (2 ** (1 / 6))) / units.nanometer
        eps = (c6 / (2 * rmin ** 6.0)) / units.kjmol
        sigma.append(sig)
        epsilon.append(eps)
    return sigma, epsilon
