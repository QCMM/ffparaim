#!/usr/bin/python3

import numpy as np

from iodata import load_one

from denspart.adapters.horton3 import prepare_input
from denspart.mbis import partition
from denspart.properties import compute_rcubed

from ffparaim import units
from ffparaim import utils

from rdkit.Chem.rdmolfiles import CanonicalRankAtoms

# Catch warnings.
np.seterr(all='warn')


class ForceFieldDerivation(object):
    """docstring for ForceFieldDerivation."""

    def __init__(self):
        self.data = None
        self.grid = None
        self.rho = None
        self.pro_model = None
        self.localgrids = None

    def load_data(self, infile):
        '''Use IOData.load_one function to parse and store
        information about coordinates for the electron density
        and molecular orbitals needed to atoms-in-molecules
        partitioning.'''

        # Load data from electron density input file.
        self.data = load_one(infile)

    def set_molgrid(self, nrad, nang):
        '''Define a molecular integration grid considering the
        number of radial and angular grid points. '''

        # Define the molecular grid.
        self.grid, self.rho = prepare_input(self.data,
                                            nrad,
                                            nang,
                                            10000)
        return

    def do_partitioning(self, method='mbis'):
        '''Apply a molecular density partitioning method to
        get a model of pro-molecular density and local grids
        for each basis function.'''

        # Apply Iterative-Hirshfeld method.
        if method == 'hi':
            print('Pending ...')
            return
        # Apply Minimal Basis Iterative Stockholder method.
        elif method == 'mbis':
            print('MBIS partitioning ...')
            # Get pro-molecular model and local grids.
            self.pro_model, self.localgrids = partition(
                self.data.atnums, self.data.atcoords, self.grid, self.rho)
        else:
            print('Invalid method')
            return
        return

    def get_charges(self):
        '''Partial atomic charges from molecular density partitioning.'''
        return self.pro_model.charges.tolist()

    def get_epol(self):
        '''Energy associated to polarize the electron density
        in a particular molecular enviroment.'''

        # Parse data from ORCA log file.
        orcalog = load_one('orca_pol_corr.out', fmt="orcalog")
        # Get energies from self-consistent-field calculation.
        scf_e = orcalog.extra['scf_energies']
        # Calculate polarization energy in kcal/mol.
        epol = (scf_e[0] - scf_e[-1]) / units.kcalmol
        return epol

    def get_rcubed(self):
        '''Third radial moment from molecular density partitioning.'''

        # Calculate third radial moment.
        rcubed = compute_rcubed(self.pro_model,
                                self.grid,
                                self.rho,
                                self.localgrids)
        return rcubed.tolist()


def symmetrize(molecule, params):
    '''Detect chemically equivalent atoms in a molecule and
    average parameters for every symmetric atom.'''

    # Create a RDKit Molecule object from Open Force Field Molecule class.
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
    '''Get Lennard-Jones parameters sigma and epsilon.
    These parameters describe the distance at which
    particle-particle potential energy is zero and
    the depth of the potential well (dispersion energy).'''

    # Create lists for store sigma and epsilon values.
    sigma, epsilon = [], []
    # Iterate for every atom in the molecule.
    for i, atom in enumerate(molecule.atoms):
        # Effective volume for atom in pro-molecule model.
        vol_aim = rcubed[i]
        # Effective volume for isolated atom.
        vol_isolated = rcubed_table[atom.atomic_number]
        # Effective volume scaling factor.
        scaling = vol_aim / vol_isolated
        # Static polarizability.
        alpha = utils.alpha_table[atom.atomic_number] * scaling
        # C6 dispersion coefficient proposed by Tkatchenko.
        c6 = utils.c6_table[atom.atomic_number] * scaling ** 2
        # Van-der-Waals radius using the Fedorov-Tkatchenko relation.
        radius = 2.54 * alpha ** (1.0 / 7.0)
        # Lennard-Jones radius.
        rmin = 2 * radius
        # Get sigma value in nanometers.
        sig = (rmin / (2 ** (1.0 / 6.0))) / units.nanometer
        # Get epsilon value in kJ/mol (Temporary for OpenMM XML system).
        eps = (c6 / (2 * rmin ** 6.0)) / units.kjmol
        # Add values to lists.
        sigma.append(sig)
        epsilon.append(eps)
    return sigma, epsilon
