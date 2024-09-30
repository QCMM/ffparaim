"""
Unit and regression test for the ffparaim package. Testing mdtools.py.
"""

# Import package, test suite, and other packages as needed

import os
import shutil
import pytest
import filecmp
import ffparaim
import pandas as pd
import ffparaim.mdtools as mdt

from openmm import System
from openmm.app import Simulation
from openmm.unit.quantity import Quantity
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
from parmed.structure import Structure
from importlib_resources import files, as_file
from numpy.testing import assert_equal, assert_allclose


def test_separate_components(tmpdir):
    os.chdir(tmpdir)
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as infile:
        mdt.separate_components(infile, 'resname MOL')
    assert filecmp.cmp(os.path.join(tmpdir, 'env.pdb'),
                       files('ffparaim.data').joinpath('env.pdb')) == 0
    assert filecmp.cmp(os.path.join(tmpdir, 'lig.pdb'),
                       files('ffparaim.data').joinpath('lig.pdb')) == 0


def test_separate_components_invalid():
    pytest.raises(TypeError, mdt.separate_components)
    pytest.raises(TypeError, mdt.separate_components, ligand_selection='resname MOL')
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as infile:
        pytest.raises(TypeError, mdt.separate_components, pdb_file=infile)
        pytest.raises(AttributeError, mdt.separate_components, infile, 0)


def test_define_molecule(tmpdir):
    os.chdir(tmpdir)
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as infile:
        mdt.separate_components(infile, 'resname MOL')
    mol = mdt.define_molecule('c1ccc(cc1)O')
    assert isinstance(mol, Molecule) is True


def test_define_molecule_invalid():
    pytest.raises(TypeError, mdt.define_molecule)
    pytest.raises(ValueError, mdt.define_molecule, 'CO', 'lig.pdb')
    with as_file(files('ffparaim.data').joinpath('lig.pdb')) as infile:
        pytest.raises(ValueError, mdt.define_molecule, 'CO', infile)


def test_define_forcefield():
    ff = mdt.define_forcefield('openff_unconstrained-2.0.0.offxml')
    assert isinstance(ff, ForceField)


def test_define_forcefield_invalid():
    pytest.raises(TypeError, mdt.define_forcefield)
    pytest.raises(OSError, mdt.define_forcefield, 'test')


def test_prepare_ligand(tmpdir):
    with as_file(files('ffparaim.data').joinpath('lig.pdb')) as lig_pdb_file:
        shutil.copy(lig_pdb_file, tmpdir)
    os.chdir(tmpdir)
    molecule = mdt.define_molecule('c1ccc(cc1)O')
    forcefield = mdt.define_forcefield('openff_unconstrained-2.0.0.offxml')
    lig_structure = mdt.prepare_ligand(molecule, forcefield)
    assert isinstance(lig_structure, Structure) is True


def test_prepare_ligand_invalid():
    pytest.raises(TypeError, mdt.prepare_ligand)
    molecule = mdt.define_molecule('c1ccc(cc1)O')
    pytest.raises(TypeError, mdt.prepare_ligand, molecule)
    forcefield = mdt.define_forcefield('openff_unconstrained-2.0.0.offxml')
    pytest.raises(FileNotFoundError,
                  mdt.prepare_ligand,
                  molecule,
                  forcefield,
                  lig_pdb_file='test.pdb')


def test_prepare_environment():
    ff_env = ['amber14-all.xml', 'amber14/tip3p.xml']
    with as_file(files('ffparaim.data').joinpath('env.pdb')) as env_pdb_file:
        env_structure = mdt.prepare_enviroment(ff_env, str(env_pdb_file))
    assert isinstance(env_structure, Structure) is True


def test_prepare_environment_invalid(tmpdir):
    pytest.raises(FileNotFoundError, mdt.prepare_enviroment)
    with as_file(files('ffparaim.data').joinpath('env.pdb')) as env_pdb_file:
        shutil.copy(env_pdb_file, tmpdir)
    os.chdir(tmpdir)
    pytest.raises(ValueError, mdt.prepare_enviroment, ff_env=[])


def test_create_system(tmpdir):
    with as_file(files('ffparaim.data').joinpath('lig.pdb')) as lig_pdb_file:
        shutil.copy(lig_pdb_file, tmpdir)
    os.chdir(tmpdir)
    molecule = mdt.define_molecule('c1ccc(cc1)O')
    forcefield = mdt.define_forcefield('openff_unconstrained-2.0.0.offxml')
    lig_structure = mdt.prepare_ligand(molecule, forcefield)
    with as_file(files('ffparaim.data').joinpath('env.pdb')) as env_pdb_file:
        env_structure = mdt.prepare_enviroment(env_pdb_file=str(env_pdb_file))
    system_structure, system = mdt.create_system(lig_structure, env_structure)
    assert isinstance(system_structure, Structure) is True
    assert isinstance(system, System) is True


def test_create_system_invalid(tmpdir):
    pytest.raises(TypeError, mdt.create_system)
    with as_file(files('ffparaim.data').joinpath('lig.pdb')) as lig_pdb_file:
        shutil.copy(lig_pdb_file, tmpdir)
    os.chdir(tmpdir)
    molecule = mdt.define_molecule('c1ccc(cc1)O')
    forcefield = mdt.define_forcefield('openff_unconstrained-2.0.0.offxml')
    lig_structure = mdt.prepare_ligand(molecule, forcefield)
    pytest.raises(TypeError, mdt.create_system, lig_structure=lig_structure)
    with as_file(files('ffparaim.data').joinpath('env.pdb')) as env_pdb_file:
        env_structure = mdt.prepare_enviroment(env_pdb_file=str(env_pdb_file))
    pytest.raises(TypeError, mdt.create_system, env_structure=env_structure)
    pytest.raises(TypeError, mdt.create_system, 'lig_structure', env_structure)
    pytest.raises(RecursionError, mdt.create_system, lig_structure, 'env_structure')
    pytest.raises(AttributeError, mdt.create_system, 'lig_structure', 'env_structure')


def test_save_serialized_system(tmpdir):
    with as_file(files('ffparaim.data').joinpath('lig.pdb')) as lig_pdb_file:
        shutil.copy(lig_pdb_file, tmpdir)
    os.chdir(tmpdir)
    molecule = mdt.define_molecule('c1ccc(cc1)O')
    forcefield = mdt.define_forcefield('openff_unconstrained-2.0.0.offxml')
    lig_structure = mdt.prepare_ligand(molecule, forcefield)
    with as_file(files('ffparaim.data').joinpath('env.pdb')) as env_pdb_file:
        env_structure = mdt.prepare_enviroment(env_pdb_file=str(env_pdb_file))
    system_structure, system = mdt.create_system(lig_structure, env_structure)
    mdt.save_serialized_system(system, 'system.xml')
    assert filecmp.cmp(os.path.join(tmpdir, 'system.xml'),
                       files('ffparaim.data').joinpath('system.xml')) == 0


def test_save_serialized_system_invalid(tmpdir):
    pytest.raises(TypeError, mdt.save_serialized_system)
    with as_file(files('ffparaim.data').joinpath('lig.pdb')) as lig_pdb_file:
        shutil.copy(lig_pdb_file, tmpdir)
    os.chdir(tmpdir)
    molecule = mdt.define_molecule('c1ccc(cc1)O')
    forcefield = mdt.define_forcefield('openff_unconstrained-2.0.0.offxml')
    lig_structure = mdt.prepare_ligand(molecule, forcefield)
    with as_file(files('ffparaim.data').joinpath('env.pdb')) as env_pdb_file:
        env_structure = mdt.prepare_enviroment(env_pdb_file=str(env_pdb_file))
    system_structure, system = mdt.create_system(lig_structure, env_structure)
    pytest.raises(TypeError, mdt.save_serialized_system, system)
    pytest.raises(TypeError, mdt.save_serialized_system, xml_file='system.xml')
    pytest.raises(ValueError, mdt.save_serialized_system, 'system', 'system.xml')


def test_prepare_off_lj():
    molecule = Molecule.from_smiles('CO')
    off_ff = mdt.define_forcefield('openff_unconstrained-2.0.0.offxml')
    sig = list(range(1, 7))
    eps = list(range(7, 13))
    smirks_dict = mdt.prepare_off_lj(molecule, off_ff, sig, eps)
    assert isinstance(smirks_dict, dict) is True
    assert list(smirks_dict.keys())[2] == '[#1:1]-[#6X4]-[#7,#8,#9,#16,#17,#35]'
    assert_equal(smirks_dict['[#6X4:1]']['sigma'][0], 1)
    assert_equal(smirks_dict['[#1:1]-[#8]']['epsilon'][0], 12)


def test_prepare_off_lj_invalid():
    pytest.raises(TypeError, mdt.prepare_off_lj)
    molecule = Molecule.from_smiles('CO')
    pytest.raises(TypeError, mdt.prepare_off_lj, molecule)
    off_ff = mdt.define_forcefield('openff_unconstrained-2.0.0.offxml')
    pytest.raises(TypeError, mdt.prepare_off_lj, molecule, off_ff)
    sig = list(range(1, 7))
    eps = list(range(7, 13))
    pytest.raises(TypeError, mdt.prepare_off_lj, molecule, off_ff, sig=sig)
    pytest.raises(TypeError, mdt.prepare_off_lj, molecule, off_ff, eps=eps)
    pytest.raises(AttributeError, mdt.prepare_off_lj, 'molecule', off_ff, sig, eps)
    pytest.raises(AttributeError, mdt.prepare_off_lj, molecule, 'off_ff', sig, eps)
    pytest.raises(IndexError, mdt.prepare_off_lj, molecule, off_ff, 'sig', eps)
    pytest.raises(TypeError, mdt.prepare_off_lj, molecule, off_ff, 0, 0)


def test_save_forcefield(tmpdir):
    with as_file(files('ffparaim.data').joinpath('lig.pdb')) as lig_pdb_file:
        shutil.copy(lig_pdb_file, tmpdir)
    os.chdir(tmpdir)
    molecule = mdt.define_molecule('c1ccc(cc1)O')
    forcefield = 'openff_unconstrained-2.0.0.offxml'
    off_ff = mdt.define_forcefield(forcefield)
    outfile = f'd-mbis_{forcefield}'
    with as_file(files('ffparaim.data').joinpath('ffparaim.dat')) as data:
        df = pd.read_csv(data)
    sig = df.sigma_nanometer.to_list()
    eps = df.epsilon_kjmol.to_list()
    smirks_dict = mdt.prepare_off_lj(molecule, off_ff, sig, eps)
    mdt.save_forcefield(off_ff, smirks_dict, outfile)
    assert filecmp.cmp(os.path.join(tmpdir, outfile),
                       files('ffparaim.data').joinpath(outfile)) == 0


def test_save_forcefield_invalid(tmpdir):
    pytest.raises(TypeError, mdt.save_forcefield)
    forcefield = 'openff_unconstrained-2.0.0.offxml'
    off_ff = mdt.define_forcefield(forcefield)
    with as_file(files('ffparaim.data').joinpath('lig.pdb')) as lig_pdb_file:
        shutil.copy(lig_pdb_file, tmpdir)
    os.chdir(tmpdir)
    molecule = mdt.define_molecule('c1ccc(cc1)O')
    with as_file(files('ffparaim.data').joinpath('ffparaim.dat')) as data:
        df = pd.read_csv(data)
    sig = df.sigma_nanometer.to_list()
    eps = df.epsilon_kjmol.to_list()
    smirks_dict = mdt.prepare_off_lj(molecule, off_ff, sig, eps)
    pytest.raises(TypeError, mdt.save_forcefield, 'off_ff', smirks_dict, 'test')
    pytest.raises(AttributeError, mdt.save_forcefield, off_ff, 'smirks_dict', 'test')
    pytest.raises(KeyError, mdt.save_forcefield, off_ff, {}, 'test')
    pytest.raises(TypeError, mdt.save_forcefield, off_ff, smirks_dict, 0)


def test_setup_simulation(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM(ligand_selection='resname MOL')
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as infile:
        molecule, lig_structure, system_structure, system = ffp.prepare('c1ccc(cc1)O', str(infile))
    positions = system_structure.positions
    system, simulation = mdt.setup_simulation(system_structure,
                                              system,
                                              positions,
                                              update=0,
                                              frames=500000,
                                              restraint_dict=None)
    assert isinstance(system, System) is True
    assert isinstance(simulation, Simulation) is True


def test_setup_simulation_invalid(tmpdir):
    pytest.raises(TypeError, mdt.setup_simulation)
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM(ligand_selection='resname MOL')
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as infile:
        molecule, lig_structure, system_structure, system = ffp.prepare('c1ccc(cc1)O', str(infile))
    pytest.raises(TypeError, mdt.setup_simulation, system_structure)
    pytest.raises(TypeError, mdt.setup_simulation, system_structure, system)
    positions = system_structure.positions
    pytest.raises(TypeError,
                  mdt.setup_simulation,
                  system_structure,
                  system,
                  positions)
    pytest.raises(TypeError,
                  mdt.setup_simulation,
                  system_structure,
                  system,
                  positions,
                  update=0)


def test_get_positions(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM(ligand_selection='resname MOL')
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as infile:
        molecule, lig_structure, system_structure, system = ffp.prepare('c1ccc(cc1)O', str(infile))
        system, simulation = mdt.setup_simulation(system_structure,
                                                  system,
                                                  positions=system_structure.positions,
                                                  update=0,
                                                  frames=500000,
                                                  restraint_dict=None)
    positions = mdt.get_positions(simulation)
    assert isinstance(positions, Quantity) is True
    assert_equal(positions.__len__(), 2551)
    assert_allclose(positions[-1][0]._value, 2.17278242, atol=0.05)
    assert positions[-1][0].unit._name == 'nanometer'
    assert filecmp.cmp(os.path.join(tmpdir, 'output.pdb'),
                       files('ffparaim.data').joinpath('output.pdb')) == 0


def test_get_positions_invalid():
    pytest.raises(TypeError, mdt.get_positions)
    pytest.raises(AttributeError, mdt.get_positions, 'simulation')


def test_image_molecule(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM(ligand_selection='resname MOL')
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as infile:
        molecule, lig_structure, system_structure, system = ffp.prepare('c1ccc(cc1)O', str(infile))
    system, simulation = mdt.setup_simulation(system_structure,
                                              system,
                                              positions=system_structure.positions,
                                              update=0,
                                              frames=500000,
                                              restraint_dict=None)
    mdt.get_positions(simulation)
    mdt.image_molecule()
    assert filecmp.cmp(os.path.join(tmpdir, 'output_recenter.pdb'),
                       files('ffparaim.data').joinpath('output_recenter.pdb')) == 0


def test_image_molecule_invalid():
    pytest.raises(OSError, mdt.image_molecule, 'test')


def test_get_atoms_idx(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM(ligand_selection='resname MOL')
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as infile:
        molecule, lig_structure, system_structure, system = ffp.prepare('c1ccc(cc1)O', str(infile))
    ligand_atom_list = mdt.get_atoms_idx(system_structure, 'resname MOL')
    assert isinstance(ligand_atom_list, list) is True
    assert_equal(ligand_atom_list[-1], 12)


def test_get_atoms_idx_invalid(tmpdir):
    os.chdir(tmpdir)
    pytest.raises(TypeError, mdt.get_atoms_idx)
    ffp = ffparaim.FFparAIM(ligand_selection='resname MOL')
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as infile:
        molecule, lig_structure, system_structure, system = ffp.prepare('c1ccc(cc1)O', str(infile))
    pytest.raises(TypeError, mdt.get_atoms_idx, system_structure)
    pytest.raises(AttributeError, mdt.get_atoms_idx, 'system_structure', 'resname MOL')
    pytest.raises(AttributeError, mdt.get_atoms_idx, system_structure, 0)


def test_update_params(tmpdir):
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM(ligand_selection='resname MOL')
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as infile:
        molecule, lig_structure, system_structure, system = ffp.prepare('c1ccc(cc1)O', str(infile))
    ligand_atom_list = mdt.get_atoms_idx(system_structure, 'resname MOL')
    with as_file(files('ffparaim.data').joinpath('ffparaim.dat')) as data:
        df = pd.read_csv(data)
    charge = df.atomic_charges.to_list()
    sigma = df.sigma_nanometer.to_list()
    epsilon = df.epsilon_kjmol.to_list()
    system = mdt.update_params(system,
                               ligand_atom_list,
                               charge,
                               sigma,
                               epsilon)
    assert isinstance(system, System) is True
    assert_allclose(system.getForce(3).getParticleParameters(0)[0]._value, charge[0])
    assert_allclose(system.getForce(3).getParticleParameters(0)[1]._value, sigma[0])
    assert_allclose(system.getForce(3).getParticleParameters(0)[2]._value, epsilon[0])


def test_update_params_invalid(tmpdir):
    pytest.raises(TypeError, mdt.update_params)
    os.chdir(tmpdir)
    ffp = ffparaim.FFparAIM(ligand_selection='resname MOL')
    with as_file(files('ffparaim.data').joinpath('solvent.pdb')) as infile:
        molecule, lig_structure, system_structure, system = ffp.prepare('c1ccc(cc1)O', str(infile))
    pytest.raises(TypeError, mdt.update_params, system)
    ligand_atom_list = mdt.get_atoms_idx(system_structure, 'resname MOL')
    pytest.raises(TypeError,
                  mdt.update_params,
                  system,
                  ligand_atom_list,
                  0)
    pytest.raises(IndexError,
                  mdt.update_params,
                  system,
                  ligand_atom_list,
                  [])
