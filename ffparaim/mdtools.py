#!/usr/bin/python3
import warnings
import mdtraj
import sys

import parmed as pmd
import numpy as np

from ffparaim.restraints import set_restraints
from openmm import openmm
from openmm import app
from openmm import unit
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from rdkit.Chem import AllChem, AddHs

# Avoid warnings.
warnings.filterwarnings('ignore')


def separate_components(pdb_file, ligand_selection):
    pdb = pmd.load_file(pdb_file)
    pmd.write_PDB(pdb[ligand_selection], 'lig.pdb')
    pmd.write_PDB(pdb[f'!{ligand_selection}'], 'env.pdb')
    return


def define_molecule(smiles, lig_pdb_file='lig.pdb'):
    template = AllChem.MolFromSmiles(smiles)
    template = AddHs(template)
    pdb = AllChem.MolFromPDBFile(lig_pdb_file, removeHs=False)
    rdmol = AllChem.AssignBondOrdersFromTemplate(template, pdb)
    return Molecule.from_rdkit(rdmol)
    # return Molecule.from_pdb_and_smiles(lig_pdb_file, smiles)


def define_forcefield(ff):
    return ForceField(ff)


def prepare_ligand(molecule, forcefield, lig_pdb_file='lig.pdb'):
    lig_pdb = pmd.load_file(lig_pdb_file)
    off_topology = Topology.from_openmm(openmm_topology=lig_pdb.topology,
                                        unique_molecules=[molecule])
    lig_system = forcefield.create_openmm_system(off_topology)
    lig_structure = pmd.openmm.load_topology(lig_pdb.topology,
                                             lig_system,
                                             xyz=lig_pdb.positions)
    return lig_structure


def prepare_enviroment(env_pdb_file='env.pdb'):
    omm_ff = app.ForceField('amber99sbildn.xml', 'tip3p.xml')
    env_pdb = app.PDBFile(env_pdb_file)
    omm_topology = env_pdb.getTopology()
    env_system = omm_ff.createSystem(omm_topology, rigidWater=False)
    env_structure = pmd.openmm.load_topology(omm_topology,
                                             env_system,
                                             xyz=env_pdb.getPositions())
    return env_structure


def create_system(lig_structure, env_structure):
    """Create OpenMM system and save it in XML format.

    Parameters
    ----------
    ff :    simtk.openmm.app.ForceField
        OpenMM ForceField class object.
    pdb :   openmm.app.PDBFile
        OpenMM PDBFile class object.
    """

    # Create system.
    system_structure = lig_structure + env_structure
    system = system_structure.createSystem(nonbondedMethod=app.PME,
                                           nonbondedCutoff=1 * unit.nanometer,
                                           constraints=app.HBonds)
    system.addForce(openmm.MonteCarloBarostat(1 * unit.bar, 298.15 * unit.kelvin))
    return system_structure, system


def save_serialized_system(system, xml_file):
    # Serialize system.
    system_serialized = openmm.XmlSerializer.serialize(system)
    # Save XML file.
    with open(xml_file, 'w') as f:
        f.write(system_serialized)


def create_smirks_dict(molecule, off_ff, sig, eps):
    # Get applied SMIRKS patterns
    applied_parameters = off_ff.label_molecules(molecule.to_topology())
    patterns = [value.smirks for value in applied_parameters[0]['vdW'].values()]
    # Create SMIRKS dict with vdW parameters.
    smirks_dict = dict()
    for index, pattern in enumerate(patterns):
        if pattern not in smirks_dict.keys():
            smirks_dict[pattern] = {'sigma': [sig[index]],
                                    'epsilon': [eps[index]]}
        else:
            smirks_dict[pattern]['sigma'].append(sig[index])
            smirks_dict[pattern]['epsilon'].append(eps[index])
    return smirks_dict


def save_forcefield(off_ff, smirks_dict, outfile):
    # Replace vdW parameters for every SMIRKS.
    for pattern in smirks_dict.keys():
        rmin_h = (np.array(smirks_dict[pattern]['sigma']).mean()) / (2 ** (5.0 / 6.0))
        eps = np.array(smirks_dict[pattern]['epsilon']).mean()
        # Define the new vdW parameters.
        vdw_handler = off_ff["vdW"].parameters[pattern]
        vdw_handler.rmin_half = round(rmin_h, 6) * unit.nanometer
        vdw_handler.epsilon = round(eps, 6) * unit.kilojoule_per_mole
    off_ff.to_file(outfile)


def setup_simulation(system_structure, system, positions, update, restraint_dict=None, ligand_atom_list=None):
    """Setup the OpenMM simulation with the current topology and the input coordinates
     or the current positions depending on the value of update.
    Standard conditions are assumed (298.15 K, 1bar).

    Parameters
    ----------
    pdb :  simtk.openmm.app.pdbfile.PDBFile
        OpenMM PDBFile object.
    system : simtk.openmm.openmm.System
        Openmm System object.
    positions :
        Current positions of atoms.
    update : int
        Value of charge update iteration.
    restraint_dict : dict

    """
    if restraint_dict is not None:
        set_restraints(system_structure.topology,
                       system,
                       positions,
                       restraint_dict,
                       ligand_atom_list)
    integrator = openmm.LangevinIntegrator(
        298.15 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    simulation = app.Simulation(system_structure.topology, system, integrator)
    simulation.reporters.append(app.StateDataReporter(
        sys.stdout, 5000, step=True, potentialEnergy=True, temperature=True, density=True))
    simulation.reporters.append(app.DCDReporter(f'traj_{update}.dcd', 50000))
    simulation.context.setPositions(positions)
    simulation.minimizeEnergy()
    return simulation


def get_positions(simulation):
    positions = simulation.context.getState(
        getPositions=True, enforcePeriodicBox=True).getPositions()
    _box_vectors = simulation.context.getState().getPeriodicBoxVectors()
    simulation.topology.setPeriodicBoxVectors(_box_vectors)
    app.PDBFile.writeFile(simulation.topology, positions,
                          open('output.pdb', 'w'))
    return positions


def image_molecule(pdbfile='output.pdb'):
    traj = mdtraj.load(pdbfile)
    try:
        traj.image_molecules()
    except ValueError:
        traj.image_molecules(anchor_molecules=[traj.top.find_molecules()[0]], other_molecules=traj.top.find_molecules()[1:])
    pos = mdtraj.utils.in_units_of(traj.xyz[0], traj._distance_unit, "angstroms")
    ucl = mdtraj.utils.in_units_of(traj.unitcell_lengths[0], traj._distance_unit, "angstroms")
    pdb_recentered = mdtraj.formats.PDBTrajectoryFile("output_recenter.pdb", "w")
    pdb_recentered.write(pos, traj.top, modelIndex=None, unitcell_lengths=ucl,
                         unitcell_angles=traj.unitcell_angles[0])


def get_atoms_idx(system_structure, residue_selection):
    """Generate an atom index lists for a residue selection.

    Parameters
    ----------
    system :    simtk.openmm.openmm.System
        Openmm System object.
    pdb :   openmm.app.PDBFile
        OpenMM PDBFile class object.
    residue_selection : str
        Amber-style mask string for residue atoms.
    """
    atom_list = system_structure.atoms
    res = system_structure[residue_selection].atoms
    atoms_idx = [idx for at in res for idx, atom in enumerate(atom_list) if (
        at.residue.number, at.name) == (atom.residue.number, atom.name)]
    return atoms_idx


def update_params(system, ligand_atoms_idx, charge=None, sigma=None, epsilon=None):
    for i, idx in enumerate(ligand_atoms_idx):
        q = system.getForces()[3].getParticleParameters(
            idx)[0] if charge is None else round(charge[i], 6)
        sig = system.getForces()[3].getParticleParameters(
            idx)[1] if sigma is None else round(sigma[i], 6)
        eps = system.getForces()[3].getParticleParameters(
            idx)[2] if epsilon is None else round(epsilon[i], 6)
        system.getForces()[3].setParticleParameters(
            idx, charge=q, sigma=sig, epsilon=eps)
    return system
