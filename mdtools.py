#!/usr/bin/python3
import warnings
import mdtraj
import sys

import parmed as pmd

from openmmtools import forces
from simtk.openmm import openmm
from simtk.openmm import app
from simtk import unit
from openforcefield.topology import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openmmforcefields.generators import GAFFTemplateGenerator


# Avoid warnings.
warnings.filterwarnings('ignore')


def get_params(smiles, forcefield):
    """Get forcefield parameters for molecule.

    Parameters
    ----------
    smiles :    str
        SMILES string for the molecule.
    """

    # Create an openforcefield Molecule object.
    molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    if forcefield.startswith('openff'):
        # Create SMIRNOFF template generator.
        template = SMIRNOFFTemplateGenerator(molecules=molecule, forcefield=forcefield)
    elif forcefield.startswith('gaff'):
        # Create GAFF template generator.
        template = GAFFTemplateGenerator(molecules=molecule, forcefield=forcefield)
    else:
        raise Exception('Invalid forcefield.')
    return template, molecule


def create_forcefield(template):
    """Create an OpenMM ForceField object from an SMIRNOFF/GAFF template generator.

    Parameters
    ----------
    template :  openmmforcefields.generators.SMIRNOFFTemplateGenerator,
                openmmforcefields.generators.GAFFFFTemplateGenerator
        OpenMM SMIRNOFF or GAFF template.
    """

    # Create an OpenMM ForceField object using AMBER14 TIP3P with compatible ions.
    ff = app.ForceField('amber/tip3p_standard.xml')
    # Register the template generator.
    ff.registerTemplateGenerator(template.generator)
    return ff


def read_pdb(pdb_file):
    """Read a PDB file.

    Parameters
    ----------
    pdb_file :  str
        PDB filename.
    """

    # Load PDB file.
    pdb = app.PDBFile(pdb_file)
    return pdb


def create_system(ff, pdb, restraint=False, ligand_atoms_idx=None,
                  receptor_atoms_idx=None):
    """Create OpenMM system and save it in XML format.

    Parameters
    ----------
    ff :    simtk.openmm.app.ForceField
        OpenMM ForceField class object.
    pdb :   openmm.app.PDBFile
        OpenMM PDBFile class object.
    """

    # Create system.
    system = ff.createSystem(pdb.topology,
                             nonbondedMethod=app.PME,
                             nonbondedCutoff=1 * unit.nanometer,
                             constraints=app.HBonds)
    system.addForce(openmm.MonteCarloBarostat(1 * unit.bar, 298 * unit.kelvin))
    if restraint:
        if ligand_atoms_idx is not None and receptor_atoms_idx is not None:
            restraint = forces.HarmonicRestraintForce(
                spring_constant=0.2 * unit.kilocalories_per_mole / unit.angstrom**2,
                restrained_atom_indices1=ligand_atoms_idx,
                restrained_atom_indices2=receptor_atoms_idx)
            system.addForce(restraint)
        else:
            raise Exception("Missing atom list to apply restraints")
    return system


def serialize_system(system, xml_file):
    # Serialize system.
    system_serialized = openmm.XmlSerializer.serialize(system)
    # Save XML file.
    with open(xml_file, 'w') as f:
        f.write(system_serialized)


def get_atoms_idx(pdb, residue_selection):
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

    struct = pmd.openmm.load_topology(pdb.topology)
    atom_list = struct.atoms
    res = struct[residue_selection].atoms
    atoms_idx = [idx for at in res for idx, atom in enumerate(atom_list) if (
        at.residue.number, at.name) == (atom.residue.number, atom.name)]
    return atoms_idx


def setup_simulation(pdb, system, positions, update):
    """Setup the OpenMM simulation with the current topology and the input coordinates
     or the current positions depending on the value of update.
    Standard conditions are assumed (298K, 1bar).

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
    restraint : bool
    ligand_atom_list : list
    receptor_atom_list : list

    """

    integrator = openmm.LangevinIntegrator(
        298 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    simulation = app.Simulation(pdb.topology, system, integrator)
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
    traj.image_molecules()
    #traj.image_molecules(anchor_molecules=[traj.top.find_molecules()[0]], other_molecules=traj.top.find_molecules()[1:])
    pos = mdtraj.utils.in_units_of(traj.xyz[0], traj._distance_unit, "angstroms")
    ucl = mdtraj.utils.in_units_of(traj.unitcell_lengths[0], traj._distance_unit, "angstroms")
    pdb_recentered = mdtraj.formats.PDBTrajectoryFile("output_recenter.pdb", "w")
    pdb_recentered.write(pos, traj.top, modelIndex=None, unitcell_lengths=ucl,
                         unitcell_angles=traj.unitcell_angles[0])


def get_polar_hydrogens(molecule):
    # check polar hydrogens involved in O-H bond.
    polar_h = molecule.chemical_environment_matches('[#1:1]-[#8]')
    # create an atom index list for polar hydrogens.
    polar_h_idx = [h[0] for h in polar_h] if len(polar_h) > 0 else []
    return polar_h_idx


def update_params(system, ligand_atoms_idx, polar_h_idx, charge=None, sigma=None, epsilon=None):
    for i, idx in enumerate(ligand_atoms_idx):
        q = system.getForces()[2].getParticleParameters(
            idx)[0] if charge is None else round(charge[i], 3)
        sig = system.getForces()[2].getParticleParameters(
            idx)[1] if sigma is None or i in polar_h_idx else round(sigma[i], 3)
        eps = system.getForces()[2].getParticleParameters(
            idx)[2] if epsilon is None or i in polar_h_idx else round(epsilon[i], 3)
        system.getForces()[2].setParticleParameters(
            idx, charge=q, sigma=sig, epsilon=eps)
    return system
