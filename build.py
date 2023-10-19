#!/usr/bin/env python

from simtk.openmm import app
from simtk import unit
import simtk.openmm as omm

def build_by_seq(seq, forcefield):
    step = 0.22  # distance between consecutive backbone beads

    topo = omm.app.Topology()
    chain = topo.addChain('X')
    atoms = []
    positions = []
    idx_offset = 0

    for i, resSymbol in enumerate(seq):
        res = topo.addResidue(resSymbol, chain)
        atoms.append(topo.addAtom('B', None, res))  # 'B' represents backbone bead
        z_position = i * step  # Update the z-position based on the step and sequence index
        positions.append([0, 0, z_position] * unit.nanometers)  # Convert position to nanometers unit
        if i > 0:  # Add bond between consecutive backbone beads
            topo.addBond(atoms[i-1], atoms[i])

    return topo, positions