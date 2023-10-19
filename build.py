#!/usr/bin/env python

from simtk.openmm import app
from simtk import unit
import simtk.openmm as omm

def build_by_seq(seq, forcefield):
    step = 0.22

    # Define a simple geometry for the backbone bead of each amino acid
    geometry = {'A': [0, 0, 0], 'R': [0, 0, 0], 'N': [0, 0, 0], 'D': [0, 0, 0], 'C': [0, 0, 0],
                'E': [0, 0, 0], 'Q': [0, 0, 0], 'G': [0, 0, 0], 'H': [0, 0, 0], 'I': [0, 0, 0],
                'L': [0, 0, 0], 'K': [0, 0, 0], 'M': [0, 0, 0], 'F': [0, 0, 0], 'P': [0, 0, 0],
                'S': [0, 0, 0], 'T': [0, 0, 0], 'W': [0, 0, 0], 'Y': [0, 0, 0], 'V': [0, 0, 0]}

    topo = omm.app.Topology()
    chain = topo.addChain('X')
    atoms = []
    positions = []
    idx_offset = 0

    for i, resSymbol in enumerate(seq):
        res = topo.addResidue(resSymbol, chain)
        atoms.append(topo.addAtom('B', None, res))  # 'B' represents backbone bead
        positions.append(geometry[resSymbol])
        if i > 0:  # Add bond between consecutive backbone beads
            topo.addBond(atoms[i-1], atoms[i])

    return topo, positions
