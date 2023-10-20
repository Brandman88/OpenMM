from simtk.openmm import app
import simtk.openmm as omm

def transform_geometry(atoms, transform):
    return [atoms[i] + transform[i] for i in range(0, 3)]

def build_by_seq(seq, forcefield):
    step = 0.22
    # Only backbone bead geometry
    geometry = {'ALA': {'B': [0, 0, 0]},
                'ARG': {'B': [0, 0, 0]},
                'ASN': {'B': [0, 0, 0]},
                'ASP': {'B': [0, 0, 0]},
                'CYS': {'B': [0, 0, 0]},
                'GLU': {'B': [0, 0, 0]},
                'GLN': {'B': [0, 0, 0]},
                'GLY': {'B': [0, 0, 0]},
                'HIS': {'B': [0, 0, 0]},
                'ILE': {'B': [0, 0, 0]},
                'LEU': {'B': [0, 0, 0]},
                'LYS': {'B': [0, 0, 0]},
                'MET': {'B': [0, 0, 0]},
                'PHE': {'B': [0, 0, 0]},
                'PRO': {'B': [0, 0, 0]},
                'SER': {'B': [0, 0, 0]},
                'THR': {'B': [0, 0, 0]},
                'TRP': {'B': [0, 0, 0]},
                'TYR': {'B': [0, 0, 0]},
                'VAL': {'B': [0, 0, 0]}}
    name_map = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
                'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
                'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
                'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}

    topo = omm.app.Topology()
    chain = topo.addChain('X')
    atoms = []
    positions = []
    idx_offset = 0
    transformStack = [[0, 0, 0]]
    last_idx = None

    for i, resSymbol in enumerate(seq):
        symbol = name_map[resSymbol]
        geometry_for_res = geometry[symbol]
        if i == 0 or i == len(seq) - 1:
            symbol = symbol + "T"

        res = topo.addResidue(symbol, chain)
        for atom in forcefield._templates[symbol].atoms:
            topo.addAtom(atom.name, forcefield._atomTypes[atom.type].element, res)
            if atom.name in geometry_for_res:
                if i % 2 == 0:
                    positions.append(geometry_for_res[atom.name])
                else:
                    new_geo = [-element for element in geometry_for_res[atom.name]]
                    positions.append(new_geo)
    transformed_positions = []
    currStackIdx = 0
    transformStack.append([0, 0, 0])
    curr_t = transformStack[0]
    for pos in positions:
        transformed_positions.append(transform_geometry(pos, curr_t))
        curr_t = [curr_t[i] + transformStack[1][i] for i in range(0, 3)]

    return topo, transformed_positions