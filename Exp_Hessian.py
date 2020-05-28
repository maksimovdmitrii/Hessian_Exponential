import numpy as np
from ase.io import read, write
from ase.constraints import Filter, FixAtoms
from ase.geometry.cell import cell_to_cellpar
from ase.optimize.precon.neighbors import estimate_nearest_neighbour_distance


import optparse
parser = optparse.OptionParser()
parser.add_option("-f", "--formatfile", dest="formatfile", default='aims',
                    help="Input format extension ")
parser.add_option("-i", "--inputfile", dest="inputfile", help="input geometry file")
parser.add_option("-o", "--outputfile", dest="outputfile", help="output Hessian file")
parser.add_option("--A", type=float, dest="A",
                    help="Parameter A")
options, args = parser.parse_args()


def vector_separation(cell_h, cell_ih, qi, qj):
    # This file is part of i-PI.
    # i-PI Copyright (C) 2014-2015 i-PI developers
    # See the "licenses" directory for full license information.


    """Calculates the vector separating two atoms.

       Note that minimum image convention is used, so only the image of
       atom j that is the shortest distance from atom i is considered.

       Also note that while this may not work if the simulation
       box is highly skewed from orthorhombic, as
       in this case it is possible to return a distance less than the
       nearest neighbour distance. However, this will not be of
       importance unless the cut-off radius is more than half the
       width of the shortest face-face distance of the simulation box,
       which should never be the case.

       Args:
          cell_h: The simulation box cell vector matrix.
          cell_ih: The inverse of the simulation box cell vector matrix.
          qi: The position vector of atom i.
          qj: The position vectors of one or many atoms j shaped as (N, 3).
       Returns:
          dij: The vectors separating atoms i and {j}.
          rij: The distances between atoms i and {j}.
    """

    sij = np.dot(cell_ih, (qi - qj).T)  # column vectors needed
    sij -= np.rint(sij)

    dij = np.dot(cell_h, sij).T         # back to i-pi shape
    rij = np.linalg.norm(dij, axis=1)

    return dij, rij


atoms = read(options.inputfile, format = options.formatfile)
N  = len(atoms)
A=options.A
r_NN = estimate_nearest_neighbour_distance(atoms)
r_cut = 2.0 * r_NN
coordinates = atoms.get_positions()
cell_h = atoms.get_cell()[:]
cell_ih = atoms.get_reciprocal_cell()[:]
mu = 1.0

hessian = np.zeros(shape = (3 * N, 3 * N))
for i in range(N-1):
    qi = coordinates[i].reshape(1,3)
    qj = coordinates[i+1:].reshape(-1,3)
    dij, rij = vector_separation(cell_h, cell_ih, qi, qj)
    
    coeff = -mu * np.exp(-A * (rij / r_NN - 1))
    mask = np.array(rij>=r_cut)
    coeff[mask] = 0

    stack = np.hstack([np.identity(3)*coef for coef in coeff])

    hessian[3 * i    , 3 *(i+1):] = stack[0]
    hessian[3 * i + 1, 3 *(i+1):] = stack[1]
    hessian[3 * i + 2, 3 *(i+1):] = stack[2]
    
hessian = hessian + hessian.T - np.diag(hessian.diagonal())
for ind in range(len(hessian)):
    hessian[ind, ind] = -np.sum(hessian[ind])

with open(options.outputfile, 'w') as hes:
    for i in hessian:
        hes.write(' '.join(["{:8.3f}".format(k) for k in i]))
        hes.write('\n')
