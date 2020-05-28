import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.optimize.precon import Exp, PreconLBFGS
from ase.io import read, write
from ase.optimize import BFGS
import os
import pandas as pd

from ase.calculators.loggingcalc import LoggingCalculator
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import optparse
parser = optparse.OptionParser()
parser.add_option("-f", "--formatfile", dest="formatfile", default='aims',
                    help="Input format extension ")
parser.add_option("-i", "--inputfile", dest="inputfile", help="input geometry file")
parser.add_option("--hessian", dest="hessian", help="Preconditionered Hessian file")
parser.add_option("-o", "--outputfile", dest="outputfile", help="output Trajectory file")
#parser.add_option("-m", "--multiplier", dest="multiplier", help="vdW multiplier for C6 coeff")
options, args = parser.parse_args()

# Create Directories
inDir = os.path.join(os.getcwd(), "in")  
hesDir = os.path.join(os.getcwd(), "Hessians")  
trajExpDir = os.path.join(os.getcwd(), "Trajectories_EXP")
trajIDDir = os.path.join(os.getcwd(), "Trajectories_ID")  
picDir = os.path.join(os.getcwd(), "Pictures")

#a0 = bulk('Cu', cubic=True)
#a0 *= [3, 3, 3]
#del a0[0]
#a0.rattle(0.1)
a0 = read(options.inputfile, format = options.formatfile)

nsteps = []
energies = []
log_calc = LoggingCalculator(EMT())

atoms = a0.copy()
atoms.set_calculator(log_calc)
opt = BFGS(atoms, trajectory=options.outputfile)
opt.H0 = pd.read_csv(options.hessian, sep='\s+', header=None)
opt.run(fmax=1e-3)
#atoms.write('relaxed.xyz', format = 'xyz')
#print('{} Done\n'.format(options.name))
#sys.exit(0)


#opt.run(fmax=1e-3)

#log_calc.plot(markers=['r-', 'b-'], energy=False, lw=2)
#plt.savefig("precon_exp.png")
#write('structure.in', a0, format = 'aims')
