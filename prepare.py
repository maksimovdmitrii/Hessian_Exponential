import os, sys, re, shutil
import numpy as np

from ase.io import read, write
from ase.calculators.lj import LennardJones
from ase.calculators.lammpslib import LAMMPSlib
from ase.build import bulk

# Create Directories
inDir = os.path.join(os.getcwd(), "in")
if not os.path.exists(inDir):
    os.mkdir(inDir)
    
hesDir = os.path.join(os.getcwd(), "Hessians")
if not os.path.exists(hesDir):
    os.mkdir(hesDir)
    
trajExpDir = os.path.join(os.getcwd(), "Trajectories_EXP")
if not os.path.exists(trajExpDir):
    os.mkdir(trajExpDir)

trajIDDir = os.path.join(os.getcwd(), "Trajectories_ID")
if not os.path.exists(trajIDDir):
    os.mkdir(trajIDDir)
    
picDir = os.path.join(os.getcwd(), "Pictures")
if not os.path.exists(picDir):
    os.mkdir(picDir)
    
# Produce FHI-aims files
for i in range(2, 10):   
    a0 = bulk('Cu', cubic=True)
    a0 *= [i, i, i]
    del a0[0]
    a0.rattle(0.05)
    write(os.path.join(inDir, "{:04d}.in".format(i)), a0, format="aims")
    
# Produce Hessians
for structure in os.listdir(inDir):
    inputfile = os.path.join(inDir, structure)
    outputfile = os.path.join(hesDir, structure.split(".")[0]+".hes")
    os.system("python Exp_Hessian.py -i {} -o {} -f aims --A 3".format(inputfile, outputfile))
    print("Hessian for {} is done".format(structure))
    



