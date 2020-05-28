import os, sys, re
import numpy as np
from ase.io import read, Trajectory
from optparse import OptionParser
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from ase.visualize.plot import plot_atoms
from matplotlib import rc, rcParams
import matplotlib.font_manager
pgf_with_rc_fonts = {"pgf.texsystem": "pdflatex"}
matplotlib.rcParams.update(pgf_with_rc_fonts)
#rc('text', usetex = True)
rcParams['mathtext.fontset'] = 'cm'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.rm'] = 'Arial'

SMALL_SIZE = 14
MEDIUM_SIZE = 14.5
BIGGER_SIZE = 16
plt.rc('font', size=MEDIUM_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE) # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)   # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = ['Tex Gyre Heros'] + plt.rcParams['font.sans-serif']
#plt.rcParams['mathtext.fontset'] = 'dejavuserif'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-r", "--run", dest="run", help="Run the test")
parser.add_option("-a", "--analyze", dest="analyze", help="Analyze the run")
(options, args) = parser.parse_args()


# Load Directories
xyzDir = os.path.join(os.getcwd(), "xyz")
inDir = os.path.join(os.getcwd(), "in")
rawDir = os.path.join(os.getcwd(), "raw")
hesDir = os.path.join(os.getcwd(), "Hessians")   
trajEXPDir = os.path.join(os.getcwd(), "Trajectories_EXP")
trajIDDir = os.path.join(os.getcwd(), "Trajectories_ID")
picDir = os.path.join(os.getcwd(), "Pictures")

names = [k.split(".")[0] for k in os.listdir(trajEXPDir)]





    
    #fig, axarr = plt.subplots(1,1)
    #axarr.set_title("{}".format(str(i)))
    #axarr.set_xlabel("Distance, $\mathrm{\AA}$")
    #axarr.set_ylabel("Distance, $\mathrm{\AA}$")
    #plot_atoms(traj[-1], axarr, rotation=('0x,0y,0z'), offset=(0,0))
    #fig.savefig(os.path.join(os.getcwd(), 'Pictures', "{}.png".format(str(i[:-5]))), dpi=200, transparent=True)
    #plt.close(fig)






fig= plt.figure(figsize=(4, 4))
#objects = ["{:0>4}".format(i) for i in sorted(names)]
#objects = range(len(names))
#y_pos = np.arange(len(objects)) * 3
y_pos = np.arange(len(names)) 
performance = [len(Trajectory(os.path.join(trajIDDir, name+".traj")))/len(Trajectory(os.path.join(trajEXPDir, name+".traj")))
   for name in sorted(names)]
#performance = [len(Trajectories[0])/float(len(i)) if len(i)<2000 else 1000 for i in Trajectories]
colors=[]
#for x in performance:
    #if x==max(performance):
        #colors.append('#45FF42')
    #elif x==min(performance):
        #colors.append('#CC3D33')
    #else:
        #colors.append('#7872A3')
        
for x in performance:
    if x==max(performance):
        colors.append('#7872A3')
    elif x==min(performance):
        colors.append('#7872A3')
    else:
        colors.append('#7872A3')
        
#for i in range(len(forces)):
    #performance.append(float(len(identity_forces))/len(forces[i]))
    #objects.append('Lindh_vdW_Exp_{}_{}'.format(str(mu_alpha[i][0]), str(mu_alpha[i][1])))

plt.bar(y_pos, performance, width=.75, color=colors, align='center', alpha=0.75)
#plt.xticks(y_pos, objects, rotation=0)

#tickmnames = [str(i) for i in sorted(names) if int(i)%10==0]
#tickpositions = [int(i) * 3  for i in tickmnames]

#plt.xticks(tickmnames, tickpositions)
plt.xticks([])
plt.xlim(-0.75, max(y_pos)+0.75)
plt.ylim(0, 6)
plt.hlines(1, -0.75, max(y_pos)+0.75, color='k', linestyle='--')
plt.hlines(np.average(performance),-1.5, max(y_pos)+1.5, color='#CC3D33', linestyle='--')
#plt.hlines(min(performance), min(y_pos)-2, max(y_pos)+2, color='k', linestyle='--')
#plt.hlines(1, min(y_pos), max(y_pos), color='k', linestyle='--')
plt.ylabel('Performance gain')
plt.xlabel('\nSystem size')
#plt.title('Comparison of Exponential\nPreconditioned Hessians\nfor different Cu bulk systems')
plt.title('Different Cu bulk systems')
#fig.autofmt_xdate()
plt.tight_layout()
plt.savefig(os.path.join(picDir, "Exp_ID.png"), dpi=300)
#plt.show()








#hessians = {
    #"python Exp_Hessian.py"   : "-g geometry.in --A 3",
    #"python2.7 Lindh_Hessian.py" : "geometry.in --full Hessian_Lindh.dat",
    #"python vdW_Hessian.py"   : "-g geometry.in -m 1",
    #}

#hesDir = os.path.join(os.getcwd(), "Hessians")

#for i in hessians:
    #os.system("cd {} && {} {}".format(hesDir, i, hessians[i]))
    #print("Done {}".format(i))
    



#run = [
       #"Identity_run.py", 
       #"precon_Exp_run.py",
       #"precon_Lindh_run.py",
       #"precon_vdW_run.py" 
       #]

#for i in run:
    #os.system("python {}".format(i))
    #print("Done {}".format(i))
