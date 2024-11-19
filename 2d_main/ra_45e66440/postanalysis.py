import numpy as np
import dedalus.public as d3
import logging
logger = logging.getLogger(__name__)
import copy
import h5py
import numpy as np
import matplotlib
import re

import matplotlib.pyplot as plt
from dedalus.extras import plot_tools

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize

import os
from os import listdir
# Parameters
Lx, Lz = 40,1
Nx, Nz = 2560, 64
Ra_M = 4.5e6
D_0 = 0
D_H = 1/3
M_0 = 0
M_H = -1
N_s2=4/3
f=0.05

Prandtl = 0.7
dealias = 3/2
stop_sim_time = 500
timestepper = d3.RK222
max_timestep = 0.125
dtype = np.float64

# Bases
coords = d3.CartesianCoordinates('x','z')
dist = d3.Distributor(coords, dtype=dtype)
xbasis = d3.RealFourier(coords['x'], size=Nx, bounds=(0, Lx), dealias=dealias)
zbasis = d3.ChebyshevT(coords['z'], size=Nz, bounds=(0, Lz), dealias=dealias)

# Fields
p = dist.Field(name='p', bases=(xbasis,zbasis))
D = dist.Field(name='D', bases=(xbasis,zbasis))
M = dist.Field(name='M', bases=(xbasis,zbasis))
u = dist.VectorField(coords, name='u', bases=(xbasis,zbasis))
uy = dist.Field(name='uy', bases=(xbasis,zbasis))
Z = dist.Field(name='Z', bases=zbasis)
tau_p = dist.Field(name='tau_p')
tau_B1 = dist.Field(name='tau_B1', bases=xbasis)
tau_B2 = dist.Field(name='tau_B2', bases=xbasis)
tau_D1 = dist.Field(name='tau_D1', bases=xbasis)
tau_D2 = dist.Field(name='tau_D2', bases=xbasis)
tau_M1 = dist.Field(name='tau_M1', bases=xbasis)
tau_M2 = dist.Field(name='tau_M2', bases=xbasis)
tau_u1 = dist.VectorField(coords,name='tau_u1', bases=xbasis)
tau_u2 = dist.VectorField(coords,name='tau_u2', bases=xbasis)
tau_u3 = dist.Field(name='tau_u3', bases=xbasis)
tau_u4 = dist.Field(name='tau_u4', bases=xbasis)

# Substitutions    
#Kuo_Bretherton Equilibrium
kappa = (Ra_M * Prandtl/((M_0-M_H)*Lz**3))**(-1/2)
nu = (Ra_M / (Prandtl*(M_0-M_H)*Lz**3))**(-1/2)
print('kappa',kappa)
print('nu',nu)
Td=Lz**2/(nu*kappa)**(1/2)
Tc=(Lz/(M_0-M_H))**(1/2)
Tr=1/f
R_0=Tr/Tc
print('R_0',R_0)


x,z = dist.local_grids(xbasis,zbasis)
Z['g']=z
Z.change_scales(3/2)

ex,ez = coords.unit_vector_fields(dist)
lift_basis = zbasis.derivative_basis(1)
lift = lambda A: d3.Lift(A, lift_basis, -1)

B_op = (np.absolute(D - M - N_s2*Z)+ M + D - N_s2*Z)/2

Max = lambda A,B: (abs(A-N_s2*Z-B)+A-N_s2*Z+B)/2
eva = lambda A: A.evaluate()

dz= lambda A: d3.Differentiate(A, coords['z'])
dx= lambda A: d3.Differentiate(A, coords['x'])

ux=u@ex
uz=u@ez
dxux=dx(ux)
dzux=dz(ux)
dxuz=dx(uz)
dzuz=dz(uz)

ux2=ux*ux
uy2=uy*uy
uz2=uz*uz

grad_u = d3.grad(u) + ez* lift(tau_u1) # First-order reduction
grad_ux = grad_u@ex # First-order reduction
grad_uz = grad_u@ez # First-order reduction
grad_uy = d3.grad(uy) + ez*lift(tau_u3)# First-order reduction 
grad_M = d3.grad(M) + ez*lift(tau_M1) # First-order reduction
grad_D = d3.grad(D) + ez*lift(tau_D1) # First-order reduction



            


def draw(taskname, name1, name2, folder_dir):

    file_paths = [os.path.join(folder_dir, file) for file in listdir(folder_dir) if os.path.isfile(os.path.join(folder_dir, file)) and file.endswith('.h5')]
    #sort by the number in the file name
    file_paths.sort(key=lambda f: int(re.sub('\D', '', f)))
    print(file_paths)
    recorded = False
    for file in file_paths:
        with h5py.File(file, mode='r') as file:
            alldata = file['tasks'][f'{taskname}']
            if recorded == False:
                max_level_old = np.max(alldata)
                min_level_old = np.min(alldata)
                recorded = True
            else:
                max_level_new = np.max(alldata)
                min_level_new = np.min(alldata)
                if max_level_new > max_level_old:
                    max_level_old = max_level_new
                if min_level_new < min_level_old:
                    min_level_old = min_level_new
    min_lev, max_lev = min_level_old, max_level_old
    levels = np.arange(min_lev-(max_lev-min_lev)/32, max_lev+(max_lev-min_lev)/32, (max_lev-min_lev)/32)
    n=0
    for file in file_paths:
        with h5py.File(file, mode='r') as file:
            task = file['tasks'][f'{taskname}']
            st = file['scales/sim_time']
            simtime = np.array(st)
            for t in range(0, len(simtime)):
                task_T = task[t,:,:].T
                
                n_rows, n_columns = task_T.shape
                x_ax = np.linspace(0, Lx, n_columns)
                y_ax = np.linspace(0, Lz, n_rows)
                X_ax, Y_ax = np.meshgrid(x_ax, y_ax)
                
                plt.contourf(X_ax, Y_ax, task_T, levels, cmap='Spectral_r')
                plt.colorbar(label=f'{name1}')
                plt.xlabel('x')
                plt.ylabel('z')
                n=n+1
                # Add time title
                title = "t="+str(st[t])
                plt.title(title)
                plt.savefig(f'{name1}/{name2}_{"%04d" % n}.png', dpi=200,bbox_inches='tight')
                matplotlib.pyplot.close()   
def draw2(taskname, name1, name2, folder_dir):
    file_paths = [os.path.join(folder_dir, file) for file in listdir(folder_dir) if os.path.isfile(os.path.join(folder_dir, file)) and file.endswith('.h5')]
    #sort by the number in the file name
    file_paths.sort(key=lambda f: int(re.sub('\D', '', f)))
    n=0
    for file in file_paths:
        with h5py.File(file, mode='r') as file:
            task = file['tasks'][f'{taskname}']
            st = file['scales/sim_time']
            simtime = np.array(st)
            for t in range(0, len(simtime)):
                task_T = task[t,:,:].T
                n_rows, n_columns = task_T.shape
                x_ax = np.linspace(0, Lx, n_columns)
                y_ax = np.linspace(0, Lz, n_rows)
                X_ax, Y_ax = np.meshgrid(x_ax, y_ax)  
                
                #min_lev, max_lev = np.min(task_T), np.max(task_T)
                #levels = np.arange(min_lev-(max_lev-min_lev)/32, max_lev+(max_lev-min_lev)/32, (max_lev-min_lev)/32)
                min_lev, max_lev = np.min(task_T), np.max(task_T)
                if max_lev - min_lev <= 0.1:
                    plt.contourf(X_ax, Y_ax, task_T, cmap='RdBu_r')
                else:
                    levels = np.arange(min_lev-(max_lev-min_lev)/16, max_lev+(max_lev-min_lev)/16, (max_lev-min_lev)/16)
                    plt.contourf(X_ax, Y_ax, task_T, levels, cmap='RdBu_r')
                plt.colorbar(label=f'{name1}')
                plt.xlabel('x')
                plt.ylabel('z')
                n=n+1
                # Add time title
                title = "t="+str(st[t])
                plt.title(title)
                plt.savefig(f'{name1}/{name2}_{"%04d" % n}.png', dpi=200,bbox_inches='tight')
                matplotlib.pyplot.close()   


if not os.path.exists('liquid_water'):    
    os.mkdir('liquid_water')           
draw("liquid water", "liquid_water", "liquid_water", "analysis_2")            
            

if not os.path.exists('moist buoyancy'):
    os.mkdir('moist buoyancy')
draw("M", "moist buoyancy", "mb", "snapshots")


if not os.path.exists('rotational_uy'):
    os.mkdir('rotational_uy')
draw2("uy", "rotational_uy", "rot_uy", "snapshots")
    

folder_dir = "snapshots"

file_paths = [os.path.join(folder_dir, file) for file in listdir(folder_dir) if os.path.isfile(os.path.join(folder_dir, file)) and file.endswith('.h5')]
#sort by the number in the file name
file_paths.sort(key=lambda f: int(re.sub('\D', '', f)))
print(file_paths)

n = 0
if not os.path.exists('rotational_buoyancy'):
    os.mkdir('rotational_buoyancy')

recorded = False
for file in file_paths:
    with h5py.File(file, mode='r') as file:
        moistbuoyancy = file['tasks']['M']
        drybuoyancy = file['tasks']['D']
        buoyancy = np.maximum(moistbuoyancy, drybuoyancy - N_s2 * z)
        if recorded == False:
            max_level_old = np.max(buoyancy)
            min_level_old = np.min(buoyancy)
            recorded = True
        else:
            max_level_new = np.max(buoyancy)
            min_level_new = np.min(buoyancy)
            if max_level_new > max_level_old:
                max_level_old = max_level_new
            if min_level_new < min_level_old:
                min_level_old = min_level_new
min_lev, max_lev = min_level_old, max_level_old
levels = np.arange(min_lev-(max_lev-min_lev)/32, max_lev+(max_lev-min_lev)/32, (max_lev-min_lev)/32)

for file in file_paths:
    with h5py.File(file, mode='r') as file:
        moistbuoyancy = file['tasks']['M']
        drybuoyancy = file['tasks']['D']
        buoyancy = np.maximum(moistbuoyancy, drybuoyancy - N_s2 * z)
        st = file['scales/sim_time']
        simtime = np.array(st)
        for t in range(0, len(simtime)):
            mb = np.transpose(moistbuoyancy[t, :, :])
            db = np.transpose(drybuoyancy[t, :, :])
            bu = np.transpose(buoyancy[t, :, :])
            
            n_rows, n_columns = bu.shape
            x_ax = np.linspace(0, Lx, n_columns)
            y_ax = np.linspace(0, Lz, n_rows)
            X_ax, Y_ax = np.meshgrid(x_ax, y_ax) 
            
            plt.contourf(X_ax, Y_ax, bu, levels, cmap='Spectral_r')
            plt.colorbar(label='buoyancy')
            plt.xlabel('x')
            plt.ylabel('z')
            n = n + 1
            # Add time title
            title = "t=" + str(simtime[t])
            plt.title(title)
            plt.savefig(f'rotational_buoyancy/rot_{"%04d" % n}.png', dpi=200, bbox_inches='tight')
            matplotlib.pyplot.close()


