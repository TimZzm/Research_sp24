import h5py
import numpy as np
import matplotlib.pyplot as plt

def get_bar_range(paths, taskname):
    recorded = False
    for file in paths:
        with h5py.File(file, mode='r') as file:
            if recorded == False:
                max_level_old = np.max(file['tasks'][taskname])
                min_level_old = np.min(file['tasks'][taskname])
                recorded = True
            else:
                max_level_new = np.max(file['tasks'][taskname])
                min_level_new = np.min(file['tasks'][taskname])
                if max_level_new > max_level_old:
                    max_level_old = max_level_new
                if min_level_new < min_level_old:
                    min_level_old = min_level_new
    return levels = np.arange(min_level_old, max_level_old, (max_level_old - min_level_old) / 40)


def produce_image(paths, name):
    if not os.path.exists(name):    
        os.mkdir(name)
    n=0
    for file in paths:
    with h5py.File(file, mode='r') as file:
        moistbuoyancy = file['tasks']['moist buoyancy']
        drybuoyancy = file['tasks']['dry buoyancy']
        ql=np.maximum(moistbuoyancy[:]-drybuoyancy[:]+N_s2*z,0)
        st = file['scales/sim_time']
        simtime = np.array(st)
        for t in range(0, len(simtime)):
            qli=np.transpose(ql[t,:,:])
            plt.contourf(qli, cmap='RdBu_r')
            plt.colorbar(label='liquid water')
            plt.xlabel('x')
            plt.ylabel('z')
            n=n+1
            # Add time title
            title = "t="+str(st[t])
            plt.title(title)
            plt.savefig(f'liquid water/liquidwater_{"%04d" % n}.png', dpi=200,bbox_inches='tight')
            matplotlib.pyplot.close() 












    


