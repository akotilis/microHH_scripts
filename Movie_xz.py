import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable

save_results = '/home/andreas/MicroHH/Results/'

file = "qlpath.xz.nc"
with nc.Dataset(file) as nc:
    print(nc.variables)
    s = nc.variables['qlpath'][:]
    x = nc.variables['x'][:]
    y = nc.variables['y'][:]
    z = nc.variables['z'][:]
    #z=10
    time = nc.variables['time'][:]
    


for i in range(len(s)): 
    fig, ax = plt.subplots(figsize=(12,8))
    X,Y = np.meshgrid(x,z)
    # level = 0,y
    data = s[i,:,0]
    vmin, vmax = np.ma.min(s), np.ma.max(s)
    level_boundaries = np.linspace(vmin,vmax,100)
    cs = ax.contourf(X,Y, data, levels=level_boundaries, vmin = vmin, vmax= vmax, extend='both',cmap='RdBu_r')
    fig.colorbar(cs)
    ax.set_title('xz cross-section of qlpath, t={time:.3f}s'.format(time=time[i]))
    ax.set_xlabel("x domain")
    ax.set_ylabel("z domain")
    fig.savefig(save_results + 'qlpath%4.4i.png'%(i), dpi='figure')
    plt.close(fig)
