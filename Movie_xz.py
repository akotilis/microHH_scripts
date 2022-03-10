import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable

save_results = '/home/andreas/MicroHH/Results/fireles_standard_case/s1_xz'

file = "s1.xz.nc"
with nc.Dataset(file) as nc:
    print(nc.variables)
    s = nc.variables['s1'][:]
    x = nc.variables['x'][:]
    y = nc.variables['y'][:]
    z = nc.variables['z'][:]
    time = nc.variables['time'][:]
    


for i in range(len(s)): 
    fig, ax = plt.subplots()
    X,Y = np.meshgrid(x,z)
    # level = 0,y
    data = s[i,:,0]
    vmin, vmax = np.ma.min(s), np.ma.max(s)
    level_boundaries = np.linspace(vmin,vmax,200)
    cs = ax.contourf(X,Y, data, levels=level_boundaries, vmin = vmin, vmax= vmax, extend='both',cmap='RdBu_r')
    fig.colorbar(cs)
    ax.set_title('xz cross-section of s1, t={time:.1f}h'.format(time=time[i]/3600))
    ax.set_xlabel("x domain")
    ax.set_ylabel("z domain")
    fig.savefig(save_results + 's1%4.4i.png'%(i), dpi='figure')
    plt.close(fig)
