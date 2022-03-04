import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
save_results = '/home/andreas/MicroHH/Results/'

file = "qlpath.xy.nc"
with nc.Dataset(file) as ncf:
	print(ncf)
	#th = ncf.variables['thl'][:]
	qlpath = ncf.variables['qlpath'][:]
	x = ncf.variables['x'][:]
	y = ncf.variables['y'][:]
	time = ncf.variables['time'][:]
	z = 10
	#vmin = 302.5
	#vmax = 303
	vmin, vmax = np.ma.min(qlpath), 0.9040#np.ma.max(qlpath)
	level_boundaries = np.linspace(vmin,vmax,200)
for i, vel in enumerate(qlpath):
		
    fig, ax = plt.subplots()
    X,Y = np.meshgrid(x,y)
    vmin, vmax = np.ma.min(qlpath), np.ma.max(qlpath) 
    level_boundaries = np.linspace(vmin,vmax,200)
    data = qlpath[i,:,:]
    cs = ax.contourf(X,Y, data, levels=level_boundaries, cmap = 'RdBu_r', vmin=vmin, vmax=vmax, extend='both')
    ax.set_title('qlpath t={time:.3f}s'.format(time=time[i]))
    ax.set_ylabel('y size of domain[m]')
    ax.set_xlabel('x size of domain[m]')
    # 		fig.colorbar(ScalarMappable(norm=cs.norm, cmap=cs.cmap),
    #     			     ticks=np.linspace(vmin, vmax, 6),
    #                             boundaries=level_boundaries,
    #                             values=(level_boundaries[:-1] + level_boundaries[1:]) / 2)
    fig.colorbar(cs)
    fig.savefig(save_results + 'qlpath%4.4i.png'%(i), dpi = 'figure')
    plt.close(fig)	
	#plt.show()
