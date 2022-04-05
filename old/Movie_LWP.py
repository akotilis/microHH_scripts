import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
save_results = '/home/andreas/MicroHH/Results/fireles_outflow/qlpath/'

file = "qlpath.xy.nc"
with nc.Dataset(file) as ncf:
	print(ncf)
	#th = ncf.variables['thl'][:]
	s = ncf.variables['qlpath'][:]
	x = ncf.variables['x'][:]
	y = ncf.variables['y'][:]
	time = ncf.variables['time'][:]
	
	vmin, vmax = np.ma.min(s), np.ma.max(s)
	level_boundaries = np.linspace(vmin,vmax,200)
for i, vel in enumerate(s):
		
    fig, ax = plt.subplots()
    X,Y = np.meshgrid(x,y)
    vmin, vmax = np.ma.min(s), np.ma.max(s) 
    level_boundaries = np.linspace(vmin,vmax,200)
    data = s[i,:,:]

    cs = ax.contourf(X,Y, data, levels=level_boundaries, cmap = 'RdBu_r', vmin=vmin, vmax=vmax, extend='both')
    ax.set_title('xy cross-section of LWP t={time:.1f}h'.format(time=time[i]/3600))
    ax.set_ylabel('y size of domain[m]')
    ax.set_xlabel('x size of domain[m]')
    
    fig.colorbar(cs)
    fig.savefig(save_results + 'LWP_outflow_%4.4i.png'%(i), dpi = 'figure')
    plt.close(fig)	
	#plt.show()
