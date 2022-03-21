import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt


stats = nc.Dataset('fireles.default.0000000.nc', 'r')
default = stats.groups['default']
thermo = stats.groups['thermo']

t = stats.variables['time'][:]
z = stats.variables['z'][:]
zh = stats.variables['zh'][:]

st = thermo.variables['thv'][:, :]
u2t = default.variables['u_2'][:, :]
v2t = default.variables['v_2'][:, :]
w2t = default.variables['w_2'][:, :]
s2t = thermo.variables['thv_2'][:, :]
sql = thermo.variables['ql'][:, :]          #Liquid water
sthl = thermo.variables['thl'][:, :]        #Liquid water potential temperature
sqt = thermo.variables['qt'][:, :]          #Total water mixing ratio
areat = default.variables['area'][:, :]
sqlpath = thermo.variables['ql_path'][:]
area = np.mean(areat[:,:], 0)
cft = thermo.variables["ql_frac"][:]
snr = thermo.variables["nr_path"][:]
ql = np.sum(areat[:,:]*sql[:,:], 0) / np.sum(areat[:,:], 0)
cf = np.sum(areat[:,:]*cft[:,:], 0) / np.sum(areat[:,:], 0)





plt.close('all')
# enable LaTeX plotting
#plt.rc('font',**{'family':'serif','serif':['Palatino']})
#plt.rc('text', usetex=True)
start, end = 0, len(t)
cft[cft==0] = np.nan
cf[cf==0] = np.nan
sql[sql==0] = np.nan

plt.figure()
for n in range(start,end):
	plt.plot(cft[n,:],z, color='#eeeeee')
plt.plot(cf, z)
plt.xlabel(r'cloud fraction [-]')
plt.ylabel(r'z [m]')

plt.figure()
for n in range(start,end):
	plt.plot(sql[n,:], z, color='#eeeeee')
plt.plot(np.nanmean(sql[:, :], 0),z, label='mean')
#plt.plot(np.std(sql[0:-1, :], 0),z, label='std')
plt.xlabel('q$_l$ [g~kg$^{-1}$]')
plt.ylabel('z [m]')
plt.title('Liquid water')
plt.legend()


plt.figure(figsize=(10,5))
#for n in range(start,end):
#	plt.plot(t[:,n], sqlpath[n,:], color='#eeeeee')
plt.plot(t, sqlpath, 'b-', label = 'liquid water')
#plt.plot(t, sqlpath1d, 'b:', label='interpolated')
plt.xlabel('time[s]')
plt.ylabel('q$_l$ [g~kg$^{-1}$]')
plt.title('LWP')
plt.legend()

plt.figure(figsize=(10,5))
plt.plot(t, snr, 'b-', label = 'number density rain')
plt.xlabel('time[s]')
plt.ylabel('n$_r$ [g~kg$^{-1}$]')
plt.title('NDR')
plt.legend()

'''
plt.figure()
plt.plot(np.mean(sqt[60:-1, :], 0), z, 'b-', label='mean over t=18000-21600s')
#plt.plot(np.std(sqt[0:-1, :], 0), z, 'r--',label='std t=18000-21600s')
plt.plot(sqt[60,:], z, 'r-', alpha=0.3, label='t={time:g}s'.format(time=t[60]))
plt.plot(sqt[-1,:], z, 'r--', alpha=0.3, label='t={time:g}s'.format(time=t[-1]))
plt.title("Total water mixing ratio")
plt.xlabel(r'$q_t [kg kg^{-1}]$')
plt.ylabel(r'$z [m]$')
plt.legend()

plt.figure()
plt.plot(np.mean(sthl[60:-1, :], 0), z, 'b-', label='mean over t=18000-21600s')
plt.plot(sthl[60,:], z, 'r-', alpha=0.3, label='t={time:g}s'.format(time=t[60]))
plt.plot(sthl[-1,:], z, 'r--', alpha=0.3, label='t={time:g}s'.format(time=t[-1]))
plt.title("Potential temperature")
plt.xlabel(r'$\theta [K]$')
plt.ylabel(r'$z [m]$')
plt.legend()


plt.figure()
for n in range(start,end):
	plt.plot(area[n,:]*100,z, color='#eeeeee')
plt.plot(aream*100, z, 
		label = 'mean={mean:.2f}% , max={max_value:.2f}%'.format(mean=np.mean(aream)*100, 			max_value=np.max(aream)*100))
plt.xlabel('area coverage [%]')
plt.ylabel('z [m]')
plt.title('Cloud coverage area')
plt.legend()

#nans, x = nan_help(sqlpath1d)
#sqlpath1d[nans] = np.interp(x(nans), x(~nans), sqlpath1d[~nans])

#ht = np.zeros(t.size)
#wstart = np.zeros(t.size)
#for n in range(t.size):
#    hindex = np.argmin(abs(sgradt[n, :] - max(sgradt[n, :])))
#    ht[n] = z[hindex]
#    wstart[n] = ((9.81/300.)*sfluxt[n, 0]*ht[n])**(1./3.)

'''
plt.show()

