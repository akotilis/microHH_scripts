import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

start = 0
end = 91
step = 2

stats = nc.Dataset('fireles.default.0000000.nc', 'r')
default = stats.groups['default']
thermo = stats.groups['thermo']

t = stats.variables['time'][:]
z = stats.variables['z'][:]
zh = stats.variables['zh'][:]

st = thermo.variables['thv'][:, :]
evisct = default.variables['evisc'][:, :]   #Eddy viscosity
u2t = default.variables['u_2'][:, :]
v2t = default.variables['v_2'][:, :]
w2t = default.variables['w_2'][:, :]
s2t = thermo.variables['thv_2'][:, :]
sturbt = thermo.variables['thl_w'][:, :]     #Turbulent flux of the Liquid water potential temperature
sdifft = thermo.variables['thv_diff'][:, :]  #Diffusive flux of the Liquid water potential temperature
sfluxt = thermo.variables['thv_flux'][:, :]  #Total flux of the Liquid water potential temperature
sgradt = thermo.variables['thv_grad'][:, :]  #Gradient of the Liquid water potential temperature
sql = thermo.variables['ql'][:, :]*1000          #Liquid water
sthl = thermo.variables['thl'][:, :]        #Liquid water potential temperature
sqt = thermo.variables['qt'][:, :]          #Total water mixing ratio
area = default.variables['area'][:, :]
sqlpath = thermo.variables['ql_path'][:]
#area = np.mean(areat[:,:], 0)
cft = thermo.variables["ql_frac"][:]
snr = thermo.variables["nr_path"][:]

aream = np.mean(area[:,:], 0)

ql = np.sum(area[:,:]*sql[:,:], 0) / np.sum(area[:,:], 0)
cf = np.sum(area[:,:]*cft[:,:], 0) / np.sum(area[:,:], 0)
areaql = area[:,:]*sql[:,:]


s = np.mean(st, axis=0)
evisc = np.mean(evisct, axis=0)

u2 = np.mean(u2t, axis=0)
v2 = np.mean(v2t, axis=0)
w2 = np.mean(w2t, axis=0)
s2 = np.mean(s2t, axis=0)

sturb = np.mean(sturbt, axis=0)
sdiff = np.mean(sdifft, axis=0)
sflux = np.mean(sfluxt, axis=0)

sqlpath1d = np.array(sqlpath).flatten()

def nan_help(y):
	return np.isnan(y), lambda z:z.nonzero()[0]
	
#nans, x = nan_help(sqlpath1d)
#sqlpath1d[nans] = np.interp(x(nans), x(~nans), sqlpath1d[~nans])

#ht = np.zeros(t.size)
#wstart = np.zeros(t.size)
#for n in range(t.size):
#    hindex = np.argmin(abs(sgradt[n, :] - max(sgradt[n, :])))
#    ht[n] = z[hindex]
#    wstart[n] = ((9.81/300.)*sfluxt[n, 0]*ht[n])**(1./3.)



#area[area==0] =  np.nan
#aream[aream==0] = np.nan
plt.close('all')
# enable LaTeX plotting
#plt.rc('font',**{'family':'serif','serif':['Palatino']})
#plt.rc('text', usetex=True)

plt.figure()
for n in range(start,end):
	plt.plot(area[n,:]*100,z, color='#eeeeee')
plt.plot(aream*100, z, 
		label = 'mean={mean:.2f}% , max={max_value:.2f}%'.format(mean=np.mean(aream)*100, 			max_value=np.max(aream)*100))
plt.xlabel('area coverage [%]')
plt.ylabel('z [m]')
plt.title('Cloud coverage area')
plt.legend()

plt.figure()
for n in range(start,end):
	plt.plot(sql[n,:], z, color='#eeeeee')
plt.plot(np.mean(sql[0:-1, :], 0),z, label='mean')
#plt.plot(np.std(sql[0:-1, :], 0),z, label='std')
plt.xlabel('q$_l$ [g~kg$^{-1}$]')
plt.ylabel('z [m]')
plt.title('Liquid water')
plt.legend()


plt.figure(figsize=(10,5))
#for n in range(start,end):
#	plt.plot(t[:,n], sqlpath[n,:], color='#eeeeee')
plt.plot(t, sqlpath, 'b-', label = 'liquid water time series')
#plt.plot(t, sqlpath1d, 'b:', label='interpolated')
plt.xlabel('time[s]')
plt.ylabel('q$_l$ [g~kg$^{-1}$]')
plt.title('')
plt.legend()

plt.figure(figsize=(10,5))
plt.plot(t, snr, 'b-', label = 'number density rain')
plt.xlabel('time[s]')
plt.ylabel('n$_r$ [g~kg$^{-1}$]')
plt.title('NDR time series')
plt.legend()


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

'''
plt.figure()
for n in range(start,end):
	plt.plot(cft[n,:],z, color='#eeeeee')
plt.plot(cf, z)
plt.xlabel(r'cloud fraction [-]')
plt.ylabel(r'z [m]')

plt.figure()
for n in range(start, end, step):
    plt.plot(st[n, :], z)
plt.xlabel(r'$\theta [K]$')
plt.ylabel(r'$z [m]$')

plt.figure()
for n in range(start, end, step):
    plt.plot(evisct[n, :], z)
plt.xlabel(r'$K_m [m^2 s^{-1}]$')
plt.ylabel(r'$z [m]$')
 
plt.figure()
for n in range(start, end, step):
    plt.plot(u2t[n, :], z)
plt.xlabel(r'$u^2 [m^2 s^{-2}]$')
plt.ylabel(r'$z [m]$')

plt.figure()
for n in range(start, end, step):
    plt.plot(w2t[n, :], zh)
plt.xlabel(r'$w^2 [m^2 s^{-2}]$')
plt.ylabel(r'$z [m]$')

plt.figure()
for n in range(start, end, step):
    plt.plot(sfluxt[n, :], zh)
plt.xlabel(r'$\overline{w\theta} [K m s^{-1}]$')
plt.ylabel(r'$z [m]$')


plt.figure()
for n in range(start, end, step):
    plt.plot(sqt[n, :], z)
plt.xlabel(r'$q_t [kg kg^{-1}]$')
plt.ylabel(r'$z [m]$')

plt.figure()
for n in range(start, end, step):
    plt.plot(sthl[n, :], z)
plt.xlabel(r'$\theta [K]$')
plt.ylabel(r'$z [m]$')

plt.figure()
for n in range(start, end, step):
    plt.plot(sql[n, :], z)
plt.xlabel(r'$q_l [kg kg^{-1}]$')
plt.ylabel(r'$z [m]$')

plt.figure()
for n in range(start, end, step):
    plt.plot(area[n, :], z)
plt.xlabel(r'$Area [\%]$')
plt.ylabel(r'$z [m]$')
'''


'''
plt.figure()

plt.plot(np.mean(sfluxt[-5:-1, :], 0), zh, 'b-')
plt.plot(np.mean(sturbt[-5:-1, :], 0), zh, 'b--')
plt.plot(np.mean(sdifft[-5:-1, :], 0), zh, 'b:')
plt.ylim(0., 1500.)
plt.xlabel(r'$\overline{w\theta} [K m s^{-1}]$')
plt.ylabel(r'$z [m]$')
'''

#plt.fill_between(sthl[61,:], sthl[72,:], where = (sthl[61,:]<sthl[72,:]))


#plt.figure(2)
#for n in range(61, end, step):
#    plt.plot(sthl[n, :], z)
#plt.xlabel(r'$\theta [K]$')
#plt.ylabel(r'$z [m]$')


'''
plt.figure()
plt.plot(t, ht)
plt.xlabel(r'$time [s]$')
plt.ylabel(r'$h [m]$')
'''
plt.show()
