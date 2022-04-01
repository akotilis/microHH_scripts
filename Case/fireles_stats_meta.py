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
sql = thermo.variables['ql'][:, :]          #Liquid water [kg kg-1]
sthl = thermo.variables['thl'][:, :]        #Liquid water potential temperature
sqt = thermo.variables['qt'][:, :]          #Total water mixing ratio
areat = default.variables['area'][:, :]
sqlpath = thermo.variables['ql_path'][:]    #[kg m-2]
area = np.mean(areat[:,:], 0)
cft = thermo.variables["ql_frac"][:]
snr = thermo.variables["nr_path"][:]
trr = thermo.variables["rr"][:]
sqr = thermo.variables["qr_path"][:]		#rain water path [kg m-2]
ql = np.sum(areat[:,:]*sql[:,:], 0) / np.sum(areat[:,:], 0)
cf = np.sum(areat[:,:]*cft[:,:], 0) / np.sum(areat[:,:], 0)

#Changing units
sql = sql*10**6 #to g m-3
trr = trr*10**3 #to mm/s
sqlpath = sqlpath*10**3 #to g m-2
sqr = sqr*10**3 #to g m-2

#conditional ql>0
cft[cft==0] = np.nan
cf[cf==0] = np.nan
sql[sql==0] = np.nan

start, end = 0, len(t)



#mean and std over the horizontal domain(time) for errorbar
sql_error = np.nanstd(sql[:,:], 0)	
sql_mean = np.nanmean(sql[:,:], 0)

#Finds the index of the max value 
ind = np.unravel_index(np.nanargmax(sql_mean), sql_mean.shape)


#Plot
plt.figure(figsize=(14,10))
#plt.figure(figsize=(10,5))
plt.subplot(221)
plt.plot(t, sqlpath, label = 'liquid water path')
#plt.plot(t, sqlpath1d, 'b:', label='interpolated')
plt.xlabel('time[s]')
plt.ylabel('q$_l$ [g m$^{-2}$]')
plt.title('LWP')
plt.legend()

#plt.figure()
plt.subplot(222)
#for n in range(start,end):
#	plt.plot(sql[n,:], z, color='#eeeeee')
plt.plot(sql_mean, z, label='mean')
#plt.plot(np.std(sql[0:-1, :], 0),z, label='std')
plt.xlabel('q$_l$ [g m$^{-3}$]')
plt.ylabel('z [m]')
plt.title('Liquid water')
plt.legend()

#plt.figure(figsize=(10,5))
plt.subplot(223)
plt.plot(t, sqr, label = 'rain water path')
plt.xlabel('time[s]')
plt.ylabel('q$_l$ [g m$^{-2}$]')
plt.title('RWP')
plt.legend()


plt.subplot(224)
#for n in range(start,end):
#	plt.plot(cft[n,:],z, color='#eeeeee')
plt.plot(cf, z)
plt.xlabel(r'cloud fraction [-]')
plt.ylabel(r'z [m]')



print("Highest value of LWP is {:.3f} g m-2 ".format(np.max(sqlpath)))
print("Highest value of RWP is {:.5f} g m-2 ".format(np.max(sqr)))
print("Highest value of ql is {:.2f} g m-3 at z={:g} m".format(np.nanmax(sql_mean), z[ind]))

'''
#plt.figure(figsize=(10,5))
plt.subplot(224)
plt.plot(t, trr, 'b-', label = 'mean surface rain rate')
plt.xlabel('time[s]')
plt.ylabel('mean surface rain rate [mm s$^{-1}$]')
plt.title('Mean surface rain rate')
plt.legend()

plt.figure()
plt.errorbar(sql_mean, z, xerr = sql_error, fmt='-', errorevery=(6,3))
plt.xlabel('q$_l$ [g m$^{-3}$]')
plt.ylabel('z [m]')
plt.title('Liquid water')


plt.figure(figsize=(10,5))
plt.plot(t, snr, 'b-', label = 'number density rain')
plt.xlabel('time[s]')
plt.ylabel('n$_r$ [kg m$^{-5}$]')
plt.title('NDR')
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


plt.figure()
for n in range(start,end):
	plt.plot(area[n,:]*100,z, color='#eeeeee')
plt.plot(aream*100, z, label = 'mean={mean:.2f}% , 
			max={max_value:.2f%'.format(mean=np.mean(aream)*100, 			max_value=np.max(aream)*100))
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

