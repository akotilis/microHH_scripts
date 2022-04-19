import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt


stats = nc.Dataset('fireles.default.0000000.nc', 'r')
default = stats.groups['default']
thermo = stats.groups['thermo']

t = stats.variables['time'][:]
z = stats.variables['z'][:]
zh = stats.variables['zh'][:]

thlf = thermo.variables['thl_flux'][:, :]  #flux of the Liq water potential temperature[K m s-1]
qtf = thermo.variables['qt_flux'][:,:] 	#flux of water mixing ratio [kg kg-1 m s-1]
rh = thermo.variables['rh'][:, :]		#relative humidity
sql = thermo.variables['ql'][:, :]          #Liquid water [kg kg-1]
thl = thermo.variables['thl'][:, :]        #Liquid water potential temperature
sqt = thermo.variables['qt'][:, :]          #Total water mixing ratio
areat = default.variables['area'][:, :]
sqlpath = thermo.variables['ql_path'][:]    #[kg m-2]
cft = thermo.variables["ql_frac"][:]
snr = thermo.variables["nr"][:, :]		#number density rain [m-3]
trr = thermo.variables["rr"][:]
sqr = thermo.variables["qr_path"][:]		#rain water path [kg m-2]
bld = thermo.variables['zi'][:]		#boundary layer depth [m]
qtbot = thermo.variables['qt_bot'][:]		#surface total water mixing ratio [kg kg-1]
thlbot = thermo.variables['thl_bot'][:] 	#surface potential temperature [K]
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
snr[snr==0] = np.nan
thlf[thlf==0] = np.nan


start, end = 0, len(t)
step = 12

#mean and std over the horizontal axis(time)
thl_mean = np.mean(thl[:,:], 0)
sql_error = np.nanstd(sql[:,:], 0)	
sql_mean = np.nanmean(sql[:,:], 0)
snr_mean = np.nanmean(snr[:,:], 0)
thlf_mean = np.nanmean(thlf[:,:], 0)
rh_mean = np.mean(rh[:,:], 0)

#Finds the index of the max value 
ind = np.unravel_index(np.nanargmax(sql_mean), sql_mean.shape)
ind1 = np.unravel_index(np.argmax(bld), bld.shape)

#max BL height and normalizing
zi = bld[ind1]
normz = z/zi
normzh= zh/zi

norm_thlf = thlf/thlf[0,0]		#normalized kinematic heat flux (thl_bot=0.3)

flux_check=np.zeros([end])
for n in range(start,end):
	flux_check[n] = (thlf[n,1]/zi)*27000
	


plt.figure()
plt.plot(np.mean(norm_thlf[:,:],0), normzh, label='normalized kinematic heat flux')
plt.ylabel('z/zi')
plt.xlabel('$θ/Qs$ [K]')
plt.title(' normalized kinematic heat flux')
plt.legend()

#hourly development of the potential temperature profile
plt.figure()
for n in range(start, 73, step):
	plt.plot(thl[n,:], z, label = 't={:g}h'.format(t[n]/3600))
plt.ylabel('z [m]')
plt.xlabel('$θ$ [K]')
plt.title('Hourly Potential temperature profile')
plt.legend()

plt.figure()
for n in range(start, end, step):
	plt.plot(norm_thlf[n,:], normzh, label = 't={:g}h'.format(t[n]/3600))
#plt.plot(norm_thlf[73,:], normzh, label= 't={:g}h'.format(t[73]/3600))
plt.ylabel('z/zi')
plt.xlabel('$θ/Qs$')
plt.title('Hourly normalized kinematic heat flux')
plt.legend()

plt.figure()
plt.plot(t, thl[:,0], label='50m temp')
plt.xlabel('time [s]')
plt.ylabel('Temperature [K]')
plt.legend()

plt.figure()
plt.plot(t, thlf[:,1], label='100m height heat flux')
plt.xlabel('time [s]')
plt.ylabel("w'θ' [K m s$^{-1}$]")
plt.legend()

plt.figure()
plt.plot(t, flux_check, label='100m height heat flux validation')
plt.xlabel('time [s]')
plt.ylabel("w'θ' [K m s$^{-1}$]")
plt.legend()


gamma = 0.0004
bowr= gamma*thlbot/qtbot
bownr = gamma*thlf[:,0]/qtf[:,0]
bowrm = np.mean(bownr)

plt.figure()
plt.plot(t, bowr, label='surface bowen ratio')
plt.xlabel('time [s]')
plt.ylabel('Bowen ratio')
plt.legend()
'''
plt.figure()
plt.plot(t, bowr[:,0], label='50m')
plt.plot(t, bowr[:,1], label='100m')
plt.plot(t, bowr[:,10], label='1000m')
plt.xlabel('time [s]')
plt.ylabel('Bowen ratio(-)')
plt.legend()


#Plot
plt.figure(figsize=(14,10))
#plt.figure(figsize=(10,5))
plt.subplot(221)
plt.plot(t, sqlpath, label = 'liquid water path')
#plt.plot(t, sqlpath1d, 'b:', label='interpolated')
plt.xlabel('time[s]')
plt.ylabel('LWP [g m$^{-2}$]')
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
plt.ylabel('q$_r$ [g m$^{-2}$]')
plt.legend()


plt.subplot(224)
#for n in range(start,end):
#	plt.plot(cft[n,:],z, color='#eeeeee')
plt.plot(cf, z)
plt.xlabel('cloud fraction [-]')
plt.ylabel('z [m]')

plt.figure()
plt.plot(thlf_mean, zh, label ='flux of temperature')
plt.xlabel('Thl_flux [K m s$^{-1}$]')
plt.ylabel('z[m]')
plt.legend()






plt.figure()
plt.plot(t, bld, label ='boundary layer depth')
plt.xlabel('time[s]')
plt.ylabel('zi [m]')
plt.legend()

plt.figure()
plt.plot(rh_mean, z)

plt.figure()
plt.plot(t, np.mean(rh[:,:], 1))

fig, ax = plt.subplots()
ax.plot(t, bld, label ='boundary layer depth', color='black')
ax.set_xlabel('time [s]')
ax.set_ylabel('zi [m]')
ax2 = ax.twinx()
ax2.plot(t, np.nanmean(thlf[:,:], 1), label='temperature flux')
ax2.set_ylabel('Thl_flux [K m s$^{-1}$]')
ax.grid()





print("Average value of LWP is {:.3f} g m-2 ".format(np.mean(sqlpath)))
print("Average value of RWP is {:e} g m-2 ".format(np.mean(sqr)))
print("Average value of ql is {:.2f} g m-3 with the highest amount at z={:g} m".format(np.nanmean(sql_mean), z[ind]))
print("Average cloud cover is {:.2f}%".format(np.nanmean(cf)*100))

plt.subplot(221)
plt.plot(snr_mean, z, label = 'number density rain')
plt.xlabel('n$_r$ [m$^{-3}$]')
plt.ylabel('z [m]')
plt.title('Nr')
plt.legend()

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

