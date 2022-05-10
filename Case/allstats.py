import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import glob
import pandas as pd

#def sql(cs):
#      stats = nc.Dataset('thl/fireles.default.0000000_'+cs+'.nc', 'r')
#      z = stats.variables['z'][:]
#      sql = stats.groups['thermo'].variables['ql'][:]*10e5 	#to g m-3
#      sql[sql==0] = np.nan
#      return np.nanmean(sql[:,:],0), z
      
	
# Record all case names of the nc files into a list
all_files=[]

for file in glob.glob('thl/*.nc', recursive=True):
	nc_file = nc.Dataset(file,'r')
	case = file[28:33]
	all_files.append(case)
	print(file)
	
all_files.sort()

all_files_qt=[]

for file in glob.glob('qt/*.nc', recursive=True):
	nc_file = nc.Dataset(file,'r')
	case = file[27:31]
	all_files_qt.append(case)
	print(file)
	
all_files_qt.sort()

print(all_files)
print(all_files_qt)

# Storing lists
ql_thl = []
cf_thl = []
th_thl = []
ql_qt = []
cf_qt = []

# Loop over the nc files 
for i in all_files:
      
      stats = nc.Dataset('thl/fireles.default.0000000_{}.nc'.format(i), 'r')
      t = stats.variables['time'][:]
      z = stats.variables['z'][:]
      zh = stats.variables['zh'][:]
      
      thl_thl = stats.groups['thermo'].variables['thl'][:, :]
      sql = stats.groups['thermo'].variables['ql'][:,:]*10e5 	#to g m-3
      cft = stats.groups['thermo'].variables["ql_frac"][:,:]
      areat = stats.groups['default'].variables['area'][:, :]
#      cf = np.sum(areat[:,:]*cft[:,:], 0) / np.sum(areat[:,:], 0)
      
      #Condition for ql>0
      cft[cft==0] = np.nan
#      cf[cf==0] = np.nan
      sql[sql==0] = np.nan
      
      sql_mean = np.nanmean(sql[:,:],0)
      cft_mean = np.nanmean(cft[:,:],0)
      
      ql_thl.append(sql_mean)
      cf_thl.append(cft_mean)
      th_thl.append(thl_thl)
      print("Average cloud cover is {:.2f}%".format(np.nanmean(cft_mean)*100))

for j in all_files_qt:
      
      stats = nc.Dataset('qt/fireles.default.0000000_{}.nc'.format(j), 'r')
      t = stats.variables['time'][:]
      z = stats.variables['z'][:]
      zh = stats.variables['zh'][:]
      
      sql_qt = stats.groups['thermo'].variables['ql'][:,:]*10e5 	#to g m-3
      cft_qt = stats.groups['thermo'].variables["ql_frac"][:, :]
      areat_qt = stats.groups['default'].variables['area'][:, :]
#      scf_qt = np.sum(areat_qt[:,:]*cft_qt[:,:], 0) / np.sum(areat_qt[:,:], 0)
      
      #Condition for ql>0
      cft_qt[cft_qt==0] = np.nan
      scf_qt[scf_qt==0] = np.nan
      sql_qt[sql_qt==0] = np.nan
      
      sql_mean_qt = np.nanmean(sql_qt[:,:],0)
      cft_qt_mean = np.nanmean(cft_qt[:,:],0)
      
      ql_qt.append(sql_mean_qt)
      cf_qt.append(cft_qt_mean)
      print("Average cloud cover is {:.2f}%".format(np.nanmean(cft_qt_mean)*100))


#Lists to 2D arrays
ql_thl = np.concatenate(ql_thl)
ql_thl = np.reshape(ql_thl, (-1,len(z))).T
cf_thl = np.concatenate(cf_thl)
cf_thl = np.reshape(cf_thl, (-1,len(z))).T
th_thl = np.concatenate(th_thl)
th_thl = np.reshape(th_thl, (-1,len(z))).T

ql_qt = np.concatenate(ql_qt)
ql_qt = np.reshape(ql_qt, (-1,len(z))).T
cf_qt = np.concatenate(cf_qt)
cf_qt = np.reshape(cf_qt, (-1,len(z))).T


# Plots
plt.rc('font',**{'family':'serif','serif':['Palatino'], 'size':16})
plt.rc('text', usetex=True)


c10 = 'C1'
c20 = 'C2'
c30 = 'C3'
c40 = 'C9'
c50 = 'C5'
#
#plt.figure()
#plt.plot(th_thl[:,0], z, 'k:') #label = 't={:g}h'.format(t[n]/3600))
#plt.plot(th_thl[:, 1], z)
#plt.ylabel('z [m]')
#plt.xlabel('$/theta$ [K]')
#plt.title('Hourly Potential temperature profile')
#plt.legend()
#plt.show()


plt.figure(figsize = (12,8))
plt.plot(cf_thl[:,0], z, 'k:', label='base', lw=2)
plt.plot(cf_thl[:,1], z, label='thl10', lw=2, color=c10)
plt.plot(cf_thl[:,2], z, label='thl20', lw=2, color=c20)
plt.plot(cf_thl[:,3], z, label='thl30', lw=2, color=c30)
plt.plot(cf_thl[:,4], z, label='thl40', lw=2, color=c40)
plt.plot(cf_thl[:,5], z, label='thl50', lw=2, color=c50)
plt.xlabel('cloud fraction [-]')
plt.ylabel('z [m]')
plt.title('Mean vertical profiles of cloud fraction ')
plt.legend(prop={'size': 14})


plt.figure(figsize = (12,8))
plt.plot(cf_qt[:,0], z, 'k:', label='base', lw=2)
plt.plot(cf_qt[:,1], z, label='qt10', lw=2, color=c10)
plt.plot(cf_qt[:,2], z, label='qt20', lw=2, color=c20)
plt.plot(cf_qt[:,3], z, label='qt30', lw=2, color=c30)
plt.plot(cf_qt[:,4], z, label='qt40', lw=2, color=c40)
plt.plot(cf_qt[:,5], z, label='qt50', lw=2, color=c50)
plt.xlabel('cloud fraction [-]')
plt.ylabel('z [m]')
plt.title('Mean vertical profiles of cloud fraction ')
plt.legend(prop={'size': 14})


plt.figure(figsize = (12,8))
plt.plot(ql_thl[:,0], z, 'k:', label='base', lw=2)
plt.plot(ql_thl[:,1], z, label='thl10', lw=2, color=c10)
plt.plot(ql_thl[:,2], z, label='thl20', lw=2, color=c20)
plt.plot(ql_thl[:,3], z, label='thl30', lw=2, color=c30)
plt.plot(ql_thl[:,4], z, label='thl40', lw=2, color=c40)
plt.plot(ql_thl[:,5], z, label='thl50', lw=2, color=c50)
plt.xlabel('q$_l$ [g m$^{-3}$]')
plt.ylabel('z [m]')
plt.title('Mean vertical profiles of liquid water ')
plt.legend(prop={'size': 14})


plt.figure(figsize = (12,8))
plt.plot(ql_qt[:,0], z, 'k:', label='base', lw=2)
plt.plot(ql_qt[:,1], z, label='qt10', lw=2, color=c10)
plt.plot(ql_qt[:,2], z, label='qt20', lw=2, color=c20)
plt.plot(ql_qt[:,3], z, label='qt30', lw=2, color=c30)
plt.plot(ql_qt[:,4], z, label='qt40', lw=2, color=c40)
plt.plot(ql_qt[:,5], z, label='qt50', lw=2, color=c50)
plt.xlabel('q$_l$ [g m$^{-3}$]')
plt.ylabel('z [m]')
plt.title('Mean vertical profiles of liquid water ')
plt.legend(prop={'size': 14})
plt.show()


#      df[str(cs)] = sql_mean
#      df[str(cs)] = pd.DataFrame(sql_mean, index=z)
#      df[str(i)] = sql(i)
#      stats = sql(i)
#      print(stats)
