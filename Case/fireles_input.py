import matplotlib.pyplot as pl
import numpy as np
import netCDF4 as nc

float_type = "f8"
# float_type = "f4"

np_dtype = np.float64 if float_type == "f8" else np.float32

# Get number of vertical levels and size from .ini file
with open('fireles.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])
        if(line.split('=')[0]=='itot'):
            itot = int(line.split('=')[1])
        if(line.split('=')[0]=='jtot'):
            jtot = int(line.split('=')[1])

dz = zsize / kmax

dthetadz = 0.003

# set the height
z  = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
u  = np.ones(np.size(z))*5
v  = np.zeros(np.size(z))
thl = np.zeros(np.size(z))
qt    = np.zeros(z.size)

s1 = np.zeros(z.size)

# linearly stratified profile
for k in range(kmax):
    thl  [k] = 300. + dthetadz*z[k]
    if (z[k] <= 1000):
        s1[k] = 0.0
        qt[k] = 7.0
    else:
        s1[k] = 0.0
        qt[k] = 7.0 - (z[k]-1000.0)*0.001
qt *= 1e-3  # kg/kg

nc_file = nc.Dataset("fireles_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", kmax)
nc_z  = nc_file.createVariable("z" , float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_u = nc_group_init.createVariable("u", float_type, ("z"))
nc_v = nc_group_init.createVariable("v", float_type, ("z"))
nc_thl = nc_group_init.createVariable("thl", float_type, ("z"))
nc_qt = nc_group_init.createVariable("qt", float_type, ("z"))
nc_s1 = nc_group_init.createVariable("s1", float_type, ("z"))

nc_z [:] = z [:]
nc_u [:] = u [:]
nc_v [:] = v [:]
nc_thl[:] = thl[:]
nc_qt[:] = qt[:]
nc_s1[:] = s1[:]

nc_file.close()

# Create time varying surface fields.
endtime = 10800*2.5
dt = 1800
nt = int((endtime / dt)+1)

thl_sbot = np.zeros((nt, jtot, itot), dtype=np_dtype)
s1_sbot = np.zeros((nt, jtot, itot), dtype=np_dtype)
qt_sbot = np.zeros((nt, jtot, itot), dtype=np_dtype)

thl_sbot[:] = 0.2
#thl_sbot[:, 11:22, 11:22] = 0.4

s1_sbot[:] = 0.0
#s1_sbot[:, 11:22, 11:22] = 1.0

qt_sbot[:] = 1e-4
#qt_sbot[:, 11:22, 11:22] = 5e-3

# Write as binary input files for MicroHH
for t in range(nt):
    thl_sbot[t,:].tofile('thl_flux.{0:07d}'.format(t*dt))
    s1_sbot[t,:].tofile('s1_flux.{0:07d}'.format(t*dt))
    qt_sbot[t,:].tofile('qt_flux.{0:07d}'.format(t*dt))
