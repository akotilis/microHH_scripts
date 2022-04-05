import matplotlib.pyplot as plt
import numpy as np
float_type = "f8"
np_dtype = np.float64 if float_type == "f8" else np.float32

itot = 128
jtot = 32
ktot = 32
xsize = 12800
ysize = 3200
zsize = 3200

endtime = 10800
dt = 1800
nt = int((endtime / dt)+1)

# s1 = np.zeros((nt, jtot, itot), dtype = np_dtype)
s1 = np.zeros([jtot, itot])     #this creates an array of jtot rows by itot columns, contains all 0s.
s1[11:22, 11:40] = 1.0 
# s1[22:28, 11:22] = 1.0

xlist = list(range(0,128))
ylist = list(range(0,32))

# X,Y = np.meshgrid(np.linspace(0, xsize, itot), np.linspace(0, ysize, jtot))
levels = 1

fig, ax = plt.subplots(figsize=(12,8))
contourf_ = ax.contourf(xlist,ylist, s1, levels=levels, vmin=0, vmax=1)
ax.set_xlabel("x domain")
ax.set_ylabel("y domain")
cbar = fig.colorbar(contourf_)


plt.show()