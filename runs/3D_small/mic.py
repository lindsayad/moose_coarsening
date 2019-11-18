import numpy as np
from netCDF4 import Dataset
import sys
import os

N=64
mic=np.zeros([N,N,N])
xx=np.arange(N)
x, y, z=np.meshgrid(xx, xx, xx)
mic[(z-N/2)*(z-N/2)+(x-3*N/4)*(x-3*N/4)+(y-3*N/4)*(y-3*N/4)<=50]=1
mic[(z-N/2)*(z-N/2)+(x-N/4)*(x-N/4)+(y-3*N/4)*(y-3*N/4)<=100]=1
mic[y<=N/3]=1


print(mic.shape)
ncdf=Dataset('mic.nc','w',type='NETCDF4')
ncdf.history="initial structure from DREAM3D"
ncdf.createDimension('x',mic.shape[2])
ncdf.createDimension('y',mic.shape[1])
ncdf.createDimension('z',mic.shape[0])
g=ncdf.createVariable('microstr','i',('z','y','x'))
dx=ncdf.createVariable('dx','d')
dy=ncdf.createVariable('dy','d')
dz=ncdf.createVariable('dz','d')
g[...]=mic
dx[...]=1.0
dy[...]=1.0
dz[...]=1.0
ncdf.close()

