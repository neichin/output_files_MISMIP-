
import netCDF4
import sys
import os
import glob
import numpy as np



path='./'
pathNewFile= path+'temp.nc'
ncfile = netCDF4.Dataset(pathNewFile,'w',format='NETCDF4')

ncfile.createDimension('nPointGL',None)
ncfile.createDimension('nTime',11)

ncfile.createVariable('time','f',('nTime'))
ncfile.createVariable('iceVolume','f',('nTime'))
ncfile.createVariable('iceVAF','f',('nTime'))
ncfile.createVariable('groundedArea','f',('nTime'))
ncfile.createVariable('xGL','f',('nPointGL','nTime'))
ncfile.createVariable('yGL','f',('nPointGL','nTime'))
ncfile.createVariable('iceThicknessGL','f',('nPointGL','nTime'))
ncfile.createVariable('uBaseGL','f',('nPointGL','nTime'))
ncfile.createVariable('vBaseGL','f',('nPointGL','nTime'))
ncfile.createVariable('vSurfaceGL','f',('nPointGL','nTime'))
ncfile.createVariable('uSurfaceGL','f',('nPointGL','nTime'))
ncfile.createVariable('uMeanGL','f',('nPointGL','nTime'))
ncfile.createVariable('vMeanGL','f',('nPointGL','nTime'))
ncfile.close()
