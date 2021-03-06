import netCDF4
import sys
import os
import glob
import numpy as np
sys.path.append('/usr/local/lib/python2.7/site-packages/')
import vtk 
from vtk.util.numpy_support import vtk_to_numpy

rhow=1.028
rhoi=0.918

ncfile = netCDF4.Dataset('/home/users/merino4i/output_files_MISMIP/output_files_MISMIP-/Ice2ra.nc','a')
xGL= ncfile.variables['xGL']
yGL= ncfile.variables['yGL']
iceThicknessGL = ncfile.variables['iceThicknessGL']
uBaseGL = ncfile.variables['uBaseGL']
vBaseGL = ncfile.variables['vBaseGL']
#vSurfaceGL = ncfile.variables['vSurfaceGL']
#uSurfaceGL = ncfile.variables['uSurfaceGL']
#uMeanGL = ncfile.variables['uMeanGL']
#vMeanGL = ncfile.variables['vMeanGL']

caseFirst='Test500m_Schoof_SSAStar_Repeated'
cases=['Test500m_Schoof_SSAStar_Repeated']
runs=['Ice2r9','Ice2a','Ice2a1','Ice2a2','Ice2a3','Ice2a4','Ice2a5','Ice2a6','Ice2a7','Ice2a8','Ice2a9',
        'Ice2a10','Ice2a11','Ice2a12','Ice2a13','Ice2a14','Ice2a15','Ice2a16','Ice2a17']
indexCases=0
time=[100.0,110.0,120.0,130.0,140.0,150.0,160.0,170.0,180.0,190.0,200.0,
        300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0]
#Integrales
Voldata=[]
VolSubdata=[]
AreaGdata=[]
XGL=[]
indexTime=0
for case in cases:
    for run in runs:
        indexFile=0
        if run==runs[0]: #First data comes from previous Run
            path='/home/users/merino4i/MISMIP+/'+caseFirst+'/'+run+'/'
            filesIce=glob.glob(path+'*.pvtu')
            filesIce.sort()
            for fileTest in filesIce:
                file1=fileTest
        else:
            
            path='/home/users/merino4i/MISMIP+/'+case+'/'+run+'/'
            filesIce=glob.glob(path+'*.pvtu')
            filesIce.sort()
            for fileTest in filesIce:
                file1=fileTest
                
        reader = vtk.vtkXMLPUnstructuredGridReader()
        print file1
        reader.SetFileName(file1)
        reader.Update()
        output=reader.GetOutput()
        Coords=vtk_to_numpy(output.GetPoints().GetData())
        PointData=output.GetPointData()
        numArrays=PointData.GetNumberOfArrays()   
        
        for i in np.arange(numArrays):
            if PointData.GetArrayName(i)=='ssavelocity':
                VarIndex=i
                break      
        Vel=vtk_to_numpy(PointData.GetArray(VarIndex))
        
        for i in np.arange(numArrays):
            if PointData.GetArrayName(i)=='groundedmask':
                VarIndex=i
                break      
        Gl=vtk_to_numpy(PointData.GetArray(VarIndex))

        for i in np.arange(numArrays):
            if PointData.GetArrayName(i)=='zb':
                VarIndex=i
                break      
        Zb=vtk_to_numpy(PointData.GetArray(VarIndex))

        for i in np.arange(numArrays):
            if PointData.GetArrayName(i)=='zs':
                VarIndex=i
                break      
        Zs=vtk_to_numpy(PointData.GetArray(VarIndex))
        
        for i in np.arange(numArrays):
            if PointData.GetArrayName(i)=='h':
                VarIndex=i
                break      
        H=vtk_to_numpy(PointData.GetArray(VarIndex))
        
        cellData=output.GetCellData()
        GeometryIDS=vtk_to_numpy(cellData.GetArray(0))
        indexGL=np.where(Gl==0.0)
        
        xGL[:,indexTime]=Coords[indexGL,0][0]
        yGL[:,indexTime]=Coords[indexGL,1][0]
        iceThicknessGL[:,indexTime]=H[indexGL]
        uBaseGL[:,indexTime]=Vel[indexGL,0][0]
        vBaseGL[:,indexTime]=Vel[indexGL,1][0]
        #uMeanGL[:,indexTime]=velGL[indexGL,0][0]
        #vMeanGL[:,indexTime]=velGL[indexGL,1][0]
        #uSurfaceGL[:,indexTime]=velGL[indexGL,0][0]
        #vSurfaceGL[:,indexTime]=velGL[indexGL,1][0]
        indexFile=indexFile+1
        indexTime=indexTime+1
            
ncfile.close()