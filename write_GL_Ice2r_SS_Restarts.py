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

ncfile = netCDF4.Dataset('/home/users/merino4i/output_files_MISMIP/output_files_MISMIP-/Ice1r_SS.nc','a')
xGL= ncfile.variables['xGL']
yGL= ncfile.variables['yGL']
iceThicknessGL = ncfile.variables['iceThicknessGL']
uBaseGL = ncfile.variables['uBaseGL']
vBaseGL = ncfile.variables['vBaseGL']
#vSurfaceGL = ncfile.variables['vSurfaceGL']
#uSurfaceGL = ncfile.variables['uSurfaceGL']
#uMeanGL = ncfile.variables['uMeanGL']
#vMeanGL = ncfile.variables['vMeanGL']

caseFirst='Conv500m_Schoof_SSAStar'
cases=['Test500m_Schoof_SSAStar_Repeated']
runs=['Run7','Ice2r','Ice2r1','Ice2r2','Ice2r3','Ice2r4','Ice2r5','Ice2r6','Ice2r7','Ice2r8','Ice2r9']
indexCases=0
time=[0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100]
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