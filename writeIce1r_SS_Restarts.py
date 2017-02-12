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

caseFirst='Conv500m_Schoof_SSAStar'
cases=['Test500m_Schoof_SSAStar']
runs=['Run7','Ice1r','Ice1r1','Ice1r2','Ice1r3','Ice1r4','Ice1r5','Ice1r6','Ice1r7','Ice1r8','Ice1r9']
indexCases=0
#Integrales
Voldata=[]
VolSubdata=[]
AreaGdata=[]
time=[0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100]
for case in cases:
    for run in runs:
        indexFile=0
        if run==runs[0]: #First data comes from previous Run
            path='/Users/imerino/Documents/These/MISMIP+/Occigen/'+caseFirst+'/'+run+'/'
            filesIce=glob.glob(path+'*.pvtu')
            for file1 in filesIce:
                file1=file1
        else:
            path='/Users/imerino/Documents/These/MISMIP+/Occigen/'+case+'/'+run+'/'
            filesIce=glob.glob(path+'*.pvtu')
            for file1 in filesIce:
                file1=file1
                
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
        pointData=output.GetPointData()
        GeometryIDS=vtk_to_numpy(cellData.GetArray(0))
        indexGEO=np.where(GeometryIDS==1)
        
        #Integrales
        #Volume
        listPoints=set()
        VolInteg=0
        for i in indexGEO[0]:
        #La primera celda de todas
            celda1=output.GetCell(i)
            ids=celda1.GetPointIds()
            if ids.GetNumberOfIds()==3:
                mean=0
                for j in np.arange(3):
                    mean=mean+H[ids.GetId(j)]/3
            VolLocal=mean*celda1.ComputeArea()
            VolInteg=VolInteg+VolLocal
            
        Voldata.append(VolInteg)
        
        #Volume_Submerged
        listPoints=set()
        VolInteg=0
        for i in indexGEO[0]:
        #La primera celda de todas
            celda1=output.GetCell(i)
            ids=celda1.GetPointIds()
            if ids.GetNumberOfIds()==3:
                mean=0
                for j in np.arange(3):
                    mean=mean+Zb[ids.GetId(j)]/3
            VolLocal=mean*celda1.ComputeArea()
            if VolLocal>0:
                VolLocal=0
            else:
                VolLocal=-VolLocal
            VolInteg=VolInteg+VolLocal
            
        VolSubdata.append(VolInteg)
        
        #Grounded Area
        listPoints=set()
        AreaInteg=0
        for i in indexGEO[0]:
        #La primera celda de todas
            celda1=output.GetCell(i)
            ids=celda1.GetPointIds()
            if ids.GetNumberOfIds()==3:
                mean=0
                for j in np.arange(3):
                    mean=mean+Gl[ids.GetId(j)]/3
            if mean>0:
                AreaLocal=celda1.ComputeArea()
            else:
                AreaLocal=0.             
            AreaInteg=AreaInteg+AreaLocal
            
        AreaGdata.append(AreaInteg)
        
        indexFile=indexFile+1

                
ncfile = netCDF4.Dataset('/Users/imerino/Documents/These/MISMIP+/Outputs/Write_Output_From_VTK/Ice1r.nc','a')
iceVolume= ncfile.variables['iceVolume']
iceVAF= ncfile.variables['iceVAF']
groundedArea = ncfile.variables['groundedArea']
timeVar = ncfile.variables['time']
iceVolume[:]=np.array(Voldata)
groundedArea[:]=np.array(AreaGdata)
iceVAF[:]=np.array(Voldata)-np.array(VolSubdata)*(rhow/rhoi)
timeVar[:]=np.array(time)
ncfile.close()