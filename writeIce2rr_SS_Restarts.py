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

caseFirst='Test500m_Schoof_SSAStar'
cases=['Test500m_Schoof_SSAStar_Repeated']
runs=['Ice2r9','Ice2r10','Ice2r11','Ice2r12','Ice2r13','Ice2r14','Ice2r15','Ice2r16','Ice2r17','Ice2r18','Ice2r19'
        'Ice2r20','Ice2r21','Ice2r22','Ice2r23','Ice2r24','Ice2r25','Ice2r26','Ice2r27']
indexCases=0
#Integrales
Voldata=[]
VolSubdata=[]
AreaGdata=[]
time=[100.0,110.0,120.0,130.0,140.0,150.0,160.0,170.0,180.0,190.0,200.0,
        300.0,400.0,500.0,600.0,700.0,800.0,900.0]
for case in cases:
    for run in runs:
        indexFile=0
        if run==runs[0]: #First data comes from previous Run
            path='/home/users/merino4i/MISMIP+/'+caseFirst+'/'+run+'/'
            filesIce=glob.glob(path+'*.pvtu')
            filesIce.sort()
            for file1 in filesIce:
                file1=file1
        else:
            path='/home/users/merino4i/MISMIP+/'+case+'/'+run+'/'
            filesIce=glob.glob(path+'*.pvtu')
            filesIce.sort()
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

                
ncfile = netCDF4.Dataset('/home/users/merino4i/output_files_MISMIP/output_files_MISMIP-/Ice2rr_SS.nc','a')
iceVolume= ncfile.variables['iceVolume']
iceVAF= ncfile.variables['iceVAF']
groundedArea = ncfile.variables['groundedArea']
timeVar = ncfile.variables['time']
iceVolume[:]=np.array(Voldata)
groundedArea[:]=np.array(AreaGdata)
iceVAF[:]=np.array(Voldata)-np.array(VolSubdata)*(rhow/rhoi)
timeVar[:]=np.array(time)
ncfile.close()

