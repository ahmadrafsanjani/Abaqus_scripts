"""
==========================================================================
ABAQUS Python script
Project: Membrane and Bending Strain Energy for Shell Elements
Description: Deformation of a square frame with Periodic Boundary Conditions 
by Ahmad Rafsanjani (https://github.com/ahmadrafsanjani)
Last update: 6 June 2017 (tested on ABAQUS/CAE 6.12-1)
==========================================================================
"""

#   ------------------------------------------------------------------------
#   How to run this script?
#   Windows command line:   abaqus cae nogui=shell_strain_energy.py
#   ABAQUS CAE: File -> Run Script...-> Browse shell_strain_energy.py  

#   ------------------------------------------------------------------------
#   Initialize simulation

from abaqus import *
import testUtils
testUtils.setBackwardCompatibility()
from abaqusConstants import *
from caeModules import *
import sketch
import part
from regionToolset import*
import load
from numpy import*
from math import *
from shutil import *
import numpy as np
from odbAccess import *
import assembly

session.journalOptions.setValues(replayGeometry=COORDINATE,
                                 recoverGeometry=COORDINATE )


#   ------------------------------------------------------------------------
#   Project 

project='shell_energy'

try:
    os.remove(project+'.lck')
    os.remove(project+'.odb')
except OSError:
    pass

model      = project+'_model'
sketch     = project+'_sketch'
material   = project+'_mat'
section    = project+'_sec'
cae_file   = project+'.cae'

#   ------------------------------------------------------------------------
#   Parameters


E=1.0e9
nu=0.3

l=1.0
t=0.05


seed_size=l/20.0

#   ------------------------------------------------------------------------
#   Model 
m = mdb.Model(name=project)

# sketch
s= m.ConstrainedSketch(name='sketch', sheetSize=2*l)
s.rectangle(point1=(-l/2,-l/2),point2=(l/2,l/2))



sp= m.ConstrainedSketch(name='partition', sheetSize=2*l)
sp.Line(point1=(0,-l/2),point2=(0,l/2))

p=m.Part(name='part_1',dimensionality=THREE_D, type=DEFORMABLE_BODY)
p.BaseShell(sketch=s)
p.PartitionFaceBySketch(faces=p.faces,sketch=sp)

p.Set(name='e_L', edges=p.edges.findAt(((-l/2,0,0),),))
p.Set(name='e_R', edges=p.edges.findAt((( l/2,0,0),),))
p.Set(name='e_T', edges=p.edges.findAt((( 0,l/2,0),),))
p.Set(name='e_B', edges=p.edges.findAt((( 0,-l/2,0),),))
p.Set(name='e_M', edges=p.edges.findAt(((0,0,0),),))



#   ------------------------------------------------------------------------
#   materials

mat= m.Material(material)
mat.Elastic(table=((E, nu), ))
mat.Density(table=((1e-10, ), ))
#   ------------------------------------------------------------------------
#   section




m.HomogeneousShellSection(idealization=NO_IDEALIZATION,
                          integrationRule=SIMPSON,
                          material=material,
                          name=section,
                          numIntPts=5,
                          poissonDefinition=DEFAULT,
                          preIntegrate=OFF,
                          temperature=GRADIENT,
                          thickness=t,
                          thicknessField='',
                          thicknessModulus=None,
                          thicknessType=UNIFORM,
                          useDensity=OFF)

p.SectionAssignment(offset=0.0,
                    offsetField='',
                    offsetType=MIDDLE_SURFACE,
                    region=Region(faces=p.faces),
                    sectionName=section,
                    thicknessAssignment=FROM_SECTION)



#   ------------------------------------------------------------------------
#   Assembly

elem_Code=S6
elem_type = mesh.ElemType(elemCode=elem_Code)

A = m.rootAssembly

I=A.Instance(dependent=OFF, name='I_1', part=m.parts['part_1'])

A.seedPartInstance(size=seed_size, deviationFactor=0.1, minSizeFactor=0.1,
                   regions=(A.instances['I_1'],))
A.setMeshControls(regions=A.instances['I_1'].faces, elemShape=TRI)
A.setElementType(regions=(A.instances['I_1'].faces,), elemTypes=(elem_type,))
A.generateMesh(regions=(A.instances['I_1'],))


ni =A.instances['I_1'].nodes



#   ------------------------------------------------------------------------
#   Step

m.StaticStep(name='Loading',
             previous='Initial',
             description='Uniform loading',
             timePeriod=1.0,
             nlgeom=ON,
             initialInc=0.001,
             minInc=1e-8,
             maxInc=0.01,
             maxNumInc=100000,
             stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
             stabilizationMagnitude=2e-4,
             adaptiveDampingRatio=0.04)

m.FieldOutputRequest(name='F-Output-1',
                     createStepName='Loading',
                     timeInterval=0.01,
                     variables=('U','S','E','SF','SE','COORD','EVOL','STH'))

m.HistoryOutputRequest(name='H-Output-1',
                       createStepName='Loading',
                       timeInterval=0.01) 
                                             
#   ------------------------------------------------------------------------
#   Displacement BC


m.DisplacementBC(name='BC_L',
                 createStepName='Loading',
                 region=A.instances['I_1'].sets['e_L'],
                 u1=0,
                 u2=0,
                 u3=0)

m.DisplacementBC(name='BC_R',
                 createStepName='Loading',
                 region=A.instances['I_1'].sets['e_R'],
                 u1=0,
                 u2=0,
                 u3=0)

m.DisplacementBC(name='BC_M',
                 createStepName='Loading',
                 region=A.instances['I_1'].sets['e_M'],
                 u3=0.1*l)




#   ------------------------------------------------------------------------
#   Job

j=mdb.Job(name=project,
          model=project,
          numCpus=4,
          numDomains=4,
          multiprocessingMode=DEFAULT,
          description='Strain Energy of Shell Elements')
j.submit()
j.waitForCompletion()
mdb.saveAs(cae_file)


#   ------------------------------------------------------------------------
#   Post processing

odb = openOdb(project+'.odb')
frames = odb.steps['Loading'].frames


session.viewports['Viewport: 1'].setValues(displayedObject=odb)

ALLSE=session.XYDataFromHistory(odb=odb,
                                name='ALLSE',
                                outputVariableName='Strain energy: ALLSE for Whole Model',
                                steps=('Loading', 'Release', ), )

N_f=len(frames)

ENERGY=np.zeros((N_f,5))


for i in range(N_f):
    EVOL_i=frames[i].fieldOutputs['EVOL']
    STH_i=frames[i].fieldOutputs['STH']
    SF_i=frames[i].fieldOutputs['SF']
    SM_i=frames[i].fieldOutputs['SM']
    SE_i=frames[i].fieldOutputs['SE']
    SK_i=frames[i].fieldOutputs['SK']

    N_e=len(SE_i.values)
    BSE=zeros(N_e)
    MSE=zeros(N_e)

    for j in range(N_e):
        EVOL=EVOL_i.values[j].data     # element volume
        STH=STH_i.values[j].data       # element thickness
        AREA=EVOL/STH                  # element area
        SK1=SK_i.values[j].data[0]
        SK2=SK_i.values[j].data[1]
        SK3=SK_i.values[j].data[2]

        SM1=SM_i.values[j].data[0]
        SM2=SM_i.values[j].data[1]
        SM3=SM_i.values[j].data[2]


        SE1=SE_i.values[j].data[0]
        SE2=SE_i.values[j].data[1]
        SE3=SE_i.values[j].data[2]        
        SE4=SE_i.values[j].data[3]
        SE5=SE_i.values[j].data[4]
        SE6=SE_i.values[j].data[5]
        
        SF1=SF_i.values[j].data[0]
        SF2=SF_i.values[j].data[1]
        SF3=SF_i.values[j].data[2]        
        SF4=SF_i.values[j].data[3]
        SF5=SF_i.values[j].data[4]
        SF6=SF_i.values[j].data[5]

        BSE[j]=0.5*AREA*(SK1*SM1+SK2*SM2+SK3*SM3+SF4*SE4+SF5*SE5)
        MSE[j]=0.5*AREA*(SE1*SF1+SE2*SF2+SE3*SF3)

    ENERGY[i][0]=ALLSE[i][0]
    ENERGY[i][1]=ALLSE[i][1]         # Total Strain Energy (Abaqus)
    ENERGY[i][2]=sum(BSE)+sum(MSE)   # Total Strain Energy (integrated)
    ENERGY[i][3]=sum(BSE)            # Total Bending Strain Energy
    ENERGY[i][4]=sum(MSE)            # Total Membrane Strain Energy
        



np.savetxt('shell_strain_energy.csv', ENERGY, fmt='%f',delimiter=',')


