"""
==========================================================================
ABAQUS Python script
Project: Square PBC
Description: Deformation of a square frame with Periodic Boundary Conditions 

by Ahmad Rafsanjani (https://github.com/ahmadrafsanjani)

Last update: 26 March 2017 (tested on ABAQUS/CAE 6.12-1)
==========================================================================
"""

#   ------------------------------------------------------------------------
#   How to run this script?
#   Windows command line:   abaqus cae nogui=square_pbc.py
#   ABAQUS CAE: File -> Run Script...-> Browse square_pbc.py  

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
import numpy
from odbAccess import *
import assembly

session.journalOptions.setValues(replayGeometry=COORDINATE,
                                 recoverGeometry=COORDINATE )


#   ------------------------------------------------------------------------
#   Project 

project='square_pbc'

try:
    os.remove(project+'.lck')
    os.remove(project+'.odb')
except OSError:
    pass

model         = project+'_model'
sketch_name   = project+'_sketch'
material_name = project+'_mat'
section_name  = project+'_sec'
cae_file      = project+'.cae'

#   ------------------------------------------------------------------------
#   Simulation Parameters

E=1.0e9   # elastic Young's modulus
nu=0.3    # Poisson's ratio"

l=1.0     # square length
t=0.1*l   # frame thickness
ts=0.01*l # plae thickness

seed_size=t/4    # mesh size


X0=l/2
Y0=l/2

xmin=-X0
xmax= X0
ymin=-Y0
ymax= Y0


time_inc=0.1  # time increment for writing output results


#   ------------------------------------------------------------------------
#   Applied deformations

# Example_1: uniform shear 
X_u1=0
X_u2=0.1*l
Y_u1=0.1*l
Y_u2=0

# Example_2: uniaxial extension (along X)
##X_u1=0.1*l
##X_u2=0
##Y_u1=0
##Y_u2=UNSET

# Example_3: uniform Bi-axial extension
##X_u1=0.1*l
##X_u2=0
##Y_u1=0
##Y_u2=0.1*l

#   ------------------------------------------------------------------------
#   Functions

def sort_node_list(node_list,axis):
    # Sorting a node list along a selected axis
    node_list.sort(key=lambda node_list: node_list.coordinates[axis], reverse=False)



def ncoord(n,i):
    # Returning the coordinate of a node in i axis 
    return round(n.coordinates[i],4)


def nappend(nset,n,i):
    nset.append(n.sequenceFromLabels(labels=(n[i].label,)))
    # append node to a node set

def periodic_bc(model,node_set1,node_set2,master,dof):
    # 1-to-1 periodic BC between two node sets for selected DOFs
    # Equal mesh on opposite boundaries is required

    num_nodes=len(node_set1)

    for i in range(num_nodes):
        n1=node_set1[i]
        n2=node_set2[i]


        n1_mask = ni.sequenceFromLabels(labels=(n1.label,))
        n2_mask = ni.sequenceFromLabels(labels=(n2.label,))
        
        n1_region = A.Set(nodes=n1_mask, name='n1_'+str(n1.label))
        n2_region = A.Set(nodes=n2_mask, name='n2_'+str(n2.label))


        for dof_i in dof:
            model.Equation(name='PBC_'+str(dof_i)+'_'+str(n1.label),
                           terms=(( 1.0 , 'n1_'+str(n1.label), dof_i),
                                  (-1.0 , 'n2_'+str(n2.label), dof_i),
                                  (-1.0, master, dof_i)))

    
#   ------------------------------------------------------------------------
#   Model & Sketches & Parts & Partitions

m = mdb.Model(name=project)

s = m.ConstrainedSketch(name=sketch_name, sheetSize=2*(X0+Y0))
s.rectangle(point1=(-l/2,-l/2),point2=(l/2,l/2))
s.rectangle(point1=(-l/2+t,-l/2+t),point2=(l/2-t,l/2-t))

sp = m.ConstrainedSketch(name=sketch_name+'partition', sheetSize=2*(X0+Y0))
sp.Line(point1=(-l/2, -l/2+t) ,point2=(l/2, -l/2+t))
sp.Line(point1=(-l/2,  l/2-t) ,point2=(l/2,  l/2-t))
sp.Line(point1=(-l/2+t, -l/2) ,point2=(-l/2+t, l/2))
sp.Line(point1=( l/2-t, -l/2) ,point2=( l/2-t, l/2))

p = m.Part(name=model, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
p.BaseShell(sketch=s)
p.PartitionFaceBySketch(faces=p.faces,sketch=sp)

px = m.Part(name='Ref_X', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
px.ReferencePoint(point=(X0, -Y0, 0))

py = m.Part(name='Ref_Y', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
py.ReferencePoint(point=(-X0, Y0, 0))

#   ------------------------------------------------------------------------
#   Material

mat= m.Material(material_name)
mat.Elastic(table=((E, nu), ))  # elastic material

#   ------------------------------------------------------------------------
#   Section

m.HomogeneousSolidSection(material=material_name, # material name
                          name=section_name,      # section name
                          thickness=ts)           # shell thickness

p.SectionAssignment(offset=0.0,
                    offsetField='',
                    offsetType=MIDDLE_SURFACE,
                    region=Region(faces=p.faces),
                    sectionName=section_name,
                    thicknessAssignment=FROM_SECTION)



#   ------------------------------------------------------------------------
#   Assembly

A = m.rootAssembly
I=A.Instance(dependent=ON, name=project, part=p)
A.makeIndependent(instances=(I, ))
A.Instance(name='Ref_X-1', part=px, dependent=ON)
A.Instance(name='Ref_Y-1', part=py, dependent=ON)
A.regenerate()

#   ------------------------------------------------------------------------
#   Mesh

elem_Code=CPE4
elem_type = mesh.ElemType(elemCode=elem_Code)

A.setElementType(regions=Region(faces=I.faces), elemTypes=(elem_type,))
A.seedEdgeBySize(edges=I.edges, constraint=FIXED, size=seed_size)
A.setMeshControls(regions=I.faces,elemShape=QUAD,technique=STRUCTURED)
A.generateMesh(regions=(I,))




#   ------------------------------------------------------------------------
#   Geometrical features

REFX=A.Set(name='RefX',referencePoints=(A.instances['Ref_X-1'].referencePoints[1],))
REFY=A.Set(name='RefY',referencePoints=(A.instances['Ref_Y-1'].referencePoints[1],))

fi = A.instances[project].faces
ei = A.instances[project].edges
ni = A.instances[project].nodes
vi = A.instances[project].vertices
eli = A.instances[project].elements

e_top=[]
e_bot=[]
e_left=[]
e_right=[]

for ej in ei:
    (xj,yj,zj)=ej.pointOn[0]

    if round(xj-xmin,4)==0:
        e_left.append(ei.findAt((ej.pointOn[0], ), ))
    elif round(xj-xmax,4)==0:
        e_right.append(ei.findAt((ej.pointOn[0], ), ))
    elif round(ymin-yj,4)==0:
        e_bot.append(ei.findAt((ej.pointOn[0], ), ))
    elif round(ymax-yj,4)==0:
        e_top.append(ei.findAt((ej.pointOn[0], ), ))


e_t=A.Set(edges=e_top,   name='e_top')
e_b=A.Set(edges=e_bot,   name='e_bot')
e_l=A.Set(edges=e_left,  name='e_left')
e_r=A.Set(edges=e_right, name='e_right')

n_top    = list(A.sets['e_top'].nodes)
n_bot    = list(A.sets['e_bot'].nodes)
n_right  = list(A.sets['e_right'].nodes)
n_left   = list(A.sets['e_left'].nodes)

sort_node_list(n_top,0)
sort_node_list(n_bot,0)
sort_node_list(n_right,1)
sort_node_list(n_left,1)

N_O=A.Set(nodes=ni.sequenceFromLabels(labels=(n_bot[0].label,)), name='Origin')

#   ------------------------------------------------------------------------
#   Periodic Boundary Conditions (PBC)

periodic_bc(m,n_top,  n_bot, 'RefY',[1,2])
periodic_bc(m,n_right,n_left,'RefX',[1,2])

#   ------------------------------------------------------------------------
#   Step

m.StaticStep(name='Loading',
             previous='Initial',
             description='Periodic deformation',
             timePeriod=1.0)


m.FieldOutputRequest(name='F-Output-1',
                     createStepName='Loading',
                     timeInterval=time_inc,
                     variables=('U','S','E','IVOL','RF','COORD'))

m.HistoryOutputRequest(name='H-Output-1',
                       createStepName='Loading',
                       timeInterval=time_inc)
                                          
#   ------------------------------------------------------------------------
#   Displacement BC

m.DisplacementBC(name='BC_O',
                 createStepName='Loading',
                 region=N_O,
                 u1=0,
                 u2=0)

m.DisplacementBC(name='BC_REFX',
                 createStepName='Loading',
                 region=REFX,
                 u1=X_u1,
                 u2=X_u2)

m.DisplacementBC(name='BC_REFY',
                 createStepName='Loading',
                 region=REFY,
                 u1=Y_u1,
                 u2=Y_u2)

#   ------------------------------------------------------------------------
#   Job

j=mdb.Job(name=project,
          model=project,
          numCpus=1,
          numDomains=1,
          multiprocessingMode=DEFAULT,
          description='Deformation of a square frame with Periodic Boundary Conditions')
j.submit()
j.waitForCompletion()
mdb.saveAs(cae_file)

