#First, we import all libraries we need for this project
#Because, we will use array, we have to import numpy, plt is for ploting our array after processing
import matplotlib.pyplot as plt
from math import *
from openmc import *
from numpy import*

#####################################
#Définition of modelisation settings#
#####################################
#We define there all properties 
#Thickness of materials 
e_plasma = 50
e_wall1 = 3 
e_vessel = 2 
e_shield = 10
e_wall2 = 1 
e_tfcoil = 20 
#Target size
H_target = 10
e_target = 3
#Plasma density
plasma_density = 0.01 #g/cm3
#Set the number of batch and particules/batch
N_batch = 50
N_particules = 1000
#mesh for plotting : precision of all
mesh_size = 1
angle_inter_plane = 20 #in degrees

#Need to be calculate
y_1 = tan(pi*angle_inter_plane/180)
size = int(300/mesh_size)

##################
#Set-up materials#
##################
tungsten=openmc.Material(name="Tungsten")
tungsten.add_element("W",1)
tungsten.set_density("g/cm3",19.3) 

ss316=openmc.Material(name="ss316")
ss316.add_element("Si",4.998E-03,"wo")
ss316.add_element("Ti",1.5E-05,"wo")
ss316.add_element("Cr",1.75E-01,"wo")
ss316.add_element("Cu",3.0E-03,"wo")
ss316.add_element("Mn",1.8E-02,"wo")
ss316.add_element("Fe",2.4814,"wo")
ss316.add_element("Ni",1.1944E-01,"wo")
ss316.add_element("Mo",2.5E-02,"wo")
ss316.add_element("C",0.0002,"wo")
ss316.set_density("g/cm3",8.0) 

B_carbide=openmc.Material(name="Boron carbide")
B_carbide.add_element("B",0.8)
B_carbide.add_element("C",0.2)
B_carbide.set_density("g/cm3",2.45)

incoloy=openmc.Material(name="incoloy 908")
incoloy=openmc.Material(4,"incoloy")
incoloy.add_element("Ni",48.0)
incoloy.add_element("Fe",40.0)
incoloy.add_element("Cr",4.0)
incoloy.add_element("Cb",3.0)
incoloy.add_element("Ti",1.5)
incoloy.add_element("Mn",1.0)
incoloy.add_element("Al",1.0)
incoloy.add_element("Co",0.5)
incoloy.add_element("Cu",0.5)
incoloy.add_element("Si",0.5)
incoloy.add_element("C",0.03)
incoloy.add_element("P",0.015)
incoloy.add_element("B",0.012)
incoloy.add_element("S",0.005)
incoloy.set_density("g/cm3",8.0)

coil_mat=openmc.Material(5,"coil")
coil_mat.add_element("Cu",98.8)
coil_mat.add_element("Zr",0.20)
coil_mat.add_element("Cr",0.9)
coil_mat.set_density("g/cm3",8.9)

plasma_H=openmc.Material(name="plasma")
plasma_H.add_nuclide("H2",0.1)
plasma_H.add_nuclide("H3",0.1)
plasma_H.set_density("g/cm3",plasma_density)

#export at XML format
mats = openmc.Materials([tungsten,ss316,B_carbide,incoloy,coil_mat,plasma_H])
isinstance(mats, list)
mats.export_to_xml()

########################################
#Create geometry and associate to cells#
#######################################
#create all cells
plasma=openmc.Cell(1,"Plasma")
wall1=openmc.Cell(2,"Wall1")
vessel=openmc.Cell(3,"Vessel")
shield=openmc.Cell(4,"Shield")
wall2=openmc.Cell(5,"wall2")
tfcoil=openmc.Cell(6,"tfcoil")
target=openmc.Cell(7,"target")
n_phage=openmc.Cell(8,"n_phage")

#create cylinder
R1 = e_tfcoil - e_target 
R2 = R1 + e_target 
R3 = R2 + e_wall2
R4 = R3 + e_shield
R5 = R4 + e_vessel
R6 = R5 + e_wall1
R7 = R6 + e_plasma
R8 = R7 + e_wall1
R9 = R8 + e_vessel
R10 = R9 + e_shield
R11 = R10 + e_wall2
R12 = R11 + e_tfcoil 

#create border
min_x = openmc.XPlane(x0=0,boundary_type='reflective')
max_x = openmc.XPlane(x0=R9+100,boundary_type='reflective')

#create border of the 20 degre wedge
x_border = openmc.YPlane(y0=0,boundary_type='reflective')
border_angle = openmc.Plane(A = y_1 ,B = 1,boundary_type='reflective')

tfcoil_plane = openmc.XPlane(x0=R1)
target_plane = openmc.XPlane(x0=R2)

cylinder_3 = openmc.ZCylinder(x0=0.0, y0=0.0, R=R3)
cylinder_4 = openmc.ZCylinder(x0=0.0, y0=0.0, R=R4)
cylinder_5 = openmc.ZCylinder(x0=0.0, y0=0.0, R=R5)
cylinder_6 = openmc.ZCylinder(x0=0.0, y0=0.0, R=R6)
cylinder_7 = openmc.ZCylinder(x0=0.0, y0=0.0, R=R7)
cylinder_8 = openmc.ZCylinder(x0=0.0, y0=0.0, R=R8)
cylinder_9 = openmc.ZCylinder(x0=0.0, y0=0.0, R=R9)
cylinder_10 = openmc.ZCylinder(x0=0.0, y0=0.0, R=R10)
cylinder_11 = openmc.ZCylinder(x0=0.0, y0=0.0, R=R11)
cylinder_12 = openmc.ZCylinder(x0=0.0, y0=0.0, R=R12)

tfcoil_1 = -x_border & +border_angle & -tfcoil_plane & +min_x
target_r = -x_border & +border_angle & +tfcoil_plane & -target_plane & +min_x
wall2_1 = -x_border & +border_angle & +target_plane & -cylinder_3 & +min_x
shield_1 = -x_border & +border_angle & +cylinder_3 & -cylinder_4 & +min_x
vessel_1 = -x_border & +border_angle & +cylinder_4 & -cylinder_5 & +min_x
wall1_1 =-x_border & +border_angle & +cylinder_5 & -cylinder_6 & +min_x
plasma_r = -x_border & +border_angle & +cylinder_6 & -cylinder_7 & +min_x
wall1_2 = -x_border & +border_angle & +cylinder_7 & -cylinder_8 & +min_x
vessel_2 = -x_border & +border_angle & +cylinder_8 & -cylinder_9 & +min_x
shield_2 = -x_border & +border_angle & +cylinder_9 & -cylinder_10 & +min_x
wall2_2 = -x_border & +border_angle & +cylinder_10 & -cylinder_11 & +min_x
tfcoil_2 = -x_border & +border_angle & +cylinder_11 & -cylinder_12 & +min_x

plasma.region = plasma_r
wall1.region = wall1_1 | wall1_2
vessel.region = vessel_1 | vessel_2
shield.region = shield_1 | shield_2
wall2.region = wall2_1 | wall2_2
tfcoil.region = tfcoil_1 | tfcoil_2
target.region = target_r
n_phage.region = -x_border & +border_angle & +cylinder_12 & -max_x & +min_x

#fill cells
plasma.fill = plasma_H
wall1.fill = tungsten
vessel.fill = ss316
shield.fill = B_carbide
wall2.fill = ss316
tfcoil.fill = incoloy
target.fill = incoloy
n_phage.fill = tungsten

#Add to the universe
universe = openmc.Universe()
universe.add_cell(plasma)
universe.add_cell(wall1)
universe.add_cell(vessel)
universe.add_cell(shield)
universe.add_cell(wall2)
universe.add_cell(tfcoil)
universe.add_cell(target)
universe.add_cell(n_phage)

#export XML file
geom = openmc.Geometry(universe)
geom.export_to_xml()

#####################
#simulation settings#
#####################

#Create the source and XML file
src = openmc.Source()
src.angle = openmc.stats.Isotropic()
src.energy = openmc.stats.Discrete([14.1e6], [1.0])
settings = openmc.Settings()
settings.run_mode = "fixed source"


settings.source = src
settings.batches = N_batch
settings.inactive = 0
settings.particles = N_particules
settings.export_to_xml()


################
#Geometry Plot#
###############
B = zeros((300,300))
for i in range (0,300,1):
    for j in range (-150,150,1):
        if (i,j,0) in plasma:
            B[j+150][i+10]=1
        elif (i,j,0) in wall1:
            B[j+150][i+10]=2
        elif (i,j,0) in vessel:
            B[j+150][i+10]=3
        elif (i,j,0) in shield:
            B[j+150][i+10]=4
        elif (i,j,0) in wall2:
            B[j+150][i+10]=5
        elif (i,j,0) in tfcoil:
            B[j+150][i+10]=6
        elif (i,j,0) in target:
            B[j+150][i+10]=7
        elif (i,j,0) in n_phage:
            B[j+150][i+10]=9

            
plt.imshow(B)
plt.savefig("geometry.pdf")

#####################
#Filters and tallies#
#####################

tallies_file=openmc.Tallies()

# Instantiate a tally mesh
mesh = openmc.Mesh(mesh_id=1)
mesh.type = 'regular'
mesh.dimension = [size,size]
mesh.lower_left = [0, -100.0]
mesh.width = [mesh_size,mesh_size]
mesh_filter = openmc.MeshFilter(mesh)

# Instantiate flux Tally in moderator and fuel
tally = openmc.Tally(name='flux')
tally.filters = [mesh_filter]
tally.scores = ['flux','heating']
tallies_file.append(tally)

# Instantiate absorption in all cells except shield
cells_abs = openmc.Tally(name='cell_absorption')
cells_abs.filters = [openmc.CellFilter([plasma,wall1,vessel,shield,wall2,tfcoil,target,n_phage])]
cells_abs.scores = ['heating']
tallies_file.append(cells_abs)

# Export to "tallies.xml"
tallies_file.export_to_xml()

#Run
openmc.run()

###################
# POST PROCESSING #
###################

# We do not know how many batches were needed to satisfy the 
# tally trigger(s), so find the statepoint file(s)
statepoints = glob.glob('statepoint.*.h5')

# Load the last statepoint file
sp = openmc.StatePoint(statepoints[-1])

#Get the flux in all fusion cylinder
flux = sp.get_tally(name='flux')
flux.get_pandas_dataframe()

#Get absorption in the cells
absorp = sp.get_tally(name='cell_absorption')
heat_deposition = absorp.get_pandas_dataframe()

#For ploting
flux_f = flux.get_slice(scores=['flux'])
flux_a = flux.get_slice(scores=['heating'])

flux_f.std_dev.shape = (size, size)
flux_f.mean.shape = (size, size)
flux_a.std_dev.shape = (size, size)
flux_a.mean.shape = (size, size)

# Determine relative error
relative_error = np.zeros_like(flux_a.std_dev)
for i in range (0,size-1):
    for j in range (0,size-1):
        if flux_a.mean[i][j]>0:
            relative_error[i][j] = flux_a.std_dev[i][j] / flux_a.mean[i][j]        
b=0            
for i in range (0,300,1):
    for j in range (-150,150,1):
        if (i,j,0) in target:
            if relative_error[j+100][i]>b:
                b=relative_error[j+100][i]
print("max error is",b)

print(heat_deposition)

for i in range (0,300,1):
    for j in range (-150,150,1):
        if (i,j,0) in target:
            print("Heating (W/sources particules) in this part is :",flux_a.mean[j+100][i]*1.6*10**-19)

#for log ploting
for i in range (0,size-1):
    for j in range (0,size-1):
        if flux_a.mean[i][j]!=0 :
            flux_a.mean[i][j]=log(flux_a.mean[i][j])   

plt.imshow(flux_a.mean)
plt.savefig("heating.pdf") 