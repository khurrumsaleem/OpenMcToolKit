"""This notebook show how we can create a model for fusion testing. We are able to change all 
shape and thickness really easily. This model wad created aiming at bringing some heating deposition
theorical mesure on SPARC reactor"""

"We have to import this module in order to make ploting, calculation and mesures"

import matplotlib.pyplot as plt
from math import *
from openmc import *
from numpy import*

"We define there all thickness and others input for the modelisation"

#####################################
#Définition of modelisation settings#
#####################################
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
N_batch = 1000
N_particules = 100000
#mesh for plotting : precision of all
mesh_size = 1
angle_inter_plane = 20 #in degrees

#Need to be calculate
y_1 = tan(pi*angle_inter_plane/180)

"We define there all materials we need, with the real proportion in SPARC"

##################
#Set-up materials#
##################
tungsten=openmc.Material(name="Tungsten")
tungsten.add_element("W",1.0)
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

" We can now create the file use by openmc"
#export at XML format
mats = openmc.Materials([tungsten,ss316,B_carbide,incoloy,coil_mat,plasma_H])
isinstance(mats, list)
mats.export_to_xml()

" Now let's move on geometry work. We want to build a 3D model, all will be explain step by step"
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
out=openmc.Cell(8,"out")

#create cylinder
R1 = e_plasma
R2 = e_plasma+e_wall1
R3 = e_plasma+e_wall1+e_vessel
R4 = e_plasma+e_wall1+e_vessel+e_shield
R5 = e_plasma+e_wall1+e_vessel+e_shield+e_wall2
R6 = e_plasma+e_wall1+e_vessel+e_shield+e_wall2+e_tfcoil

D1 = R2-R1
D2 = R3-R1 
D3 = R4-R1 
D4 = R5-R1 
D5 = R6-R1

#create border
border_x_1 = D5
border_x_2 = R6
min_x=openmc.XPlane(x0=-border_x_1,boundary_type='reflective')
max_x=openmc.XPlane(x0=border_x_2,boundary_type='reflective')
min_y=openmc.YPlane(y0=-200,boundary_type='reflective')
max_y=openmc.YPlane(y0=200,boundary_type='reflective')

#create border of the 20 degre wedge
border_angle_1 = openmc.Plane(A = y_1/2 ,C = 1,D = -D5*y_1/2-1,boundary_type='reflective')
border_angle_2 = openmc.Plane(A = -y_1/2 ,C = 1,D = D5*y_1/2+1,boundary_type='reflective')

"Because we want to mesure the heat deposition on a tiny area of the solenoid, we define a new cell call target"
trait_target1 = openmc.XPlane(x0=-D4-e_target)
trait_target2 = openmc.YPlane(y0=H_target/2)
trait_target3 = openmc.YPlane(y0=-H_target/2)

plasma_cylinder = openmc.ZCylinder(x0=0.0, y0=0.0, R=R1)
wall1_cylinder = openmc.ZCylinder(x0=0.0, y0=0.0, R=R2)
vessel_cylinder = openmc.ZCylinder(x0=0.0, y0=0.0, R=R3)
shield_cylinder = openmc.ZCylinder(x0=0.0, y0=0.0, R=R4)
wall2_cylinder = openmc.ZCylinder(x0=0.0, y0=0.0, R=R5)
tfcoil_cylinder = openmc.ZCylinder(x0=0.0, y0=0.0, R=R6)

med = openmc.XPlane(x0=0)
trait_haut = openmc.YPlane(y0=R1)
trait_bas = openmc.YPlane(y0=-R1)

circle1_1 = openmc.ZCylinder(x0=0, y0=R1, R=D1)
circle1_2 = openmc.ZCylinder(x0=0, y0=R1, R=D2)
circle1_3 = openmc.ZCylinder(x0=0, y0=R1, R=D3)
circle1_4 = openmc.ZCylinder(x0=0, y0=R1, R=D4)
circle1_5 = openmc.ZCylinder(x0=0, y0=R1, R=D5)

circle2_1 = openmc.ZCylinder(x0=0, y0=-R1, R=D1)
circle2_2 = openmc.ZCylinder(x0=0, y0=-R1, R=D2)
circle2_3 = openmc.ZCylinder(x0=0, y0=-R1, R=D3)
circle2_4 = openmc.ZCylinder(x0=0, y0=-R1, R=D4)
circle2_5 = openmc.ZCylinder(x0=0, y0=-R1, R=D5)

pla1 = openmc.XPlane(x0=-D1)
pla2 = openmc.XPlane(x0=-D2)
pla3 = openmc.XPlane(x0=-D3)
pla4 = openmc.XPlane(x0=-D4)
pla5 = openmc.XPlane(x0=-D5)

target_zone = +trait_target1 & -trait_target2 & +trait_target3 & -pla4
no_target_zone = (-trait_target1 & -trait_target2 & +trait_target3) | +trait_target2 | -trait_target3

zone1_1 = -plasma_cylinder & +med
zone1_2 = -wall1_cylinder & +med & +plasma_cylinder
zone1_3 = -vessel_cylinder & +med & +wall1_cylinder
zone1_4 = -shield_cylinder & +med & +vessel_cylinder
zone1_5 = -wall2_cylinder & +med & +shield_cylinder
zone1_6 = -tfcoil_cylinder & +med & +wall2_cylinder

zone2_2 = -med & +trait_haut & -circle1_1 
zone2_3 = -med & +trait_haut & -circle1_2 & +circle1_1
zone2_4 = -med & +trait_haut & -circle1_3 & +circle1_2
zone2_5 = -med & +trait_haut & -circle1_4 & +circle1_3
zone2_6 = -med & +trait_haut & -circle1_5 & +circle1_4

zone3_2 = -med & -trait_bas & -circle2_1 
zone3_3 = -med & -trait_bas & -circle2_2 & +circle2_1
zone3_4 = -med & -trait_bas & -circle2_3 & +circle2_2
zone3_5 = -med & -trait_bas & -circle2_4 & +circle2_3
zone3_6 = -med & -trait_bas & -circle2_5 & +circle2_4


zone4_2 = +trait_bas & -trait_haut & +pla1 & -med
zone4_3 = +trait_bas & -trait_haut & +pla2 & -pla1
zone4_4 = +trait_bas & -trait_haut & +pla3 & -pla2
zone4_5 = +trait_bas & -trait_haut & +pla4 & -pla3
zone4_6 = +trait_bas & -trait_haut & +pla5 & -pla4 & no_target_zone


out1 = +tfcoil_cylinder & +med
out2 = +circle1_5 & -med & +trait_haut
out3 = +circle2_5 & -med & -trait_bas
out4 = -pla5

#Associate region
plasma.region = zone1_1 & +border_angle_1 & -border_angle_2 & +min_x & -max_x & +min_y & -max_y
wall1.region = (zone1_2 | zone2_2 | zone3_2 | zone4_2) & +border_angle_1 & -border_angle_2 & +min_x & -max_x & +min_y & -max_y
vessel.region = (zone1_3 | zone2_3 | zone3_3 | zone4_3) & +border_angle_1 & -border_angle_2 & +min_x & -max_x & +min_y & -max_y
shield.region = (zone1_4 | zone2_4 | zone3_4 | zone4_4) & +border_angle_1 & -border_angle_2 & +min_x & -max_x & +min_y & -max_y
wall2.region = (zone1_5 | zone2_5 | zone3_5 | zone4_5) & +border_angle_1 & -border_angle_2 & +min_x & -max_x & +min_y & -max_y
tfcoil.region = (zone1_6 | zone2_6 | zone3_6 | zone4_6) & +border_angle_1 & -border_angle_2 & +min_x & -max_x & +min_y & -max_y
target.region = target_zone & +border_angle_1 & -border_angle_2 & +min_x & -max_x & +min_y & -max_y
out.region = (out1 | out2 | out3 | out4) & +border_angle_1 & -border_angle_2 & +min_x & -max_x & +min_y & -max_y

#fill cells
plasma.fill = plasma_H
wall1.fill = tungsten
vessel.fill = ss316
shield.fill = B_carbide
wall2.fill = ss316
tfcoil.fill = incoloy


#Add to the universe
universe = openmc.Universe()
universe.add_cell(plasma)
universe.add_cell(wall1)
universe.add_cell(vessel)
universe.add_cell(shield)
universe.add_cell(wall2)
universe.add_cell(tfcoil)
universe.add_cell(target)
universe.add_cell(out)

#export XML file
geom = openmc.Geometry(universe)
geom.export_to_xml()

#####################
#simulation settings#
#####################

" The source must be define as a fusion source"
" That's why we build an isotropic source with 14.1Mev energy neutron"
" It is different than a fission source so we have to define a fixed source"
#Create the source and XML file
point = openmc.stats.Point((0, 0, 0))
src = openmc.Source(space=point)
src.angle = openmc.stats.Isotropic()
src.energy = openmc.stats.Discrete([14.1e6], [1.0])
settings = openmc.Settings()
settings.run_mode = "fixed source"

"With the geometry and materials finished, we now just need to define simulation parameters."
"Here, inactive batch are useless"
settings.source = src
settings.batches = N_batch
settings.inactive = 0
settings.particles = N_particules
settings.export_to_xml()

################
#Geometry Plot#
###############
"In order to watch if our geometry is good, we can plot it in an array by using the method of 'in cell'"
"This method allow us to plot the geometry in a space even if we don't know if the geometry will run"
"With this method we can plot just one cells, combinaison of many cells or just all cells, it is really usefull when we want to watch if it is working well"

"Here, we plot the xy face in z=0"
#xy plot
B = zeros((400,border_x_1+border_x_2))
for i in range (-border_x_1,border_x_2,1):
    for j in range (-200,200,1):
        if (i,j,0) in plasma:
            B[j+200][i+border_x_1]=1
        elif (i,j,0) in wall1:
            B[j+200][i+border_x_1]=2
        elif (i,j,0) in vessel:
            B[j+200][i+border_x_1]=3
        elif (i,j,0) in shield:
            B[j+200][i+border_x_1]=4
        elif (i,j,0) in wall2:
            B[j+200][i+border_x_1]=5
        elif (i,j,0) in tfcoil:
            B[j+200][i+border_x_1]=6
        elif (i,j,0) in target:
            B[j+200][i+border_x_1]=7

"Here, we plot the xz face in y=0"             
#xz plot
C = zeros((200,border_x_1+border_x_2))
for i in range (-border_x_1,border_x_2,1):
    for j in range (-100,100,1):
        if (i,0,j) in plasma :
            C[j+100][i+border_x_1]=1
        elif (i,0,j) in wall1:
            C[j+100][i+border_x_1]=2
        elif (i,0,j) in vessel:
            C[j+100][i+border_x_1]=3
        elif (i,0,j) in shield:
            C[j+100][i+border_x_1]=4
        elif (i,0,j) in wall2:
            C[j+100][i+border_x_1]=5
        elif (i,0,j) in tfcoil:
            C[j+100][i+border_x_1]=6
        elif (i,0,j) in target:
            C[j+100][i+border_x_1]=7
            
"Ones our array are build, we can use the subplot of plt in order to watch our two geometry"
fig = plt.subplot(121)
fig.imshow(B)
fig2 = plt.subplot(122)
fig2.imshow(C)
plt.savefig("geometry_xy_xz") 

#####################
#Filters and tallies#
#####################

"We did all this work with one goal : mesure the heat deposition in some cells"
"Now, will setup all tallies to be able to make mesure"

"First we create the file"
tallies_file=openmc.Tallies()

"Than we have to instantiate our 3D size in order to make a mesh"
size_x = int(border_x_1+border_x_2/mesh_size)
size_y = int(400/mesh_size)
size_z = int(200/mesh_size)

"At the beginning of our code we define the mesh_size : more it will be low, more our mesure will be"
"precise in the space. But we have to be carefull on the uncertainty!"
# Instantiate a tally mesh
mesh = openmc.Mesh(mesh_id=1)
mesh.type = 'regular'
mesh.dimension = [size_x,size_y,size_z]
mesh.lower_left = [-border_x_1, -200.0, -100.0]
mesh.width = [mesh_size,mesh_size,mesh_size]
mesh_filter = openmc.MeshFilter(mesh)

"First, we want to mesure the flux and the heat deposition. We create this tally for that"
"This tally is using our mesh, to mesure the heat deposition in 'box' sizing by mesh_size"
# Instantiate flux Tally in moderator and fuel
tally = openmc.Tally(name='flux')
tally.filters = [mesh_filter]
tally.scores = ['flux','heating']
tallies_file.append(tally)

"Than, we wants to mesure the mean of heat deposition in the cells in order to have a macroscopic idea of the values"
# Instantiate absorption in all cells except shield
cells_abs = openmc.Tally(name='cell_absorption')
cells_abs.filters = [openmc.CellFilter([plasma,wall1,vessel,shield,wall2,tfcoil,target,out])]
cells_abs.scores = ['heating']
tallies_file.append(cells_abs)

"To end the tallies, we just have to export in xml"
# Export to "tallies.xml"
tallies_file.export_to_xml()

"All is now complete, we can run !"
#Run
openmc.run()


###################
# POST PROCESSING #
###################
"We do not know how many batches were needed to satisfy the tally trigger(s), so find the statepoint file(s)"
"This command allow us to find the good statepoint file and load it"
statepoints = glob.glob('statepoint.*.h5')
# Load the last statepoint file
sp = openmc.StatePoint(statepoints[-1])

"Than, we want to get the tallies named flux who contain the mesh of heat depostion"
#Get the flux in all fusion cylinder
flux = sp.get_tally(name='flux')
flux.get_pandas_dataframe()

"We get the others tallies to know the heat deposition in all cells"
#Get absorption in the cells
absorp = sp.get_tally(name='cell_absorption')
heat_deposition = absorp.get_pandas_dataframe()

"Than, we have to get the specific mesh tallies of heating"
#For ploting
flux_a = flux.get_slice(scores=['heating'])
"We have to reshape this array with the mesh we instantiate"
flux_a.std_dev.shape = (size_x,size_y,size_z)
flux_a.mean.shape = (size_x,size_y,size_z)

"This program all us to know the maximum relative error in the target in order to know if we need to increase batches and particules number to reduce it" 
"So, to do that, we are going to mesure all relative error in target and print the maximum"
# Determine max relative error
b=0
for i in range (0,size_x-1):
    for j in range (0,size_y-1):
        for k in range (0,size_z-1):
            if (i-size_x,j-size_y,k-size_z) in target:
                if flux_a.mean[i][j][k]>0:
                    relative_error = (flux_a.std_dev[i][j][k] / flux_a.mean[i][j][k])   
                    if relative_error>b:
                        b=relative_error
print("Maximum relative error is",b)


"In order to make the plot more visible, we plot something like the log of the real value"

#for ploting xy
"For ploting an xy slice of our 3D model, we choose to be center on z axis"
D = zeros((400,border_x_1+border_x_2))
for i in range (0,border_x_1+border_x_2-1):
    for j in range (0,400-1):
        if flux_a.mean[i][j][100]!=0 :
            D[j][i]= log(flux_a.mean[i][j][100])

#for ploting xz
"For ploting an xz slice of our 3D model, we choose to be center on y axis"
E = zeros((200,border_x_1+border_x_2))
for i in range (0,border_x_1+border_x_2-1):
    for j in range (0,200-1):
        if flux_a.mean[i][200][j]!=0 :
            E[j][i] = log(flux_a.mean[i][200][j])
        
            
"We can now print the heat deposition of all cells in order to match the mean value in an array"
print(heat_deposition)

"We can plot the two slice of our 3D model to watch the heat deposition in the whole reactor"
fig = plt.subplot(121)
fig.imshow(D)
fig2 = plt.subplot(122)
fig2.imshow(E)
plt.savefig("heating.pdf") 



