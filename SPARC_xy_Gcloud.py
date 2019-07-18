import matplotlib.pyplot as plt
from math import *
from openmc import *
from numpy import*

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
N_batch = 50
N_particules = 1000
#mesh for plotting : precision of all
mesh_size = 1



size = int(600/mesh_size)

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
air_in=openmc.Cell(8,"air_in")
n_phage=openmc.Cell(9,"n_phage")

#create border
min_x=openmc.XPlane(x0=-300,boundary_type='reflective')
max_x=openmc.XPlane(x0=300,boundary_type='reflective')
min_y=openmc.YPlane(y0=-300,boundary_type='reflective')
max_y=openmc.YPlane(y0=300,boundary_type='reflective')

#create cylinder
R1 = e_plasma
R2 = e_plasma+e_wall1
R3 = e_plasma+e_wall1+e_vessel
R4 = e_plasma+e_wall1+e_vessel+e_shield
R5 = e_plasma+e_wall1+e_vessel+e_shield+e_wall2
R6 = R5 + e_tfcoil
center_circle = R5/2


plasma_cylinder = openmc.ZCylinder(x0=0.0, y0=0.0, R=R1)
wall1_cylinder = openmc.ZCylinder(x0=0.0, y0=0.0, R=R2)
vessel_cylinder = openmc.ZCylinder(x0=0.0, y0=0.0, R=R3)
shield_cylinder = openmc.ZCylinder(x0=0.0, y0=0.0, R=R4)
wall2_cylinder = openmc.ZCylinder(x0=0.0, y0=0.0, R=R5)

R_small_circle_1 = R5-center_circle
R_small_circle_2 = R5-center_circle+e_tfcoil
R_big_circle_1 = center_circle+R5
R_big_circle_2 = center_circle+R5+e_tfcoil

circle1_1 = openmc.ZCylinder(x0=-center_circle, y0=R5, R=R_small_circle_1)
circle1_2 = openmc.ZCylinder(x0=-center_circle, y0=R5, R=R_small_circle_2)
circle1_3 = openmc.ZCylinder(x0=-center_circle, y0=-R5, R=R_small_circle_1)
circle1_4 = openmc.ZCylinder(x0=-center_circle, y0=-R5, R=R_small_circle_2)
circle2_1 = openmc.ZCylinder(x0=-center_circle, y0=0.0, R=R_big_circle_1)
circle2_2 = openmc.ZCylinder(x0=-center_circle, y0=0.0, R=R_big_circle_2)


sol_in = openmc.XPlane(x0=-R5)
sol_out = openmc.XPlane(x0=-R6)
trait_1 = openmc.XPlane(x0=-center_circle)
trait_2 = openmc.YPlane(y0=+R5)
trait_3 = openmc.YPlane(y0=-R5)
trait_target1 = openmc.XPlane(x0=-R5-e_target)
trait_target2 = openmc.YPlane(y0=H_target/2)
trait_target3 = openmc.YPlane(y0=-H_target/2)

zone1 = +circle2_1 & -circle2_2 & +trait_1
zone2 = +circle1_1 & -circle1_2 & -trait_1 & +trait_2
zone3 = -sol_in & +sol_out & -trait_2 & +trait_3
zone4 = +circle1_3 & -circle1_4 & -trait_1 & -trait_3

inter_1 = -circle2_1 & +wall2_cylinder & +trait_1
inter_2 = +sol_in & +wall2_cylinder & -trait_2 & +trait_3 & -trait_1
inter_3 = -circle1_1 & +sol_in & +trait_2 & -trait_1 & +wall2_cylinder
inter_4 = -circle1_3 & +sol_in & -trait_3 & -trait_1 & +wall2_cylinder
inter = (inter_1 | inter_2 | inter_3 | inter_4) 
exter = +circle2_2 & +trait_1 | +circle1_2 & +sol_out & -trait_1 & +trait_2 | +circle1_4 & +sol_out & -trait_1 & -trait_3 | -sol_out

#Associate region
plasma.region = -plasma_cylinder
wall1.region = +plasma_cylinder & -wall1_cylinder
vessel.region = +wall1_cylinder & -vessel_cylinder
shield.region = +vessel_cylinder & -shield_cylinder
wall2.region = +shield_cylinder & -wall2_cylinder
tfcoil.region = zone1 | zone2 | zone3 | zone4 
target.region = +trait_target1 & -trait_target2 & +trait_target3 & -sol_in
air_in.region = inter
n_phage.region = exter & +min_x & -max_x & +min_y & -max_y

#fill cells
plasma.fill = plasma_H
wall1.fill = tungsten
vessel.fill = ss316
shield.fill = B_carbide
wall2.fill = ss316
tfcoil.fill = incoloy
target.fill = incoloy
#air_in.fill = air
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
universe.add_cell(air_in)
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
B = zeros((600,600))
for i in range (-300,300,1):
    for j in range (-300,300,1):
        if (i,j,0) in plasma:
            B[j+300][i+300]=1
        elif (i,j,0) in wall1:
            B[j+300][i+300]=2
        elif (i,j,0) in vessel:
            B[j+300][i+300]=3
        elif (i,j,0) in shield:
            B[j+300][i+300]=4
        elif (i,j,0) in wall2:
            B[j+300][i+300]=5
        elif (i,j,0) in tfcoil:
            B[j+300][i+300]=6
        elif (i,j,0) in target:
            B[j+300][i+300]=7
        elif (i,j,0) in air_in:
            B[j+300][i+300]=0
        elif (i,j,0) in n_phage:
            B[j+300][i+300]=9

            
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
mesh.lower_left = [-300.0, -300.0]
mesh.width = [mesh_size,mesh_size]
mesh_filter = openmc.MeshFilter(mesh)

# Instantiate flux Tally in moderator and fuel
tally = openmc.Tally(name='flux')
tally.filters = [mesh_filter]
tally.scores = ['flux','heating']
tallies_file.append(tally)

# Instantiate absorption in all cells except shield
cells_abs = openmc.Tally(name='cell_absorption')
cells_abs.filters = [openmc.CellFilter([plasma,wall1,vessel,shield,wall2,tfcoil,target,air_in,n_phage])]
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
B = zeros((600,600))
for i in range (-300,300,1):
    for j in range (-300,300,1):
        if (i,j,0) in target:
            if relative_error[j+300][i+300]>b:
                b=relative_error[j+300][i+300]
print("max error is",b)

print(heat_deposition)

for i in range (-300,300,1):
    for j in range (-300,300,1):
        if (i,j,0) in target:
            print("Heating (W/sources particules) in this part is :",flux_a.mean[j+100][i]*1.6*10**-19)

#for ploting
for i in range (0,size-1):
    for j in range (0,size-1):
        if flux_a.mean[i][j]!=0 :
            flux_a.mean[i][j]=log(flux_a.mean[i][j])
    
plt.imshow(flux_a.mean)
plt.savefig("heating.pdf")   

