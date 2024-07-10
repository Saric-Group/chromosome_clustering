#Create two crystal surfaces with regularly spaced grafted polymer chains.
#The grafted portions of each surface face each other.

#The crystal surface is assumed to be a square lattice.
#The polymer chains are also arranged on a square lattice.
#The polymer chains are semiflexible (angular potential).

#This version reads the exact charge distribution of the polymers from an input file.
#The charge distribution is assumed to be the same for every polymer.
#Each row of the file is expected to contain an integer charge value.

#This version also creates a certain number of free polymers.

import os, sys, time
from numpy import *
from numpy.random import *
#import matplotlib.pyplot as plt

bond_length=0.96
surface_density0=1.

npoly0=100
polymer_length=24

npoly_free=24
polymer_length_free=100
#poly_charge_free=-3.00

#surface_charge=-0.50
#nspace=6

dz_surf=1e-3
dz_free=1.
z0_poly_free=dz_free+dz_surf

#B5 bead charge:
#b5_charge=float(loadtxt("b5_charge.dat"))
b5_charge=0 #This value is overwritten by the "set charge" command

#Mean Ki-67 charge excluding B5 bead:
#poly_charge=float(loadtxt("poly_charge.dat"))
poly_charge=0 #This value is overwritten by the "set charge" command

#surface_charge=float(loadtxt("surface_charge.dat"))
surface_charge=0 #This value is overwritten by the "set charge" command

#poly_charge_free=float(loadtxt("poly_charge_free.dat"))
poly_charge_free=0 #This value is overwritten by the "set charge" command

n_b5=20 #Counting from the binding site, the B5 bead is the 20th bead

#nspace=int(loadtxt("nspace.dat"))
box_z=float(loadtxt("box_z.dat"))
#nspace=int(loadtxt("nspace.dat"))
nspace=6

box_z_ext=box_z+2*dz_surf #needed to use fixed boundary conditions in z direction
theta=arcsin((box_z-2*dz_free)/(bond_length*polymer_length_free))

#box_z=float(input("Enter box lenght in z direction:"))
#npoly0=int(input("Enter n. of grafted polymers:"))
#polymer_length=int(input("Enter grafted polymer length:"))
#surface_density0=float(input("Enter surface density:"))
#surface_charge=float(input("Enter elementary surface charge (with sign!):"))
#nspace=int(input("Enter desired lattice spacing between two grafted polymers (in x or y direction):"))

lattice=1/sqrt(surface_density0)
lattice_graft=(nspace+1)*lattice
npoly_side=round(sqrt(npoly0))
npoly=npoly_side**2
nsurface_side=npoly_side*(nspace+1)
nsurface=nsurface_side**2
box=lattice*nsurface_side
surface_density=nsurface/(box**2)

nsurface_tot=2*nsurface
npoly_tot=2*npoly

#NB the grafted polymers have 1 more bond and 1 more angle since they're attached to the surface

natoms= nsurface_tot + npoly_tot*polymer_length + npoly_free*polymer_length_free
nbonds= npoly_tot*polymer_length + npoly_free*(polymer_length_free-1)
nangles= npoly_tot*(polymer_length-1) + npoly_free*(polymer_length_free-2)

#A vector we will use to print the coordinates of the polymers
vec=array([0,cos(theta),sin(theta)])

with open('surfaces_freepoly_Lz%.2f_npoly%d_Npoly%d_npolyfree%d_Npolyfree%d_nspace%d.lammpsdata'%(box_z,npoly,polymer_length,npoly_free,polymer_length_free,nspace),'w') as fout:
	fout.write("LAMMPS data file\n\n%d atoms\n4 atom types\n%d bonds\n1 bond types\n%d angles\n1 angle types\n\n0.0 %f xlo xhi\n0.0 %f ylo yhi\n0.0 %f zlo zhi\n\nMasses\n\n1 1\n2 1\n3 1\n4 1\n\nAtoms\n\n"%(natoms,nbonds,nangles,box,box,box_z_ext))
	
	atom_id=0

	#Print 1st layer (fixed surface -- type 1) excluding binding sites
	for ix in range(nsurface_side):
		for iy in range(nsurface_side):
			if((ix%(nspace+1)!=0 or iy%(nspace+1)!=0)):
				atom_id+=1
				x=ix*lattice
				y=iy*lattice
				fout.write("%d 1 1 %f %.15f %.15f %.15f\n"%(atom_id,surface_charge,x,y,dz_surf))

	#Print 2nd layer (movable surface -- type 4) excluding binding sites
	for ix in range(nsurface_side):
		for iy in range(nsurface_side):
			if((ix%(nspace+1)!=0 or iy%(nspace+1)!=0)):
				atom_id+=1
				x=ix*lattice
				y=iy*lattice
				fout.write("%d 1 1 %f %.15f %.15f %.15f\n"%(atom_id,surface_charge,x,y,box_z+dz_surf))

	bonds=[]
	angles=[]
	mol_id=0

	#Print 1st layer grafted polymers, including binding sites
	for ix in range(npoly_side):
		for iy in range(npoly_side):
			#Print binding site (type 1)
			atom_id+=1
			mol_id+=1
			x=ix*lattice_graft
			y=iy*lattice_graft
			r0=array([x,y,dz_surf])
			fout.write("%d 1 1 %f %.15f %.15f %.15f\n"%(atom_id,surface_charge,x,y,dz_surf))
		
			#Print rest of the chain (type 2)
			for n in range(1,polymer_length+1):
				atom_id+=1
				r=r0+bond_length*n*vec
				x=r[0]
				y=r[1]
				z=r[2]				
				bonds.append((atom_id,atom_id-1))
				if(n>1):
					angles.append((atom_id-2,atom_id-1,atom_id))
				if(n!=n_b5):
					fout.write("%d %d 2 %f %.15f %.15f %.15f\n"%(atom_id,mol_id,poly_charge,x,y,z))
				else:
					fout.write("%d %d 4 %f %.15f %.15f %.15f\n"%(atom_id,mol_id,b5_charge,x,y,z))

	#Print 2nd layer grafted polymers, including binding sites
	for ix in range(npoly_side):
		for iy in range(npoly_side):
			#Print binding site (type 1)
			atom_id+=1
			mol_id+=1
			x=ix*lattice_graft
			y=iy*lattice_graft
			r0=array([x,y,box_z+dz_surf])
			fout.write("%d 1 1 %f %.15f %.15f %.15f\n"%(atom_id,surface_charge,x,y,box_z+dz_surf))
		
			#Print rest of the chain (type 2)
			for n in range(1,polymer_length+1):
				atom_id+=1
				r=r0-bond_length*n*vec
				x=r[0]
				y=r[1]
				z=r[2]				
				bonds.append((atom_id,atom_id-1))
				if(n>1):
					angles.append((atom_id-2,atom_id-1,atom_id))
				if(n!=n_b5):
					fout.write("%d %d 2 %f %.15f %.15f %.15f\n"%(atom_id,mol_id,poly_charge,x,y,z))
				else:
					fout.write("%d %d 4 %f %.15f %.15f %.15f\n"%(atom_id,mol_id,b5_charge,x,y,z))

	#Print free polymers
	for poly_id in range(npoly_free):
		atom_id+=1
		mol_id+=1
		z=z0_poly_free
		x=random()*box
		y=random()*box
		r0=array([x,y,z])
		fout.write("%d %d 3 %f %.15f %.15f %.15f\n"%(atom_id,mol_id,poly_charge_free,x,y,z))
		for n in range(1,polymer_length_free):
			atom_id+=1
			r=r0+bond_length*n*vec
			x=r[0]
			y=r[1]
			z=r[2]
			bonds.append((atom_id,atom_id-1))
			if(n>1):
				angles.append((atom_id-2,atom_id-1,atom_id))
			fout.write("%d %d 3 %f %.15f %.15f %.15f\n"%(atom_id,mol_id,poly_charge_free,x,y,z))

	if(len(bonds)!=nbonds):
		print("WARNING: Number of bonds different from expected number of bonds (have %d, expect %d)"%(len(bonds),nbonds))

	#Print bonds
	fout.write("\nBonds\n\n")
	for nb, bond in enumerate(bonds):
		fout.write("%d 1 %d %d\n"%(nb+1,bond[0],bond[1]))

	#Print angles
	fout.write("\nAngles\n\n")
	for na, angle in enumerate(angles):
		fout.write("%d 1 %d %d %d\n"%(na+1,angle[0],angle[1],angle[2]))




