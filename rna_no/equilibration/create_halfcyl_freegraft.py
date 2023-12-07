#Creates two polymer-grafted cylinders in a cubic box with additional free polymers.
#In this version, we also generate some free polymers of the same type as the grafted ones.

import os, sys, time
from numpy import *
from numpy.random import *
#import matplotlib.pyplot as plt

def rotate(vector,angle): #Rotate vector by angle
	r=array([cos(angle),-sin(angle),sin(angle),cos(angle)])
	rot=reshape(r,(2,2)) #2d rotation matrix
	return rot.dot(vector)
	
bond_length=0.96
nspace=5
lattice_graft=float(nspace+1)
surface_mass=0.1 #mass of cylinder surface particles

polymer_length=24 		#length of grafted polymers
npoly_free_graft=100    #number of free polymers (of same type as the grafted ones)

#B5 bead charge:
#b5_charge=float(loadtxt("b5_charge.dat"))
b5_charge=0 #This value is overwritten by the "set charge" command

#Mean Ki-67 charge excluding B5 bead:
#poly_charge=float(loadtxt("poly_charge.dat"))
poly_charge=0 #This value is overwritten by the "set charge" command

#surface_charge=float(loadtxt("surface_charge.dat"))
surface_charge=0 #This value is overwritten by the "set charge" command

n_b5=20 #Counting from the binding site, the B5 bead is the 20th bead

#n_graft=int(input("Enter n. of grafted polymers per ring:"))
#n_grafted_layers=int(input("Enter n. of layers of grafted polymers:"))

n_graft=18
n_grafted_layers=5

theta_graft=2*pi/n_graft
radius=lattice_graft/theta_graft
theta_surf=theta_graft/(nspace+1)
lattice_surf=theta_surf*radius

natoms_ring=(nspace+1)*n_graft 				#total number of atoms per ring
n_rings=(nspace+1)*n_grafted_layers+nspace 	#number of rings
n_grafted_poly=n_grafted_layers*(n_graft-2) #number of grafted polymers
height=n_grafted_layers*lattice_graft 		#cylinder height

box_yz=height*2.5 							#box size in y and z directions
box_x=height*4  					 		#box size in x direction
center_distance=2*radius+polymer_length+2 	#distance between centers of the two cylinders

natoms_surf= 2*((int((natoms_ring+1)/2.-1)*n_rings))
natoms_grafted_poly= (n_grafted_poly*polymer_length)
natoms= natoms_surf + natoms_grafted_poly + npoly_free_graft*polymer_length
nbonds= natoms_grafted_poly + npoly_free_graft*(polymer_length-1)
nangles= n_grafted_poly*(polymer_length-1) + npoly_free_graft*(polymer_length-2)

with open('halfcyl_freegraft_ngraft%d_nlayers%d_npolyfreegraft%d.lammpsdata'%(n_graft,n_grafted_layers,npoly_free_graft),'w') as fout:
	fout.write("LAMMPS data file\n\n%d atoms\n5 atom types\n%d bonds\n1 bond types\n%d angles\n1 angle types\n\n0.0 %f xlo xhi\n0.0 %f ylo yhi\n0.0 %f zlo zhi\n\nMasses\n\n1 %f\n2 1\n3 1\n4 1\n5 %f\n\nAtoms\n\n"%(natoms,nbonds,nangles,box_x,box_yz,box_yz,surface_mass,surface_mass))

	atom_id=0
	mol_id=0
	bonds=[]
	angles=[]

	#Generate first cylinder
	zstart=(box_yz-height)/2.

	for k in range(1,n_rings+1):
		x0_cyl1=(box_x-center_distance)/2.
		x0=x0_cyl1
		y0=box_yz/2.
		z0=zstart+k*lattice_surf
		r0=array([x0,y0,z0])

		#Create 1 layer (of cylinder's surface)
		vec=array([0,-1])
		for n in range(1,int((natoms_ring+1)/2.)):
			if( ((n%(nspace+1))!=0) or k%(nspace+1)!=0):
				atom_id+=1
				rotvec=rotate(vec,n*theta_surf)
				r=r0+radius*array([rotvec[0],rotvec[1],0])
				fout.write("%d 1 1 %f %.15f %.15f %.15f\n"%(atom_id,surface_charge,r[0],r[1],r[2]))

			#elif( ((n%(nspace+1))==0) and k%(nspace+1)==0):
			else:
				#Print binding site
				mol_id+=1
				atom_id+=1
				rotvec=rotate(vec,n*theta_surf)
				r=r0+radius*array([rotvec[0],rotvec[1],0])
				fout.write("%d 1 1 %f %.15f %.15f %.15f\n"%(atom_id,surface_charge,r[0],r[1],r[2]))

				#Print rest of chain
				for j in range(1,polymer_length+1):
					atom_id+=1
					bonds.append((atom_id,atom_id-1))
					if(j>1):
						angles.append((atom_id-2,atom_id-1,atom_id))
					r=r0+(radius+j*bond_length)*array([rotvec[0],rotvec[1],0])
					if(j!=n_b5):
						fout.write("%d %d 2 %f %.15f %.15f %.15f\n"%(atom_id,mol_id,poly_charge,r[0],r[1],r[2]))
					else:
						fout.write("%d %d 4 %f %.15f %.15f %.15f\n"%(atom_id,mol_id,b5_charge,r[0],r[1],r[2]))

	#Generate second cylinder
	zstart=(box_yz-height)/2.

	for k in range(1,n_rings+1):
		x0_cyl2=(box_x+center_distance)/2.
		x0=x0_cyl2		
		y0=box_yz/2.
		z0=zstart+k*lattice_surf
		r0=array([x0,y0,z0])

		#Create 1 layer (of cylinder's surface)
		vec=array([0,1])
		for n in range(1,int((natoms_ring+1)/2.)):
			if( ((n%(nspace+1))!=0) or k%(nspace+1)!=0):
				atom_id+=1
				rotvec=rotate(vec,n*theta_surf)
				r=r0+radius*array([rotvec[0],rotvec[1],0])
				fout.write("%d 1 5 %f %.15f %.15f %.15f\n"%(atom_id,surface_charge,r[0],r[1],r[2]))

			#elif( ((n%(nspace+1))==0) and k%(nspace+1)==0):
			else:
				#Print binding site
				mol_id+=1
				atom_id+=1
				rotvec=rotate(vec,n*theta_surf)
				r=r0+radius*array([rotvec[0],rotvec[1],0])
				fout.write("%d 1 5 %f %.15f %.15f %.15f\n"%(atom_id,surface_charge,r[0],r[1],r[2]))

				#Print rest of chain
				for j in range(1,polymer_length+1):
					atom_id+=1
					bonds.append((atom_id,atom_id-1))
					if(j>1):
						angles.append((atom_id-2,atom_id-1,atom_id))
					r=r0+(radius+j*bond_length)*array([rotvec[0],rotvec[1],0])
					if(j!=n_b5):
						fout.write("%d %d 2 %f %.15f %.15f %.15f\n"%(atom_id,mol_id,poly_charge,r[0],r[1],r[2]))
					else:
						fout.write("%d %d 4 %f %.15f %.15f %.15f\n"%(atom_id,mol_id,b5_charge,r[0],r[1],r[2]))

	#Print free polymers
#	for poly_id in range(npoly_free):
#		atom_id+=1
#		mol_id+=1
#		z=random()*box_yz
#		x=random()*box_x
#		while(((x>x0_cyl1) and (x<x0_cyl1+(radius+0.5))) or (((x>x0_cyl2-(radius+0.5)) and (x<x0_cyl2)))): #Pick again if inside cylinder
#			x=random()*box_x
#		y=random()*box_yz
#		while(((y>box_yz/4.-(radius+0.5)) and (y<box_yz/4.+(radius+0.5))) or (((y>3*box_yz/4.-(radius+0.5)) and (y<3*box_yz/4.+(radius+0.5))))): #Pick again if inside cylinder
#			y=random()*box_yz
#		r0=array([x,y,z])
#		fout.write("%d %d 3 %f %.15f %.15f %.15f\n"%(atom_id,mol_id,poly_charge_free,x,y,z))
#		for n in range(1,polymer_length_free):
#			atom_id+=1
#			r=r0+bond_length*n*array([0,pi/4.,pi/4.])
#			x=r[0]
#			y=r[1]
#			z=r[2]
#			bonds.append((atom_id,atom_id-1))
#			if(n>1):
#				angles.append((atom_id-2,atom_id-1,atom_id))
#			fout.write("%d %d 3 %f %.15f %.15f %.15f\n"%(atom_id,mol_id,poly_charge_free,x,y,z))

	#Print free polymers (of same type as the grafted ones)
	for poly_id in range(npoly_free_graft):
		atom_id+=1
		mol_id+=1
		z=random()*box_yz
		x=random()*box_x
		while(((x>x0_cyl1) and (x<x0_cyl1+(radius+0.5))) or (((x>x0_cyl2-(radius+0.5)) and (x<x0_cyl2)))): #Pick again if inside cylinder
			x=random()*box_x
		y=random()*box_yz
#		while(((y>box_yz/4.-(radius+0.5)) and (y<box_yz/4.+(radius+0.5))) or (((y>3*box_yz/4.-(radius+0.5)) and (y<3*box_yz/4.+(radius+0.5))))): #Pick again if inside cylinder
#			y=random()*box_yz
		r0=array([x,y,z])
		fout.write("%d %d 2 %f %.15f %.15f %.15f\n"%(atom_id,mol_id,poly_charge,x,y,z))
		for n in range(1,polymer_length):
			atom_id+=1
			r=r0+bond_length*n*array([0,pi/4.,pi/4.])
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

	if(len(bonds)!=nbonds):
		print("WARNING: Number of bonds different from expected number of bonds (have %d, expect %d)"%(len(bonds),nbonds))

	if(len(angles)!=nangles):
		print("WARNING: Number of angles different from expected number of angles (have %d, expect %d)"%(len(angles),nangles))

	#Print bonds
	fout.write("\nBonds\n\n")
	for nb, bond in enumerate(bonds):
		fout.write("%d 1 %d %d\n"%(nb+1,bond[0],bond[1]))

	#Print angles
	fout.write("\nAngles\n\n")
	for na, angle in enumerate(angles):
		fout.write("%d 1 %d %d %d\n"%(na+1,angle[0],angle[1],angle[2]))









