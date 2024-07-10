#Creates free Ki-67 and RNA polymers.

import os, sys, time
from numpy import *
from numpy.random import *
#import matplotlib.pyplot as plt

def rotate(vector,angle): #Rotate vector by angle
	r=array([cos(angle),-sin(angle),sin(angle),cos(angle)])
	rot=reshape(r,(2,2)) #2d rotation matrix
	return rot.dot(vector)
	
bond_length=0.96

polymer_length=24 		#length of grafted polymers
npoly_free_graft=400    #number of free polymers (of same type as the grafted ones)

npoly_free=100			#number of free polymers
polymer_length_free=100 #length of free polymers

#B5 bead charge:
#b5_charge=float(loadtxt("b5_charge.dat"))
b5_charge=0 #This value is overwritten by the "set charge" command

#Mean Ki-67 charge excluding B5 bead:
#poly_charge=float(loadtxt("poly_charge.dat"))
poly_charge=0 #This value is overwritten by the "set charge" command

#poly_charge_free=float(loadtxt("poly_charge_free.dat"))
poly_charge_free=0 #This value is overwritten by the "set charge" command

n_b5=20 #Counting from the binding site, the B5 bead is the 20th bead

natoms= npoly_free*polymer_length_free + npoly_free_graft*polymer_length
nbonds= npoly_free*(polymer_length_free-1) + npoly_free_graft*(polymer_length-1)
nangles= npoly_free*(polymer_length_free-2) + npoly_free_graft*(polymer_length-2)

rho=float(loadtxt("rho.dat"))
box=(natoms/rho)**(1/3.)

with open('freepoly_freegraft_npolyfreegraft%d_npolyfree%d_Npolyfree%d.lammpsdata'%(npoly_free_graft,npoly_free,polymer_length_free),'w') as fout:
	fout.write("LAMMPS data file\n\n%d atoms\n4 atom types\n%d bonds\n1 bond types\n%d angles\n1 angle types\n\n0.0 %f xlo xhi\n0.0 %f ylo yhi\n0.0 %f zlo zhi\n\nMasses\n\n1 1\n2 1\n3 1\n4 1\n\nAtoms\n\n"%(natoms,nbonds,nangles,box,box,box))

	atom_id=0
	mol_id=0
	bonds=[]
	angles=[]

	#Print free polymers
	for poly_id in range(npoly_free):
		atom_id+=1
		mol_id+=1
		r_old=rand(3)*box
		fout.write("%d %d 3 %f %.15f %.15f %.15f\n"%(atom_id,mol_id,poly_charge_free,r_old[0],r_old[1],r_old[2]))
		for n in range(1,polymer_length_free):
			atom_id+=1
			vec=uniform(-1,1,3)
			vec_norm=vec/sqrt(dot(vec,vec))
			r=r_old+bond_length*vec_norm
			x=r[0]
			y=r[1]
			z=r[2]
			bonds.append((atom_id,atom_id-1))
			if(n>1):
				angles.append((atom_id-2,atom_id-1,atom_id))
			fout.write("%d %d 3 %f %.15f %.15f %.15f\n"%(atom_id,mol_id,poly_charge_free,x,y,z))
			r_old=r

	#Print free polymers (of same type as the grafted ones)
	for poly_id in range(npoly_free_graft):
		atom_id+=1
		mol_id+=1
		r_old=rand(3)*box
		fout.write("%d %d 2 %f %.15f %.15f %.15f\n"%(atom_id,mol_id,poly_charge_free,r_old[0],r_old[1],r_old[2]))
		for n in range(1,polymer_length):
			atom_id+=1
			vec=uniform(-1,1,3)
			vec_norm=vec/sqrt(dot(vec,vec))
			r=r_old+bond_length*vec_norm
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
			r_old=r

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









