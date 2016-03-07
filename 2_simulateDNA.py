import os
import numpy
import time
import cPickle as pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpi4py import MPI
from DNAclass import DNA

def mkdir(dirName):
	if (os.path.exists(dirName) == False):
		os.makedirs(dirName)
		
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

####################################################
####################################################
# USER INPUTS
numDNA = 10
numPerturb = 1e5
dList = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
deltaList = [28,35,39,45,49,53,57,60,64,67]
BList = [50,50,50,50,50,50,50,50,50,50]
MList = [200,200,200,200,200,200,200,200,200,200]
####################################################
####################################################


####################################################
####################################################
#SIMULATING DNA BASED ON THE DIFFERENT d AND delta VALUES
#AT THE LAST STAGE THESE VALUES WILL BE STORED IN A DICTIONARY
#1. X,Y,Z COORDINATES OF THE DNA
#2. TOTAL NUMBER OF PERTURBATION TRIALS, FEASIBLE TRIALS, AND ACCEPTED TRIALS
#3. TANGENT VECTOR CORRELATION
#4. BENDING ANGLE LIST
#5. END TO END DISTANCE
#6. TIME REQUIRED FOR THE SIMULATION AND COMPUTATIONS
####################################################
####################################################
for B,M,d,delta in zip(BList,MList,dList,deltaList):
	if (rank == 0):
		mkdir(str(B)+'_'+str(M)+'_'+str(d)+'_'+str(delta))
		
for B,M,d,delta in zip(BList,MList,dList,deltaList):
	for DNAcounter in range(1,int(numDNA)+1):
		if ((DNAcounter-1)%size == rank):
			DNAdict = {}
			tic = time.time()
			hp = DNA(B,M,d)
			while (hp.feasibleCounter < numPerturb):
				hp.copyDNA()
				hp.perturb(mode=numpy.random.randint(2), delta=delta)
				hp.acceptReject()
			hp.tangentCorrelation()
			hp.bendingAngle()
			hp.end2endDistance()
			toc = time.time()
			
			DNAdict['x'] = hp.x0.astype('float32')
			DNAdict['y'] = hp.y0.astype('float32')
			DNAdict['z'] = hp.z0.astype('float32')
			DNAdict['acceptCounter'] = hp.acceptCounter
			DNAdict['feasibleCounter'] = hp.feasibleCounter
			DNAdict['totalCounter'] = hp.totalCounter
			DNAdict['tangentCorrList'] = hp.corrList.astype('float32')
			DNAdict['bendingAngleList'] = hp.bendingAngleList.astype('float32')
			DNAdict['end2end'] = hp.end2end
			DNAdict['B'] = hp.B
			DNAdict['M'] = hp.M
			DNAdict['d'] = hp.d
			DNAdict['time'] = toc-tic
			pickle.dump(DNAdict, open(str(B)+'_'+str(M)+'_'+str(d)+'_'+str(delta)+'/'+str(DNAcounter).zfill(len(str(int(numDNA)))), 'wb'))
			del hp, DNAdict
			if (rank == 0):
				print d,delta,DNAcounter,toc-tic
####################################################
####################################################


####################################################
####################################################
# COMBINE ALL DICTIONARIES INTO SINGLE HDF5 FILE
####################################################
####################################################
if (rank == 0):
	for B,M,d,delta in zip(BList,MList,dList,deltaList):
		h5 = h5py.File(str(B)+'_'+str(M)+'_'+str(d)+'_'+str(delta)+'/DNA.h5py', 'w')
		for i in range(1,numDNA+1):
			fileName = str(B)+'_'+str(M)+'_'+str(d)+'_'+str(delta)+'/'+str(i).zfill(len(str(int(numDNA))))
			DNAdict = pickle.load(open(fileName,'rb'))
			DNA = h5.create_group(str(i).zfill(len(str(int(numDNA)))))
			DNA.create_dataset('x', data=DNAdict['x'], compression='gzip')
			DNA.create_dataset('y', data=DNAdict['y'], compression='gzip')
			DNA.create_dataset('z', data=DNAdict['z'], compression='gzip')
			DNA.create_dataset('acceptCounter', data=DNAdict['acceptCounter'])
			DNA.create_dataset('feasibleCounter', data=DNAdict['feasibleCounter'])
			DNA.create_dataset('totalCounter', data=DNAdict['totalCounter'])
			DNA.create_dataset('tangentCorrList', data=DNAdict['tangentCorrList'], compression='gzip')
			DNA.create_dataset('bendingAngleList', data=DNAdict['bendingAngleList'], compression='gzip')
			DNA.create_dataset('end2end', data=DNAdict['end2end'])
			DNA.create_dataset('B', data=DNAdict['B'])
			DNA.create_dataset('M', data=DNAdict['M'])
			DNA.create_dataset('d', data=DNAdict['d'])
			DNA.create_dataset('time', data=DNAdict['time'])
			del DNAdict
			os.remove(fileName)
		h5.close()
####################################################
####################################################
