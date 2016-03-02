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
B = 50
M = 200

numDNA = 1e5
numPerturb = 1e5
dList = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
deltaList = [28,35,39,45,49,53,57,60,64,67]
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
for d, delta in zip(dList,deltaList):
	if (rank == 0):
		mkdir(str(d)+'_'+str(delta))
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
			
			DNAdict['x'] = hp.x0
			DNAdict['y'] = hp.y0
			DNAdict['z'] = hp.z0
			DNAdict['acceptCounter'] = hp.acceptCounter
			DNAdict['feasibleCounter'] = hp.feasibleCounter
			DNAdict['totalCounter'] = hp.totalCounter
			DNAdict['correlationDict'] = hp.corrDict
			DNAdict['bendingAngleList'] = hp.bendingAngleList
			DNAdict['end2end'] = hp.end2end
			DNAdict['B'] = hp.B
			DNAdict['M'] = hp.M
			DNAdict['d'] = hp.d
			DNAdict['time'] = toc-tic
			pickle.dump(DNAdict, open(str(d)+'_'+str(delta)+'/'+str(DNAcounter).zfill(6), 'wb'))
			del hp, DNAdict
			
			if (rank == 0):
				print d,delta,DNAcounter,toc-tic
####################################################
####################################################
