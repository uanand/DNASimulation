import numpy
import cPickle as pickle
import matplotlib.pyplot as plt
import os
from matplotlib import cm
from mpi4py import MPI
from DNAclass import DNA

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

dList = numpy.arange(0.5,5.1,0.1)
BList = numpy.zeros(numpy.size(dList),dtype='int'); BList[:] = 50
MList = numpy.zeros(numpy.size(dList),dtype='int'); MList[:] = 200
deltaList = range(0,181,1)
numPerturb = 1e5

counter = 0
outFile = open('acceptProb_'+str(rank).zfill(2)+'.txt', 'wb')

for mode in ['2d','3d']:
	if (mode == '2d'):
		modeInt = 2
	else:
		modeInt = 3
	for B,M,d in zip(BList,MList,dList):
		if (rank == 0):
			print mode, B,M,d
		for delta in deltaList:
			if (counter%size == rank):
				if (delta>0):
					hp = DNA(B,M,d,mode=mode)
					while (hp.feasibleCounter < numPerturb):
						hp.copyDNA()
						hp.perturb(perturbMode=numpy.random.randint(2), delta=delta)
						hp.acceptReject()
					acceptProb = 1.0*hp.acceptCounter/hp.feasibleCounter
					del hp
				else:
					acceptProb = 1.0
				outFile.write("%d %d %d %.1f %d %f %d\n" %(modeInt, B, M, d, delta, acceptProb, numPerturb))
			counter += 1
outFile.close()

comm.Barrier()
if (rank == 0):
	for i in range(size):
		data = numpy.loadtxt('acceptProb_'+str(i).zfill(2)+'.txt', delimiter=' ', skiprows=0)
		print data.shape
		if (i == 0):
			DATA = data.copy()
		else:
			DATA = numpy.concatenate((DATA,data))
		os.remove('acceptProb_'+str(i).zfill(2)+'.txt')
	print DATA.shape
	sort = numpy.lexsort((DATA[:,4],DATA[:,3],DATA[:,2],DATA[:,1],DATA[:,0]))
	
	outFile = open('acceptProb.txt', 'wb')
	outFile.write("2D/3D B M d delta AcceptanceProbability NumberOfPerturbations\n")
	for i in sort:
		outFile.write("%d %d %d %.1f %d %f %d\n" %(DATA[i,0],DATA[i,1],DATA[i,2],DATA[i,3],DATA[i,4],DATA[i,5],DATA[i,6]))
	outFile.close()
	
##################################################
# PLOTTING THE PROBABILTY OF ACCEPTANCE FOR DIFFERENT
# d VALUES AND DELTA
##################################################
dList = numpy.arange(0.5,5.1,0.5)
BList = numpy.zeros(numpy.size(dList),dtype='int'); BList[:] = 50
MList = numpy.zeros(numpy.size(dList),dtype='int'); MList[:] = 200

comm.Barrier()
if (rank == 0):
	DATA = numpy.loadtxt('acceptProb.txt', delimiter=' ', skiprows=1)
	
	start,stop,numLines = 0.0,1.0,dList.size
	cm_subsection = numpy.linspace(start,stop,numLines)
	colors = [cm.jet(x) for x in cm_subsection]
	cbarImg = numpy.zeros([len(colors),2])
	
	for mode in ['2d','3d']:
		fig = plt.figure(figsize=(2.5,2))
		ax = fig.add_axes([0,0,1,1])
		if (mode == '2d'):
			modeInt = 2
		else:
			modeInt = 3
		for B,M,d,i,color in zip(BList,MList,dList,range(dList.size),colors):
			data = DATA[DATA[:,0] == modeInt]
			data = data[data[:,1] == B]
			data = data[data[:,2] == M]
			data = data[data[:,3] == d]
			x,y = data[:,4],data[:,5]
			ax.plot(x,y,color=color)
			cbarImg[i,:] = dList.size-1-i
		
		ax.axhline(y=0.5)
		ax.set_xlabel(r'$\alpha$ (deg)')
		ax.set_ylabel('Acceptance Probability')
		ax.set_xlim(numpy.min(x),numpy.max(x))
		ax.set_ylim(0,1)
		
		ax2 = fig.add_axes([1.1,0.75,0.04,0.15])
		ax2.imshow(cbarImg, cmap='jet')
		ax2.set_xticks([]); ax2.set_yticks([])
		ax2.set_xlabel('d='+str(dList[0]))
		ax2.set_title('d='+str(dList[-1]))
		plt.savefig('Probability_'+mode+'.png',format='png')
		plt.savefig('Probability_'+mode+'.pdf',format='pdf')
		plt.close()
##################################################
