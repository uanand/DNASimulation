import numpy
import cPickle as pickle
import matplotlib.pyplot as plt
from matplotlib import cm
from mpi4py import MPI
from DNAclass import DNA

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

dList = numpy.arange(0.5,5.1,0.5)
BList = numpy.zeros(numpy.size(dList),dtype='int'); BList[:] = 50
MList = numpy.zeros(numpy.size(dList),dtype='int'); MList[:] = 200
deltaList = range(0,181,20)
numPerturb = 1e3

counter = 0
outFile = open('acceptProb_'+str(rank).zfill(2)+'.txt', 'wb')

for mode in ['2d','3d']:
	if (mode == '2d'):
		modeInt = 2
	else:
		modeInt = 3
	#probDict = {}
	for B,M,d in zip(BList,MList,dList):
		#probDict[d] = []
		for delta in deltaList:
			if (counter%size == rank):
				
				if (delta>0):
					hp = DNA(B,M,d,mode=mode)
					while (hp.feasibleCounter < numPerturb):
						hp.copyDNA()
						hp.perturb(perturbMode=numpy.random.randint(2), delta=delta)
						hp.acceptReject()
					acceptProb = 1.0*hp.acceptCounter/hp.feasibleCounter
					#probDict[d].append(1.0*hp.acceptCounter/hp.feasibleCounter)
					del hp
				else:
					#probDict[d].append(1.0)
					acceptProb = 1.0
				print modeInt, B, M, d, delta, acceptProb
				outFile.write("%d %d %d %f %d %f\n" %(modeInt, B, M, d, delta, acceptProb))
				#outFile.write("Hello world\n")
			counter += 1
	#probDict['deltaList'] = deltaList
	#pickle.dump(probDict, open('probDict_'+mode, 'wb'))
	
outFile.close()

##################################################
# PLOTTING THE PROBABILTY OF ACCEPTANCE FOR DIFFERENT
# d VALUES AND DELTA
##################################################
#for mode in ['2d','3d']:
	#probDict = pickle.load(open('probDict','rb'))
	#keys = probDict.keys(); keys.remove('deltaList')
	#keys.sort()
	
	#start,stop,numLines = 0.0,1.0,len(keys)
	#cm_subsection = numpy.linspace(start,stop,numLines)
	#colors = [cm.jet(x) for x in cm_subsection]
	#cbarImg = numpy.zeros([len(colors),2])
	
	#fig = plt.figure(figsize=(3,2.5))
	#ax = fig.add_axes([0,0,1,1])
	#minKey,maxKey = min(keys),max(keys)
	#for key,color,i in zip(keys,colors,range(len(keys))):
		#ax.plot(probDict['deltaList'],probDict[key],color=color)
		#cbarImg[i,:] = len(keys)-1-i
	#ax.axhline(y=0.5)#,xmin=min(probDict['deltaList']),xmax=max(probDict['deltaList']))
	#ax.set_xlabel(r'$\alpha$ (deg)')
	#ax.set_ylabel('Acceptance Probability')
	##ax.text(100,0.9,r'Perturbation ~ U(-$\delta$,$\delta$)')
	#ax.set_xlim(min(probDict['deltaList']),max(probDict['deltaList']))
	#ax.set_ylim(0,1)
	
	#ax2 = fig.add_axes([1.1,0.75,0.04,0.15])
	#ax2.imshow(cbarImg, cmap='jet')
	#ax2.set_xticks([]); ax2.set_yticks([])
	#ax2.set_xlabel('d='+str(keys[0]))
	#ax2.set_title('d='+str(keys[-1]))
	#plt.savefig('Probability_'+mode+'.png',format='png')
	#plt.savefig('Probability_'+mode+'.pdf',format='pdf')
	#plt.show()
##################################################
