import numpy
import time
import cPickle as pickle
import matplotlib.pyplot as plt
from matplotlib import cm
from DNAclass import DNA

BList = [50,50,50,50,50,50,50,50,50,50]
MList = [200,200,200,200,200,200,200,200,200,200]
dList = numpy.arange(0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0)
deltaList = range(0,181,2)
numPerturb = 1e5
probDict = {}

for B,M,d in zip(BList,MList,dList):
	probDict[d] = []
	for delta in deltaList:
		if (delta>0):
			tic = time.time()
			hp = DNA(B,M,d)
			while (hp.feasibleCounter < numPerturb):
				hp.copyDNA()
				hp.perturb(mode=numpy.random.randint(2), delta=delta)
				hp.acceptReject()
			probDict[d].append(1.0*hp.acceptCounter/hp.feasibleCounter)
			del hp
			toc = time.time()
			print d, delta, toc-tic
		else:
			probDict[d].append(1.0)
probDict['deltaList'] = deltaList
pickle.dump(probDict, open('probDict', 'wb'))

##################################################
# PLOTTING THE PROBABILTY OF ACCEPTANCE FOR DIFFERENT
# d VALUES AND DELTA
##################################################
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
#ax.set_xlabel(r'$\delta$ (deg)')
#ax.set_ylabel('Acceptance Probability')
#ax.text(100,0.9,r'Perturbation ~ U(-$\delta$,$\delta$)')
#ax.set_xlim(min(probDict['deltaList']),max(probDict['deltaList']))
#ax.set_ylim(0,1)

#ax2 = fig.add_axes([1.1,0.75,0.04,0.15])
#ax2.imshow(cbarImg, cmap='jet')
#ax2.set_xticks([]); ax2.set_yticks([])
#ax2.set_xlabel('d='+str(keys[0]))
#ax2.set_title('d='+str(keys[-1]))
#plt.savefig('Probability.png',format='png')
#plt.savefig('Probability.pdf',format='pdf')
#plt.show()
##################################################
