import os
import numpy
import time
import cPickle as pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import optimize
from DNAclass import DNA

####################################################
####################################################
# USER INPUTS
numDNA = 1e5
numPerturb = 1e5
dList = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
deltaList = [28,35,39,45,49,53,57,60,64,67]
BList = [50,50,50,50,50,50,50,50,50,50]
MList = [200,200,200,200,200,200,200,200,200,200]
####################################################
####################################################


def tangentCorr(x,lamda):
	f = numpy.exp(-x/lamda)
	return f

####################################################
####################################################
#FINDING THE TANGENT CORRELATION USING THE DNA DICTIONARIES
#ABOVE. THE CORRELATION VALUES WILL BE AVERAGED OVER ALL THE
#DIFFERENT DNA.
####################################################
####################################################
for B,M,d,delta in zip(BList,MList,dList,deltaList):
	h5 = h5py.File(str(B)+'_'+str(M)+'_'+str(d)+'_'+str(delta)+'/DNA.h5py', 'r')
	x,corrList,thetaList,end2endList,corrDict,analysisDict = range(M),[],[],[],{},{}
	for key in h5.keys():
		if (key == h5.keys()[0]):
			corrArr = h5[key]['tangentCorrList'].value
		else:
			corrArr += h5[key]['tangentCorrList'].value
		for theta in h5[key]['bendingAngleList'].value:
			thetaList.append(theta)
		end2endList.append(h5[key]['end2end'].value)
	
	#DNAdict = pickle.load(open(str(dList[0])+'_'+str(deltaList[0])+'/1','rb'))
	#for key in DNAdict['correlationDict'].keys():
		#x.append(key)
		#corrDict[key] = []
	#for i in range(1,numDNA+1):
		#print i,B,M,d,delta
		#DNAdict = pickle.load(open(str(d)+'_'+str(delta)+'/'+str(i),'rb'))
		#for key in DNAdict['correlationDict'].keys():
			#corrDict[key].append(numpy.mean(DNAdict['correlationDict'][key]))
		#for theta in DNAdict['bendingAngleList']:
			#thetaList.append(theta)
		#end2endList.append(DNAdict['end2end'])
	#for key in DNAdict['correlationDict'].keys():
		#corrList.append(numpy.mean(corrDict[key]))
	analysisDict['x'] = x
	analysisDict['correlation'] = corrArr/len(h5.keys())
	analysisDict['theta'] = thetaList
	analysisDict['end2end'] = end2endList
	pickle.dump(analysisDict, open(str(B)+'_'+str(M)+'_'+str(d)+'_'+str(delta)+'/analysisDict', 'wb'))
	h5.close()
####################################################
####################################################
	
	
####################################################
####################################################
#PLOTTING THE FIGURES - VECTOR CORRELATION	
fig = plt.figure(figsize=(2.5,1.5))
ax = fig.add_axes([0,0,1,1])
for B,M,d,delta in zip(BList,MList,dList,deltaList):
	analysisDict = pickle.load(open(str(d)+'_'+str(delta)+'/analysisDict','rb'))
	
	#INITIAL GUESS
	lamda = 1
	#FITTING WITH EXPONENTIAL DECAY FUNCTION
	flag = True
	[params, pcov] = optimize.curve_fit(tangentCorr, analysisDict['x'], analysisDict['correlation'], [lamda])
	lamda = params[0]*d
	
	X = numpy.linspace(0,M-1,1000)
	Y = numpy.exp(-X/params[0])
	print params[0], params[0]*d
	ax.plot(analysisDict['x'],analysisDict['correlation'],label='Correlation',lw=2)
	ax.plot(X,Y,label='Fit')
	ax.set_xlabel(r'$\Delta$(L)')
	ax.set_ylabel('Correlation')
	ax.set_xlim(0,M)
	ax.set_ylim(0,1)
	plt.savefig(str(d)+'_'+str(delta)+'/figCorr.png',format='png')
	plt.savefig(str(d)+'_'+str(delta)+'/figCorr.pdf',format='pdf')
	plt.cla()
	
	ax.hist(analysisDict['theta'], bins=range(0,21,1), normed=True)
	ax.set_xlabel('Bending angle (deg)')
	ax.set_ylabel('Probability')
	plt.savefig(str(d)+'_'+str(delta)+'/figBendingAngle.png',format='png')
	plt.savefig(str(d)+'_'+str(delta)+'/figBendingAngle.pdf',format='pdf')
	plt.cla()
	
	print numpy.mean(numpy.asarray(analysisDict['end2end'])**2)
	print numpy.var(analysisDict['end2end'])
	print 2*lamda*100 - 2*lamda**2*(1-numpy.exp(-100/lamda))
	
plt.close()
####################################################
####################################################


#####################################################
#####################################################
##CALCULATING THE BENDING ANGLE DISTRIBUTION FOR ALL TEH DNA
##THE CORRESPONDING FIT FUNCTION WILL BE DONE LATER
#####################################################
#####################################################
##y = []
##for i in range(1,numDNA+1):
	##print i
	##DNAdict = pickle.load(open('0.5_28/'+str(i),'rb'))
	##for theta in DNAdict['bendingAngleList']:
		##y.append(theta)
##plt.figure()
##plt.hist(y, bins=range(0,41,5), normed=True)
##plt.show()
#####################################################
#####################################################


#####################################################
#####################################################	
##CALCULATING THE END TO END DISTANCE DISTRIBUTION FOR ALL THE DNA.
##THE CORRESPONDING FIT FUNCTION WILL BE DONE LATER
#####################################################
#####################################################
#y = []
#for i in range(1,numDNA+1):
	#print i
	#DNAdict = pickle.load(open('0.5_28//'+str(i),'rb'))
	#y.append(DNAdict['end2end'])
#plt.figure()
#plt.hist(y,bins=50)
#plt.show()
#####################################################
#####################################################


#####################################################
#####################################################
##PLOTTING ALL THE FINAL DNA TOGETHER AFTER GOING
##THROUGH numPerturb PERTURBATIONS
#####################################################
#####################################################
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#for i in range(1,numDNA+1):
	#DNAdict = pickle.load(open('1.0_30/'+str(i),'rb'))
	#ax.plot(DNAdict['x'],DNAdict['y'],DNAdict['z'],color='k')
	##ax.set_xlabel('X'), ax.set_ylabel('Y'), ax.set_zlabel('Z')
#ax.plot([0],[0],[0],marker='o',markersize=3,color='r')
#ax.set_axis_off()
#plt.show()
#####################################################
#####################################################
