import os
import numpy
import time
import h5py
import cPickle as pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import optimize
from DNAclass import DNA

####################################################
####################################################
# USER INPUTS
numDNA = 1e4
numPerturb = 1e5
dList = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
deltaList = [27,34,40,45,50,53,57,61,64,67]
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
#DATA = numpy.loadtxt('acceptProb.txt', delimiter=' ', skiprows=1)
#dList = numpy.unique(DATA[:,3])
#BList = numpy.zeros(numpy.size(dList),dtype='int'); BList[:] = 50
#MList = numpy.zeros(numpy.size(dList),dtype='int'); MList[:] = 200

#for mode in ['3d']:  #['2d','3d']:
    #print mode
    #if (mode == '2d'):
        #modeInt = 2
    #else:
        #modeInt = 3
        
    #deltaList = []
    #for B,M,d in zip(BList,MList,dList):
        #data = DATA[DATA[:,0] == modeInt]
        #data = data[data[:,1] == B]
        #data = data[data[:,2] == M]
        #data = data[data[:,3] == d]
        #index = numpy.argmin(numpy.abs(data[:,5]-0.5))
        #deltaList.append(int(data[index,4]))
            
    #for B,M,d,delta in zip(BList,MList,dList,deltaList):
        #print B,M,d,delta
        #h5 = h5py.File(mode+'/'+str(B)+'_'+str(M)+'_'+str(d)+'_'+str(delta)+'/DNA.h5py', 'r')
        #x,corrList,thetaList,end2endList,corrDict,analysisDict = range(M-1),[],[],[],{},{}
        #for key in h5.keys():
            #print mode,B,M,d,delta,key
            #if (key == h5.keys()[0]):
                #corrArr = h5[key]['tangentCorrList'].value
            #else:
                #corrArr += h5[key]['tangentCorrList'].value
            #for theta in h5[key]['bendingAngleList'].value:
                #thetaList.append(theta)
            #end2endList.append(h5[key]['end2end'].value)
        #analysisDict['x'] = x
        #analysisDict['correlation'] = corrArr/len(h5.keys())
        #analysisDict['theta'] = thetaList
        #analysisDict['end2end'] = end2endList
        #pickle.dump(analysisDict, open(mode+'/'+str(B)+'_'+str(M)+'_'+str(d)+'_'+str(delta)+'/analysisDict', 'wb'))
        #h5.close()
####################################################
####################################################


####################################################
####################################################
#PLOTTING THE FIGURES - TANGENT CORRELATION
#DATA = numpy.loadtxt('acceptProb.txt', delimiter=' ', skiprows=1)
#dList = numpy.unique(DATA[:,3])
#BList = numpy.zeros(numpy.size(dList),dtype='int'); BList[:] = 50
#MList = numpy.zeros(numpy.size(dList),dtype='int'); MList[:] = 200
#fig = plt.figure(figsize=(2.5,1.5))
#ax = fig.add_axes([0,0,1,1])
#for mode in ['3d']:
    #if (mode == '2d'):
        #modeInt = 2
    #else:
        #modeInt = 3
        
    #deltaList = []
    #for B,M,d in zip(BList,MList,dList):
        #data = DATA[DATA[:,0] == modeInt]
        #data = data[data[:,1] == B]
        #data = data[data[:,2] == M]
        #data = data[data[:,3] == d]
        #index = numpy.argmin(numpy.abs(data[:,5]-0.5))
        #deltaList.append(int(data[index,4]))
        
    #theoryend2endDistList,end2endDistList,DNALength = [],[],[]
    #for B,M,d,delta in zip(BList,MList,dList,deltaList):
        #analysisDict = pickle.load(open(mode+'/'+str(B)+'_'+str(M)+'_'+str(d)+'_'+str(delta)+'/analysisDict','rb'))
        
        ##INITIAL GUESS
        #lamda = 1
        
        ##FITTING WITH EXPONENTIAL DECAY FUNCTION
        #flag = True
        #[params, pcov] = optimize.curve_fit(tangentCorr, analysisDict['x'], analysisDict['correlation'], [lamda])
        #lamda = params[0]*d
        #X = numpy.linspace(0,M-1,1000)
        #Y = numpy.exp(-X/params[0])
        #print params[0], params[0]*d
        #ax.plot(analysisDict['x'],analysisDict['correlation'],label='Correlation',lw=2)
        #ax.plot(X,Y,label='Fit')
        #ax.set_xlabel(r'$\Delta$(L)')
        #ax.set_ylabel('Correlation')
        #ax.set_xlim(0,M)
        #ax.set_ylim(0,1)
        #plt.savefig(mode+'/'+str(B)+'_'+str(M)+'_'+str(d)+'_'+str(delta)+'/figCorr.png',format='png')
        #plt.savefig(mode+'/'+str(B)+'_'+str(M)+'_'+str(d)+'_'+str(delta)+'/figCorr.pdf',format='pdf')
        #plt.cla()
        
        #hist, bin_edges = numpy.histogram(analysisDict['theta'], bins=range(0,91,1), normed=True)
        #bin_center = bin_edges[:-1] + (bin_edges[1:]-bin_edges[:-1])/2.0
        #ax.plot(bin_center,hist)
        #ax.set_xlabel('Bending angle (deg)')
        #ax.set_ylabel('Probability')
        #plt.savefig(mode+'/'+str(B)+'_'+str(M)+'_'+str(d)+'_'+str(delta)+'/figBendingAngle.png',format='png')
        #plt.savefig(mode+'/'+str(B)+'_'+str(M)+'_'+str(d)+'_'+str(delta)+'/figBendingAngle.pdf',format='pdf')
        #plt.cla()
        
        #end2endDistList.append(numpy.mean(numpy.asarray(analysisDict['end2end'])**2))
        #theoryend2endDistList.append(2*lamda*M*d - 2*lamda**2*(1-numpy.exp(-M*d/lamda)))
        #DNALength.append(M*d)
    #ax.plot(DNALength,end2endDistList,lw=2)
    #ax.plot(DNALength,theoryend2endDistList,lw=1)
    #ax.set_xlabel('Chain length (nm)')
    #ax.set_ylabel(r'Mean of square ($nm^2$)')
    #plt.savefig(mode+'/end2endVar.png',format='png')
    #plt.savefig(mode+'/end2endVar.pdf',format='pdf')
    
    #plt.close()
####################################################
####################################################


#####################################################
#####################################################
##CALCULATING THE BENDING ANGLE DISTRIBUTION FOR ALL THE DNA
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
#####################################################    ##CALCULATING THE END TO END DISTANCE DISTRIBUTION FOR ALL THE DNA.
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
fig = plt.figure(figsize=(4,4))
ax = fig.gca(projection='3d')

dList = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
deltaList = [27,34,40,45,50,53,57,61,64,67]
BList = [50,50,50,50,50,50,50,50,50,50]
MList = [200,200,200,200,200,200,200,200,200,200]

B,M,d,delta = BList[1],MList[1],dList[1],deltaList[1]
mode='3d'

h5 = h5py.File(mode+'/'+str(B)+'_'+str(M)+'_'+str(d)+'_'+str(delta)+'/DNA.h5py', 'r')
for i in range(1,101):
    x,y,z = h5[str(i).zfill(5)]['x'].value,h5[str(i).zfill(5)]['z'].value,h5[str(i).zfill(5)]['z'].value
    ax.plot(x,y,z,color='k')
ax.plot([0],[0],[0],linestyle='none',marker='o',color='r',fillstyle='full', markersize=5, markeredgewidth=0)
ax.set_xlabel('X (nm)')
ax.set_ylabel('Y (nm)')
ax.set_zlabel('Z (nm)')
ax.set_xticks([-200,-100,0,100,200])
ax.set_yticks([-200,-100,0,100,200])
ax.set_zticks([-200,-100,0,100,200])
plt.savefig(mode+'/sampleDNA.png',format='png')
plt.savefig(mode+'/sampleDNA.pdf',format='pdf')
plt.close()
h5.close()
    
    
#for mode in ['3d']:  #['2d','3d']:
    #print mode
    #if (mode == '2d'):
        #modeInt = 2
    #else:
        #modeInt = 3
        
    #for B,M,d in zip(BList,MList,dList):
        #data = DATA[DATA[:,0] == modeInt]
        #data = data[data[:,1] == B]
        #data = data[data[:,2] == M]
        #data = data[data[:,3] == d]
        #index = numpy.argmin(numpy.abs(data[:,5]-0.5))
        #deltaList.append(int(data[index,4]))
            
    #for B,M,d,delta in zip(BList,MList,dList,deltaList):
        #print B,M,d,delta
        #h5 = h5py.File(mode+'/'+str(B)+'_'+str(M)+'_'+str(d)+'_'+str(delta)+'/DNA.h5py', 'r')
        #x,corrList,thetaList,end2endList,corrDict,analysisDict = range(M-1),[],[],[],{},{}
        #for key in h5.keys():
            #print mode,B,M,d,delta,key
            #if (key == h5.keys()[0]):
                #corrArr = h5[key]['tangentCorrList'].value
            #else:
                #corrArr += h5[key]['tangentCorrList'].value
            #for theta in h5[key]['bendingAngleList'].value:
                #thetaList.append(theta)
            #end2endList.append(h5[key]['end2end'].value)
        #analysisDict['x'] = x
        #analysisDict['correlation'] = corrArr/len(h5.keys())
        #analysisDict['theta'] = thetaList
        #analysisDict['end2end'] = end2endList
        #pickle.dump(analysisDict, open(mode+'/'+str(B)+'_'+str(M)+'_'+str(d)+'_'+str(delta)+'/analysisDict', 'wb'))
        #h5.close()


#for i in range(1,20+1):
    #DNAdict = pickle.load(open('1.0_30/'+str(i),'rb'))
    #ax.plot(DNAdict['x'],DNAdict['y'],DNAdict['z'],color='k')
    ##ax.set_xlabel('X'), ax.set_ylabel('Y'), ax.set_zlabel('Z')
#ax.plot([0],[0],[0],marker='o',markersize=3,color='r')
#ax.set_axis_off()
#plt.show()
#####################################################
#####################################################
