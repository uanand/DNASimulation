import numpy
from numpy import sin, cos, arccos, rad2deg, deg2rad, sqrt

class DNA(object):
	##############################################################
	def __init__(self,B,M,d):
		self.B = B
		self.M = M
		self.d = d
		self.bCoeff = 1.0*self.B/self.d
		self.totalCounter = 0
		self.feasibleCounter = 0
		self.acceptCounter = 0
		self.x0,self.y0,self.z0 = numpy.linspace(0,self.d*(self.M-1),self.M),numpy.zeros(self.M),numpy.zeros(self.M)
		self.e0 = 0
		self.eps = 1e-20
		
	def rotate(self,x,y,z,a,b,c,u,v,w,alpha):
		x,y,z = (a*(v**2+w**2) - u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(alpha)) + x*cos(alpha) + (-c*v+b*w-w*y+v*z)*sin(alpha), (b*(u**2+w**2) - v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(alpha)) + y*cos(alpha) + (c*u-a*w+w*x-u*z)*sin(alpha), (c*(u**2+v**2) - w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(alpha)) + z*cos(alpha) + (-b*u+a*v-v*x+u*y)*sin(alpha)
		return x,y,z
		
	def normalize(self,x,y,z):
		length = sqrt(x**2+y**2+z**2)
		return x/length,y/length,z/length
		
	def copyDNA(self):
		self.x1,self.y1,self.z1 = self.x0.copy(),self.y0.copy(),self.z0.copy()
		self.e1 = self.e0
		
	def compareDNA(self):
		if (numpy.max(numpy.abs(numpy.concatenate((self.x1,self.y1,self.z1)) - numpy.concatenate((self.x0,self.y0,self.z0)))) <= self.eps):
			return 0
		else:
			return 1
		
	def updateDNA(self):
		self.x0,self.y0,self.z0 = self.x1.copy(),self.y1.copy(),self.z1.copy()
		self.e0 = self.e1.copy()
		
	def energy(self):
		lengthSegment = sqrt((self.x0[1:]-self.x0[:-1])**2 + (self.y0[1:]-self.y0[:-1])**2 + (self.z0[1:]-self.z0[:-1])**2)
		tangentX, tangentY, tangentZ = (self.x0[1:]-self.x0[:-1])/lengthSegment, (self.y0[1:]-self.y0[:-1])/lengthSegment, (self.z0[1:]-self.z0[:-1])/lengthSegment
		self.e0 = self.bCoeff/2*numpy.sum((tangentX[1:]-tangentX[:-1])**2 + (tangentY[1:]-tangentY[:-1])**2 + (tangentZ[1:]-tangentZ[:-1])**2)
		
		lengthSegment = sqrt((self.x1[1:]-self.x1[:-1])**2 + (self.y1[1:]-self.y1[:-1])**2 + (self.z1[1:]-self.z1[:-1])**2)
		tangentX, tangentY, tangentZ = (self.x1[1:]-self.x1[:-1])/lengthSegment, (self.y1[1:]-self.y1[:-1])/lengthSegment, (self.z1[1:]-self.z1[:-1])/lengthSegment
		self.e1 = self.bCoeff/2*numpy.sum((tangentX[1:]-tangentX[:-1])**2 + (tangentY[1:]-tangentY[:-1])**2 + (tangentZ[1:]-tangentZ[:-1])**2)
		
	def perturb(self,mode=0,delta=20):
		if (mode==0):
			i = numpy.random.randint(1,self.M-1)
			theta = numpy.random.rand()*numpy.pi
			phi = numpy.random.rand()*numpy.pi*2
			alpha = numpy.deg2rad(numpy.random.rand()*delta*2-delta)
			a,b,c = self.x1[i],self.y1[i],self.z1[i]
			u,v,w = sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)
			self.x1[i+1:],self.y1[i+1:],self.z1[i+1:] = self.rotate(self.x1[i+1:],self.y1[i+1:],self.z1[i+1:],a,b,c,u,v,w,alpha)
			
		elif (mode==1):
			i = numpy.random.randint(0,self.M-2)
			j = numpy.random.randint(i+2,self.M)
			if (i==0 and j==self.M-1):
				if (numpy.random.randint(2) == 0):
					i+=1
				else:
					j-=1
			theta = numpy.random.rand()*numpy.pi
			phi = numpy.random.rand()*numpy.pi*2
			alpha = numpy.deg2rad(numpy.random.rand()*delta*2-delta)
			a,b,c = self.x1[i],self.y1[i],self.z1[i]
			u,v,w = self.normalize(self.x1[j]-self.x1[i],self.y1[j]-self.y1[i],self.z1[j]-self.z1[i])
			self.x1[i+1:j],self.y1[i+1:j],self.z1[i+1:j] = self.rotate(self.x1[i+1:j],self.y1[i+1:j],self.z1[i+1:j],a,b,c,u,v,w,alpha)
			
	def acceptReject(self):
		lengthSegment = sqrt((self.x1[1:]-self.x1[:-1])**2 + (self.y1[1:]-self.y1[:-1])**2 + (self.z1[1:]-self.z1[:-1])**2)
		tangentX, tangentY, tangentZ = (self.x1[1:]-self.x1[:-1])/lengthSegment, (self.y1[1:]-self.y1[:-1])/lengthSegment, (self.z1[1:]-self.z1[:-1])/lengthSegment
		self.e1 = self.bCoeff/2*numpy.sum((tangentX[1:]-tangentX[:-1])**2 + (tangentY[1:]-tangentY[:-1])**2 + (tangentZ[1:]-tangentZ[:-1])**2)
		
		if (self.compareDNA()==0):
			self.totalCounter+=1
		else:
			self.totalCounter+=1
			self.feasibleCounter+=1
			if (self.e1 < self.e0):
				self.updateDNA()
				self.acceptCounter+=1
			elif (numpy.random.rand() <= numpy.exp(-(self.e1-self.e0))):
				self.updateDNA()
				self.acceptCounter+=1
				
	def tangentCorrelation(self):
		self.corrDict, self.corrList = {}, []
		self.lengthSegment = sqrt((self.x0[1:]-self.x0[:-1])**2 + (self.y0[1:]-self.y0[:-1])**2 + (self.z0[1:]-self.z0[:-1])**2)
		self.tangentX, self.tangentY, self.tangentZ = (self.x0[1:]-self.x0[:-1])/self.lengthSegment, (self.y0[1:]-self.y0[:-1])/self.lengthSegment, (self.z0[1:]-self.z0[:-1])/self.lengthSegment
		for i in range(self.M-1):
			self.corrDict[i] = []
			for j,k in zip(range(self.M-1-i),range(i,self.M-1)):
				self.corrDict[i].append(self.tangentX[j]*self.tangentX[k] + self.tangentY[j]*self.tangentY[k] + self.tangentZ[j]*self.tangentZ[k])
		for i in range(self.M-1):
			self.corrList.append(numpy.mean(self.corrDict[i])
		self.corrList = numpy.asarray(self.corrList)
			#corrDict[key].append(numpy.mean(DNAdict['correlationDict'][key]))
			
	def bendingAngle(self):
		self.bendingAngleList = []
		self.lengthSegment = sqrt((self.x0[1:]-self.x0[:-1])**2 + (self.y0[1:]-self.y0[:-1])**2 + (self.z0[1:]-self.z0[:-1])**2)
		self.tangentX, self.tangentY, self.tangentZ = (self.x0[1:]-self.x0[:-1])/self.lengthSegment, (self.y0[1:]-self.y0[:-1])/self.lengthSegment, (self.z0[1:]-self.z0[:-1])/self.lengthSegment
		for i,j in zip(range(self.M-2),range(1,self.M-1)):
			self.bendingAngleList.append(rad2deg(arccos(self.tangentX[i]*self.tangentX[j] + self.tangentY[i]*self.tangentY[j] + self.tangentZ[i]*self.tangentZ[j])))
			
	def end2endDistance(self):
		self.end2end = sqrt((self.x0[0]-self.x0[-1])**2 + (self.y0[0]-self.y0[-1])**2 + (self.z0[0]-self.z0[-1])**2)
		
