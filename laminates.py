# -*- coding: utf-8 -*-
"""
laminates.py
@author: Martin Loeland
@date: 04.03.17
"""
import numpy as np
import plotting as plot
from numpy import cos as cos
from numpy import sin as sin
import csv

class Ply:
	def __init__(self, E1, E2, G12, v12, thi, ang):
		self.E1 = E1
		self.E2 = E2
		self.G12 = G12
		self.v12 = v12
		self.v21 = v12*E2/E1
		self.t = thi
		self.ang = ang
		self.degraded = False
		self.nb = None
		
		self.stiffness()
	
	def stiffness(self):		
		self.S = self.findS() #Compliance in ply 
		self.Q = np.linalg.inv(self.S) #Stiffness in ply 
		self.Qk = self.rotateQ(self.Q, self.ang) #Stiffness in laminate
		
	def findS(self):
		S = np.zeros((3,3))
		S[0,0] = 1.0/self.E1
		S[0,1] = -self.v12/self.E1
		S[1,0] = -self.v12/self.E1
		S[1,1] = 1.0/self.E2
		S[2,2] = 1.0/self.G12
		return S
	
	def rotateQ(self, Q, ang):
		# Qk = TsInv*Q*Te
		return np.linalg.inv(Ts(ang)).dot(Q).dot(Te(ang))
		
	def setStrength(self,xt,xc,yt,yc,s12):
		self.xt = xt
		self.xc = xc
		self.yt = yt
		self.yc = yc
		self.s12 = s12
		
	def setDegradingFactors(self, E1d, E2d, G12d, v12d):
		self.E1d = E1d
		self.E2d = E2d
		self.G12d = G12d
		self.v12d = v12d
		
	def degrade(self):
		# Degrade material properties
		self.E1 *= self.E1d
		self.E2 *= self.E2d
		self.G12 *= self.G12d
		self.v12 *= self.v12d
		self.v21 = self.v12*self.E2/self.E1
		
		# Calculate new stiffness
		self.stiffness()
		
		self.degraded = True
		
class Laminate:
	def __init__(self, plyList, half=False, odd=False):
		if half:
			# Mirror the laminate to produce whole laminate
			self.plys = self.mirrorLaminate(plyList, odd)
		else:
			self.plys = plyList
			
		#Declare which number in the stack each ply is
		for i in range(len(plyList)):
			plyList[i].nb = i
		
		# Laminate thickness
		self.t = self.totalThi(self.plys)

		self.stiffness()
		self.fe = None #Exposure factor
	
	def stiffness(self):
		# Find submatrixes
		self.A = self.findSubmat('A')
		self.B = self.findSubmat('B')
		self.D = self.findSubmat('D')
		# Stiffness and compliance
		self.Q = self.findQ()
		self.S = np.linalg.inv(self.Q)
		
	def __str__(self):
		str = "["+repr(self.plys[0].ang)
		for ply in self.plys[1:]:
			str += ("/"+repr(ply.ang))
		str += "]"
		return str
	

	def totalThi(self, plyList):
		t = 0
		for ply in plyList:
			t += ply.t
		return t
		
	def mirrorLaminate(self, plyList, odd):
		revList = plyList[::-1]
		if odd: #odd number, dont copy middle element
			revList = revList[1:]
		return plyList + revList
	
	def findSubmat(self,case):
		X = np.zeros((3,3))
		h_sum = 0
		for ply in self.plys:
			h_sum += ply.t
			hk = -self.t/2+h_sum
			hj = hk - ply.t
			if case == 'A':
				X += ply.Qk * (hk - hj)
			elif case == 'B':
				X += 0.5 * ply.Qk * (hk**2 - hj**2)
			elif case == 'D':
				X += 1.0/3.0 * ply.Qk * (hk**3 - hj**3)
		return X
	
	def findQ(self):
		Q = np.zeros((6,6))
		Q[0:3,0:3]=self.A
		Q[0:3,3:7]=self.B
		Q[3:7,0:3]=self.B
		Q[3:7,3:7]=self.D
		return Q

class Loadcase():
	def __init__(self, laminate, load, maxStrain=None):
		self.l = laminate #Laminate class
		self.f = np.asarray(load)
		self.mStrain = np.asarray(maxStrain) #max strain
		self.strain = np.dot(laminate.S,self.f) #calculated strain
		
	def validLoad(self):
		self.fe = 0 #Strain exposure factor
		valid = True
		# Check that the current load does not exceed the requirement
		for requirement, strain in zip(self.mStrain,self.strain):
			if requirement: #Only check for non-zero values
				if strain >= requirement:
					valid = False
				if strain/requirement > self.fe:
					self.fe = strain/requirement
		return valid
		
def checkForFailure(laminate, strainVec):
	strainVec = strainVec[:3] # Only include normal directions
	ff = 0 #Start with no failure
	ultimateFail = False
	h_sum = 0
	kappa = curvature(laminate, strainVec)
	for ply in laminate.plys:
		# h-coordinates (z)
		h_sum += ply.t
		hk = -laminate.t/2+h_sum
		hj = hk - ply.t
		
		# Ply strain, top and bottom, lam direction, e=e0+zk
		e_ply_top = strainVec + hk*kappa
		e_ply_bot = strainVec + hj*kappa
		
		# Ply stress, top and bottom, lam direction, s = Q*e
		s_ply_top = ply.Qk.dot(e_ply_top)
		s_ply_bot = ply.Qk.dot(e_ply_bot)
		
		# Ply stress, top and bottom, fiber direction, s' = Ts*s 
		s_ply_top_f = Ts(ply.ang).dot(s_ply_top)
		s_ply_bot_f = Ts(ply.ang).dot(s_ply_bot)
		
		# Check if stress is above strength of ply
		failureDir = aboveStrength(s_ply_top_f, s_ply_bot_f, ply)
		if failureDir and not ply.degraded:			
			ply.degrade() # Reduce strength of this ply
		if failureDir == 1 and not ultimateFail:
			ultimateFail = True
		elif failureDir != 0:
			ff = failureDir
	if ff:
		# Update the reduced stiffness of the laminate
		laminate.stiffness()
		
	if ultimateFail:
		return 1
	else:
		return ff # in which direction failure happened
	
def curvature(laminate, strainVec):
	'''
	FIX this if support for curvature is wanted
	when calculating strain in plies.
	Has no effect if no moment is applied and the
	laminate is symmetric and balanced
	'''
	return np.zeros(3)
	
def reactionForce(laminate, strainVec):
	# N = Q*e
	return laminate.Q.dot(strainVec)
	
def strainFromForce(laminate, forceVec):
	# e=S*N
	return laminate.S.dot(forceVec)
	
def aboveStrength(sTop, sBot, ply):
	# Using max stress criteria
	failureDir = 0
	exposure1dir = max(sTop[0]/ply.xt, sBot[0]/ply.xt,-sTop[0]/ply.xc, -sBot[0]/ply.xc)
	exposure2dir = max(sTop[1]/ply.yt, sBot[1]/ply.yt,-sTop[1]/ply.yc, -sBot[1]/ply.yc)
	exposure12dir = max(abs(sTop[2])/ply.s12,abs(sBot[2])/ply.s12)

	if exposure1dir >= 1:
		failureDir = 1
	elif exposure2dir >= 1:
		failureDir = 2
	elif exposure12dir >= 1:
		failureDir = 3
	
	return failureDir

def Ts(theta):
	a = np.radians(theta)
	return np.array([[cos(a)**2, sin(a)**2, 2*cos(a)*sin(a)],
					 [sin(a)**2, cos(a)**2, -2*cos(a)*sin(a)],
					 [-cos(a)*sin(a), cos(a)*sin(a), cos(a)**2-sin(a)**2]])

def Te(theta):
	a = np.radians(theta)
	return np.array([[cos(a)**2, sin(a)**2, cos(a)*sin(a)],
					 [sin(a)**2, cos(a)**2, -cos(a)*sin(a)],
					 [-2*cos(a)*sin(a), 2*cos(a)*sin(a), cos(a)**2-sin(a)**2]])
					 
def plyListFromCSV(filename):
	plyList = []
	with open(filename, 'rb') as f:
		reader = csv.DictReader(f,delimiter = ',')
		for row in reader:
			plyList.append(Ply(
			float(row['E1']),float(row['E2']),float(row['G12']),
			float(row['v12']),float(row['thi']),float(row['ang'])))
	f.close()
	return plyList
	
def plotReactionForce(laminate, direction, filename=None, show=False):
	plot.reactionForce(laminate, direction, filename, show)