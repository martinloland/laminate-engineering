# -*- coding: utf-8 -*-
"""
plotting.py
@author: Martin Loeland
@date: 06.03.17
"""
import numpy as np
import matplotlib.pyplot as plt
import laminates as lam

def reactionForce(laminate, direction, filename=None, show=False):
	deltaStrain = np.zeros(6)
	deltaStrain[direction] = 0.0001
	forceVec = np.zeros(6)
	force = []
	strain = []
	
	while True:
		#strain from force e=S*N
		strainVec = lam.strainFromForce(laminate, forceVec)
		
		#Check for failure
		failDir = lam.checkForFailure(laminate, strainVec)
		if  failDir == 1:
			# Ultimate fiber failure
			break
		
		# Increment strain
		strainVec += deltaStrain
		strain.append(strainVec[direction])
		
		#Find forces
		force.append(lam.reactionForce(laminate, strainVec)[direction])
		forceVec[direction] = force[-1]		

	## Plot here
	if direction == 0:
		cor = "x"
	elif direction == 1:
		cor = "y"
	elif direction == 2:
		cor = "xy"
	
	plt.figure()	
	plt.plot(strain, force)

	plt.xlabel('Strain $e_{'+cor+'}^0$')
	plt.ylabel('Force $N_{'+cor+'}$')
	plt.title('Laminate load with degradation')
	plt.grid(True)
	
	if filename: # Save figure
		plt.savefig(filename+".png")
		plt.savefig(filename+".pdf")
		plt.savefig(filename+'.pgf')		
	if show:
		plt.show()	