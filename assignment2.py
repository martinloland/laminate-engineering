# -*- coding: utf-8 -*-
"""
assignment2.py
@author: Martin Loeland
@date: 04.03.17
"""
import numpy as np
import copy
import laminates as lam
np.set_printoptions(precision=3)

def problem1():
	header(1)	
	E1 = 132000
	E2 = 8100
	G12 = 4700
	v12 = 0.27
	thi = 0.5
	
	# p,e: load and max strain vector respectively, [N/mm]	
	p1 = [300,0,0,0,0,0]
	p2 = [0,0,150,0,0,0]
	p3 = [0,0,0,500,0,0]
	p4 = [0,0,0,0,100,0]
	e1 = [0.001,0,0,0,0,0]
	e2 = [0,0,0.002,0,0,0]
	e3 = [0,0,0,0.0005,0,0]
	e4 = [0,0,0,0,0.0005,0]
	p = [p1, p2, p3, p4]
	e = [e1, e2, e3, e4]
	
	nbTested = 0
	validLayups = []
	angles = [0, 45, -45, 90]
	anglesBal = [0,90]
	
	# Test all possible layups
	for ang1 in angles:
		for ang2 in angles:
			for ang3 in angles:
				for ang4 in angles:
					for ang5 in angles:
						plyList = []
						plyList.append(lam.Ply(E1, E2, G12, v12, thi, ang1))  #6
						plyList.append(lam.Ply(E1, E2, G12, v12, thi, ang2))  #1
						plyList.append(lam.Ply(E1, E2, G12, v12, thi, ang3)) #2
						plyList.append(lam.Ply(E1, E2, G12, v12, thi, ang4))  #3
						plyList.append(lam.Ply(E1, E2, G12, v12, thi, ang5)) #4
						l = lam.Laminate(plyList, half=True, odd=False)
				
						valid = True
						for load, maxStrain in zip(p,e):
							f = lam.Loadcase(l,load,maxStrain)
							if not f.validLoad():
								valid = False
						if valid:
							l.fe = f.fe
							validLayups.append(l)
						nbTested += 1
	
	# Print result
	if len(validLayups):
		print('\n{} laminates with {} layers were tested, {} meets the requirements:\n'
			.format(nbTested,len(validLayups[0].plys),len(validLayups)))
		for layup in validLayups:
			print("{:35s} fe: {:.2}".format(str(layup),layup.fe))
	else:
		print('\n{} laminates were tested, none found valid\n'.format(nbTested))
		
def problem2():
	header(2)
	
	plyList = lam.plyListFromCSV('layup.csv')
	for ply in plyList:
		ply.setStrength(990, 760, 35, 160, 55)
		ply.setDegradingFactors(1, 0.1, 0.1, 0.5)
	
	# Make three identical copies of laminate
	l = lam.Laminate(plyList)
	g = copy.deepcopy(l)
	h = copy.deepcopy(l)
	
	# Plot reactionforce in x/y/xy directions and save figure
	lam.plotReactionForce(l, 0, '..\img\plot1')
	lam.plotReactionForce(g, 1, '..\img\plot2')
	lam.plotReactionForce(h, 2, '..\img\plot3')

def header(nb):
	print '\n','#'*10,'\n','PROBLEM ',nb,'\n','#'*10
	
if __name__ == '__main__':
    problem1()
    problem2()