from PyQt5.QtCore import QThread
import PyQt5.QtCore as QtCore
from sympy import *
import sympy as sy
import numpy as np
import math
from numpy.linalg import inv,pinv

class OptimizeThread(QThread):

	iterationDone = QtCore.Signal(object)
	designOptimized = QtCore.Signal(object)
	def __init__(self, design, controls):
		QThread.__init__(self)
		self.design = list(design)
		self.optimizedDesign = self.design

		self.objectiveSequence = []
		#self.variables = list(variables)
		self.nodes = 0
		self.beams = 0
		# Controls = [damping, max Iterations, max Stress, stiffness, cross section type, density]
		self.damping = controls[0]
		self.maxIterations = controls[1]
		self.maxStress = controls[2]
		# 1 = "Rectangular - Equal Thickness":
		# 2 = "Rectangular":
		# 3 = "Rectangular - Hollow":
		# 4 = "Square":
		# 5 = "Square - Hollow":
		# 6 = "Round":
		# 7 = "Round - Hollow":
		self.crossSection = controls[3]
		self.density = controls[4]
	def __del__(self):
		self.wait()

	def runOptimization(self, design):
		# design = [nodeArray, beamArray]
		# Nodes = [X,Y,Fix X, Fix Y, Rx, Ry,Applied Force, Force Angle]
		# Beams = [From, To, Dim1, Dim2,stress]
		self.nodes = design[0]
		self.beams = design[1]
		# Design parameter list creation
		designInitialPosition = []
		for i in range(0,len(self.nodes)):
			if self.nodes[i][2] == 0:
				designInitialPosition.append(self.nodes[i][0])
			if self.nodes[i][3] == 0:
				designInitialPosition.append(self.nodes[i][1])

		# Constraint functions
		# Tensile stress < TMax
		# Compressive stress < CCrit
		# (2*nodes) x (beams + reactions )
		beamForces = self.trussForceAnalysis(self.nodes,self.beams)
		#print(beamForces)
		#print(expression)
		#print(nodePositionVars)
		#print(designInitialPosition)
		rp = 1
		rpMax = 1000000

		optimumParameters = list(designInitialPosition)
		n = 0.
		nMax = 34.
		while rp < rpMax:
			n = n+1.0
			(weight,optimumParameters) = self.numericalMin(optimumParameters,
				rp=rp,
				echo=False,
				damping=self.damping,
				epsilon=0.00001,
				nMax=self.maxIterations,
				alpha = [0,0.01,0.02],
				printResults=False)
			
			# Return to main thread
			optimizedNodes = []
			for node in self.nodes:
				optimizedNodes.append([float(x) for x in node])
			
			j = 0
			for i in range(0,len(self.nodes)):
				if self.nodes[i][2] == 0:
					optimizedNodes[i][0] = float(optimumParameters[j])
					j = j + 1
				if self.nodes[i][3] == 0:
					optimizedNodes[i][1] = float(optimumParameters[j])
					j = j + 1
			beamForces = self.trussForceAnalysis(optimizedNodes,self.beams)
			for i in range(0,len(self.beams)):
				self.beams[i][4] = beamForces[i]
			progress = float(n/nMax)
			self.optimizedDesign = [optimizedNodes,self.beams,progress]
			self.iterationDone.emit(self.optimizedDesign)
			rp = rp*1.5
			self.objectiveSequence.append(weight)

		return self.optimizedDesign

	def trussForceAnalysis(self,nodes,beams):
		reactions = 0
		for i in range(0,len(nodes)):
			if nodes[i][4] != 0:
				reactions = reactions + 1
			if nodes[i][5] != 0:
				reactions = reactions + 1

		M = np.zeros((2*len(nodes),len(beams) + reactions))
		F = np.zeros(2*len(nodes))

		for i in range(0,len(beams)):
			fromNode = beams[i][0]
			toNode = beams[i][1]
			fromX = nodes[fromNode][0]
			toX = nodes[toNode][0]
			fromY = nodes[fromNode][1]
			toY = nodes[toNode][1]
			dx = toX - fromX
			dy = toY - fromY
			length = np.sqrt(dx**2 + dy**2)
			M[2*fromNode,i] = dx/length
			M[2*fromNode+1,i] = dy/length
			M[2*toNode,i] = -dx/length
			M[2*toNode+1,i] = -dy/length

		j = 0
		for i in range(0,len(nodes)):
			# X reactions
			if nodes[i][4] != 0:
				M[2*i,len(beams) + j] = nodes[i][4]
				j = j + 1
			if nodes[i][5] != 0:
				M[2*i+1,len(beams) + j] = nodes[i][5]
				j = j + 1
			if nodes[i][6] != 0:
				F[2*i] = nodes[i][6]*np.cos(math.radians(nodes[i][7])) 
				F[2*i+1] = nodes[i][6]*np.sin(math.radians(nodes[i][7])) 
		MMatrix = np.matrix(M)
		Minv = pinv(MMatrix)

		beamForces = Minv.dot(F)
		beamForces = beamForces.tolist()
		return beamForces[0]

	def evaluateObjectiveFunction(self,design,crossSection=1,fMax = 10,fMin = 10,rp=1):
		
		theseNodes = list(self.nodes)
		j = 0
		for i in range(0,len(self.nodes)):
			if self.nodes[i][2] == 0:
				theseNodes[i][0] = float(design[j])
				j = j + 1
			if self.nodes[i][3] == 0:
				theseNodes[i][1] = float(design[j])
				j = j + 1
		
		theseBeams = list(self.beams)

		forces =  self.trussForceAnalysis(theseNodes,theseBeams)
		# TODO: switch this depending on cross section parameter
		fMax = self.maxStress
		fMin = self.maxStress
		area = 1
		Ival = 1

		epsilon = -0.2*np.sqrt(1/rp)

		objective = 0

		for i in range(0,len(theseBeams)):
			fromNode = theseBeams[i][0]
			toNode = theseBeams[i][1]
			# X Component
			xFrom = theseNodes[fromNode][0]
			xTo = theseNodes[toNode][0]

			# Y Component
			yFrom = theseNodes[fromNode][1]
			yTo = theseNodes[toNode][1]

			objective = objective + area*self.density*self.vLen([xFrom,yFrom],[xTo,yTo])

			if forces[i] > 0:
				# Maximum tensile stress constraint
				maxConstraintValue = forces[i] - fMax
				if maxConstraintValue > epsilon:
					inequalityModifier =  - (2*epsilon - maxConstraintValue)/epsilon**2                  
				else:
					inequalityModifier =  - 1/maxConstraintValue
			else:
				# Maximum compressive stress/buckling constraint
				minConstraintValue = -forces[i] - fMin
				if minConstraintValue > epsilon:
					inequalityModifier =  - (2*epsilon - minConstraintValue)/epsilon**2                  
				else:
					inequalityModifier =  - 1/minConstraintValue
			
			objective = objective + inequalityModifier/rp

		return objective

	def vLen(self,point1,point2):
		# takes in two points in the format [x,y] and returns the float of vector length
		dx = point1[0] - point2[0]
		dy = point1[1] - point2[1]

		length = np.sqrt(dx*dx + dy*dy)
		return length

	def getNumGradient(self,design,rp = 1,delta = 0.0001):
		testDesign = list(design)
		for v in range(0,len(design)):
			newVariableValue = design[v] - delta/2.0
			testDesign[v] = newVariableValue

		startingValue = self.evaluateObjectiveFunction(testDesign, rp = rp)
		slopeList = []
		for v in range(0,len(design)):
			newVariableValue = design[v] + delta/2.0
			testDesign = list(design)
			testDesign[v] = newVariableValue
			testObjectiveValue = self.evaluateObjectiveFunction(testDesign, rp = rp)
			slopeList.append((testObjectiveValue-startingValue)/delta)

		return slopeList

	def threePointQuadraticApprox(self,x,y):
		C2 = (((y[2]-y[0])/(x[2]-x[0])) - ((y[1]-y[0])/(x[1]-x[0])))/(x[2]-x[1])
		C1 = (y[1] - y[0])/(x[1]-x[0]) - C2*(x[0]+x[1])
		C0 = y[0] - C1*x[0] - C2*x[0]**2

		return [C0,C1,C2]

	def minimizeParabola(self,c):
		# Inputs: Coefficients for polynomial equation according to the form C0 + C1*x + C2*x^2...
		# Outputs: Values of x and y where y is minimized
		minX = -c[1]/(2*c[2])

		minY = self.getValueOfPoly(c,minX)
		return (minX,minY)

	def getValueOfPoly(self,c,x):
		#Inputs: Coefficients for polynomial equation according to the form C0 + C1*x + C2*x^2...
		#Inputs: x - value to get value at
		constantQuantity = len(c)

		if constantQuantity == 1:
			# Flat line
			y = c[0]
		elif constantQuantity == 2:
			# Linear
			y = c[0] + c[1] * x
		elif constantQuantity == 3:
			# Quadratic
			y = c[0] + c[1]*x + c[2]*x**2
		elif constantQuantity == 4:
			# Cubic
			y = c[0] + c[1]*x + c[2]*x**2 + c[3]*x**3
		else:
			print("Polynomial could not be calculated. Check getValueOfPoly function.")
			y = 99999999

		return y

	def numericalMin(self,design,rp=1,echo=False,damping=0.1,epsilon=0.001,nMax=100,alpha = [0,0.01,0.02],printResults=False):

		# Constants
		alphaStarMax = 10
		deltaMax = 0.25 + 100/rp
		i = 0
		oscillationEpsilon = epsilon/100

		# Loop
		shouldContinue = True
		position = list(design)
		objectiveValue = self.evaluateObjectiveFunction(position,rp=rp)
		lastPosition = position
		if echo == True:
			headerString = "Iteration\t"
			headerString += "Design\t\t\t"
			#headerString += "Gradient\t"        
			headerString += "F(x)"
			print(headerString)

		while shouldContinue == True:
			i = i+1
		
			#print("Total Iterations should be  %i" %(nMax))
			#print("Iteration %i" %(i))
			# Get gradient at position
			# print("About to get gradient")

			slopeList = self.getNumGradient(design,rp=rp)
			#print("Slope values from getGradient")
			#print(slopeList)
			# print("About to fit polynomial")
			# Get three points in that direction at intervals of 0.5,1,2
			functionValues = [objectiveValue]
			for alphaValue in alpha:
				if alphaValue != alpha[0]:
					testLocation = []
					for oldPosition, slope in zip(position,slopeList):
						testLocation.append(oldPosition-slope*alphaValue)
					functionValues.append(self.evaluateObjectiveFunction(testLocation,rp=rp))


			# Fit parabola to curve
			C = self.threePointQuadraticApprox(alpha, functionValues)
			# Check parabola is concave up
			# Calculate alpha that gives minimum
			alphaStar = 0.0
			if C[2] < 0:
				if echo == True:
					print("Fitted parabola is concave down. Minimum alpha value is not bounded.")
				alphaStar = 0.1
			else:
				(alphaStar,bestY) = self.minimizeParabola(C)

			if alphaStar > alphaStarMax:
				alphaStar = alphaStarMax
			# Move to position of calculated alpha
			newPosition = []
			for oldPosition, slope in zip(position,slopeList):
				moveDist = slope*damping*alphaStar
				if abs(moveDist) < deltaMax:
					newPosition.append(oldPosition-moveDist)
				else:
					newPosition.append(oldPosition-moveDist/abs(moveDist)*deltaMax)

			twoPositionsAgo = lastPosition
			lastPosition = position
			position = newPosition
			objectiveValueLast = objectiveValue
			objectiveValue = self.evaluateObjectiveFunction(position,rp=rp)

			# Print current iteration results
			if echo == True:
				resultsString = "%i        \t" %(i)
				resultsString += "{}\t".format(position)
				resultsString += "%2.6f" % (objectiveValue)
				print(resultsString)

			# Check convergence
			deltaObjective = objectiveValueLast - objectiveValue
			variableChanges = [abs(x - y) for x, y in zip(position,lastPosition)]
			deltaVar = float(max(variableChanges))
			#print("Delta Objective = %2.4f" % (float(deltaObjective)))
			if abs(float(deltaObjective)) < epsilon and i > 1 and deltaVar < epsilon:
				shouldContinue = False
				if printResults == True:
					print("Local Optimium found")

			# Check oscillation
			oscillationDelta = [abs(x - y) for x, y in zip(position,twoPositionsAgo)]
			oscillationDelta = max(oscillationDelta)
				
			if oscillationDelta < oscillationEpsilon:
				damping = damping*0.75

			#print("About to check iteration maximum")
			if i > nMax:
				if printResults == True:
					print("Function timed out. Returning final result")
				shouldContinue = False

		if printResults==True:
			print("#### - Results - ####")
			for variableValue in position:
				print("%2.6f" % (variableValue))
			print("F = %2.6f" % (objectiveValue))

		return (objectiveValue, position)


	def run(self):
		self.design = self.runOptimization(self.design)
		self.designOptimized.emit(self.design)
		
		self.sleep(2)
		