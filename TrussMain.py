#!/bin/env/python
# David Wright
# Copyright 2017
# Written for Python 3.5.2



#Import 
import sys
from PyQt5.QtWidgets import (QWidget, QTreeView, QMessageBox, QHBoxLayout, 
							 QFileDialog, QLabel, QSlider, QCheckBox, 
							 QLineEdit, QVBoxLayout, QApplication, QPushButton,
							 QTableWidget, QTableWidgetItem,QSizePolicy,
							 QGridLayout,QGroupBox, QMainWindow,QAction)
from PyQt5.QtCore import Qt, QTimer, QCoreApplication
import TrussUI
import TrussOptimize
from scipy import stats
from sympy import *
import statistics
import numpy as np
import os
import cmath
import math
import Optimization as opt
import FunctionApproximation as approx
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import rcParams
import matplotlib.mlab as mlab
from PyQt5.QtCore import QThread
from PyQt5 import QtCore, QtGui

class MainWindow(QMainWindow, TrussUI.Ui_MainWindow):
	def __init__(self):
		super(self.__class__, self).__init__()
		self.setupUi(self)
		self.runButton.clicked.connect(self.startOptimization)
		# Initialize variables
		self.design = [1,2,3,4]
		self.damping = 1.0
		self.iterations = 0
		self.maxIterations = 100
		self.mechanismStartAngle = 0.017

		self.dampingSlider.valueChanged.connect(self.dampingChanged)
		self.data_table.itemSelectionChanged.connect(self.selectedTableItem)
		self.data_table.cellChanged.connect(self.cellChanged)

	def startOptimization(self):
		# We have a list of subreddits which we use to create a new getPostsThread
		# instance and we pass that list to the thread
		self.optimizeThread = TrussOptimize.OptimizeThread(self.design)

		# Next we need to connect the events from that thread to functions we want
		# to be run when those signals get fired
		self.optimizeThread.iterationDone.connect(self.iterationDone)
		self.optimizeThread.designOptimized.connect(self.designOptimized)

		# We have all the events we need connected we can start the thread
		self.optimizeThread.start()
		# At this point we want to allow user to stop/terminate the thread
		# so we enable that button
		#self.btn_stop.setEnabled(True)
		# And we connect the click of that button to the built in
		# terminate method that all QThread instances have
		#self.btn_stop.clicked.connect(self.optimizeThread.terminate)
		# We don't want to enable user to start another thread while this one is
		# running so we disable the start button.
		#self.btn_start.setEnabled(False)

	def load_data(self):
		#Write this function
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		fileName, _ = QFileDialog.getOpenFileName(self,
			"QFileDialog.getOpenFileName()", 
			"","All Files (*);;Python Files (*.py)", 
			options=options)
		if fileName:
			f= open(fileName,"r")
			if f.mode == 'r':
				contents =f.read()
				# Do stuff with contents
	
	def save_data(self):
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		fileName, _ = QFileDialog.getSaveFileName(self,
			"QFileDialog.getSaveFileName()","",
			"All Files (*);;Text Files (*.txt)", 
			options=options)
		# Check to make sure ends in .txt
		if fileName:
			f= open(fileName,"w+")
			f.close() 

	def dampingChanged(self,value):
		self.damping = float(value)/1000
		self.dampingLabel.setText("Damping = %1.2f"%(self.damping))

	def cellChanged(self,row,column):
		"""
		if row == len(self.xTargets):
			if column == 1:
				cell = self.data_table.item(row, column)
				cellText = cell.text()
				cellValue = float(cellText)
				self.thetaTargets.append(cellValue)
				self.xTargets.append(0)
				self.yTargets.append(0)
				self.exact.append(False)
			elif column == 2:
				cell = self.data_table.item(row, column)
				cellText = cell.text()
				cellValue = float(cellText)
				self.xTargets.append(cellValue) 
				self.thetaTargets.append(0)
				self.ytargets.append(0)
				self.exact.append(False)
			elif column == 3:
				cell = self.data_table.item(row, column)
				cellText = cell.text()
				cellValue = float(cellText)
				self.thetaTargets.append(0)
				self.xTargets.append(0)
				self.yTargets.append(cellValue)            
				self.exact.append(False)
		else:
			if column == 1:
				cell = self.data_table.item(row, column)
				cellText = cell.text()
				cellValue = float(cellText)
				self.thetaTargets[row] = cellValue
			elif column == 2:
				cell = self.data_table.item(row, column)
				cellText = cell.text()
				cellValue = float(cellText)
				self.xTargets[row] = cellValue
			elif column == 3:
				cell = self.data_table.item(row, column)
				cellText = cell.text()
				cellValue = float(cellText)
				self.yTargets[row] = cellValue
		self.data_table.setRowCount(len(self.thetaTargets)+1)
		targets = [self.thetaTargets,self.xTargets,self.yTargets]
		path = [self.pathX,self.pathY]
		self.graph_canvas.plotFourBar(targets, self.mechanismPoints, path)
		"""
		print('Cell Changed')

	def selectedTableItem(self):
		#print("\n")
		for currentQTableWidgetItem in self.data_table.selectedItems():
			#print(currentQTableWidgetItem.row(), currentQTableWidgetItem.column(), currentQTableWidgetItem.text())
			if currentQTableWidgetItem.column() == 0:
				cellValue = self.exact[currentQTableWidgetItem.row()]
				if cellValue == True:
					self.exact[currentQTableWidgetItem.row()] = not self.exact[currentQTableWidgetItem.row()]
				elif cellValue == False:
					if self.exactSelected < 3:
						self.exact[currentQTableWidgetItem.row()] = not self.exact[currentQTableWidgetItem.row()]
			else:
				self.data_table.editItem(currentQTableWidgetItem)

			if self.exact[currentQTableWidgetItem.row()] == True:
				self.data_table.setItem(currentQTableWidgetItem.row(),0, QTableWidgetItem('\u2714'))
			else:
				self.data_table.setItem(currentQTableWidgetItem.row(),0, QTableWidgetItem(''))
		self.exactSelected = sum(self.exact)        
		self.data_table.clearSelection()
		#self.redrawTable()

	def redrawTable(self):
		print('Redrawing Table')

	def redrawResultsLabels(self):
		self.iterationLabel.setText("Iteration %i"%(self.damping))
		#self.length1Label.setText("Length 1 = %2.4f"%(self.mechanismLengths[1]))
		#self.length2Label.setText("Length 2 = %2.4f"%(self.mechanismLengths[2]))
		#self.length3Label.setText("Length 3 = %2.4f"%(self.mechanismLengths[3]))
		#self.length4Label.setText("Length 4 = %2.4f"%(self.mechanismLengths[4]))
		#self.length5Label.setText("Length 5 = %2.4f"%(self.mechanismLengths[5]))
		#self.base1Label = QLabel("Base 1 Location  = (%2.4f, %2.4f)"%(self.mechanismBases[0][0],self.mechanismBases[0][1]),self)
		#self.base2Label = QLabel("Base 2 Location  = (%2.4f, %2.4f)"%(self.mechanismBases[1][0],self.mechanismBases[1][1]),self)

	def parseDesign(self,design):
		print('going to parse design')

	def packageDesign(self):
		print('going to create design list')

	def iterationDone(self,design):
		print(design)

	def designOptimized(self,design):
		print(design)
class FourBarOptimizer(QWidget):

	def __init__(self,optimizeThread):
		super().__init__()
		self.optimizeThread = optimizeThread
		# Upon startup, run a user interface routine
		self.init_ui()
			  
	def init_ui(self):
		

		#self.thetaTargets = [0,90,180,270]
		#self.xTargets = [1,4,6,3]
		#self.yTargets = [1,0,1,2]
		self.thetaTargets = [0,90]
		self.xTargets = [1,4]
		self.yTargets = [1,0]
		self.targets = [self.thetaTargets,self.xTargets,self.yTargets]
		self.exact = [True,True,False,True]
		self.exactSelected = sum(self.exact)
		self.fixedPointCount = min(len(self.thetaTargets),len(self.xTargets),len(self.yTargets))
		# Mechanism encoding: [z1, z2, z3, z4, z5, x1,y1,x2,y2]
		# or more compact: [[lengths],[base1],[base2],]
		# or point based: [[Base1X,Base1Y],[Base2X,base2Y],[point3],[point4],[point5]]
		self.initialMechanismBases = [[0.06,0.04],[1.5,-0.02]]
		baselength = self.vLen(self.initialMechanismBases[0],self.initialMechanismBases[1])
		#self.initialMechanismLengths = [baselength,0.3,0.75,0.6,1,1.5]
		self.initialMechanismLengths = [baselength,0.25, 1.1,0.83,1.25,1.46]
		self.initialAlpha = 0.75
		
		# Four Bar Properties
		self.mechanismBases = self.initialMechanismBases
		self.mechanismLengths = self.initialMechanismLengths
		self.mechanismPoints = self.calculateFourBarPoint(self.initialAlpha)

		self.plotAngles = np.linspace(0.0,6.28318530718,num=20).tolist()
		self.plotAngles.append(self.plotAngles[0])
		self.pathX = [0]*len(self.plotAngles)
		self.pathY = [0]*len(self.plotAngles)

		#Builds GUI
		# Four input tables: 
		# one with initial node coordinates (3xn: Node,X,Y)
		# one with node connectivity (3xbeams: Element, From Node, To Node)
		# one with reactions locations (3xreactions: Reaction, Node, Direction)
		# and one with external loading (3xForces: On Node, Force, Angle)
		
		# Dynamic plot updates with triangles for reactions, lines for beams and filled circles for nodes, and arrow for applied forces
		# Checks: all nodes have at least one member connectivity
		# 
		# Objective function: Sum(Area[i]*length[i])
		# subject to:
		# max(stresses) < maxStress
		# Any locaton constraints, such as: (generated checklist?)
		# Node[1][0] = 1

		self.setGeometry(200,30,1000,700)
	

	def calculateActualPath(self,angleList):
		xActual = [0]*len(angleList)
		yActual = [0]*len(angleList)
		for i in range(0,len(angleList)):
			points = self.calculateFourBarPoint(angleList[i])
			xActual[i] = points[4][0]
			yActual[i] = points[4][1]
		self.pathX = xActual
		self.pathY = yActual

	def fourBarExpression(self):
		# variables = b1x,b1y,b2x,b2y,l0,l1,l2,l3,l4,l5,alpha
		b1x,b1y,b2x,b2y,l0,l1,l2,l3,l4,l5,alpha, a0 = symbols('b1x,b1y,b2x,b2y,l0,l1,l2,l3,l4,l5,alpha,a0')
		p1x = b1x+l1*cos(alpha + a0)
		p1y = b1y+l1*sin(alpha + a0)
		baseAngle = atan((b2y-b1y)/(b2x-b1x))

		x3 = l0*cos(baseAngle) - l1*cos(alpha + a0)
		y3 = l0*sin(baseAngle) - l1*sin(alpha + a0)
		theta3ArccosValue = (x3**2 + y3**2 + (-l3)**2 - l5**2)/2 * (-l3) * sqrt(x3**2 + y3**2)
		theta3 = atan2(y3,x3) + acos(theta3ArccosValue)
		theta5 = atan2((y3-(-l3)*sin(theta3))/l5, (x3-(-l3)*cos(theta3))/l5)

		dyadAngle1 = acos((l5**2 + l2**2 - l4**2)/(2*l5*l2))
		Cx = p1x+l2*cos(theta5 + dyadAngle1)
		Cy = p1y+l2*sin(theta5 + dyadAngle1)

		return (Cx,Cy)

	def calculateFourBarPoint(self,alphaInput):        
		# lengths is in format [base, z1,z2,z3,z4,z5]
		# bases are in format [[base1X,base1Y],[base2X,base2Y]]
		# alpha is a single float for crank angle from horizontal

		bases = list(self.mechanismBases)
		length = list(self.mechanismLengths)
		alpha = alphaInput + self.mechanismStartAngle

		point1 = [0,0]
		point1[0] = bases[0][0] + length[1]*np.cos(alpha)
		point1[1] = bases[0][1] + length[1]*np.sin(alpha)
		
		baseRun = np.float64(bases[1][0] - bases[0][0])
		baseRise = np.float64(bases[1][1]-bases[0][1])
		print(baseRun)
		print(baseRise)
		baseAngle = 0
		if baseRise == 0:
			if baseRun > 0:
				baseAngle = 0
			elif baseRun < 0:
				baseAngle = pi()
		elif baseRun == 0:
			if baseRise > 0:
				baseAngle = pi()/2
			elif baseRise < 0:
				baseAngle = 3*pi()/2
		else:
			baseAngle = np.arctan(baseRise/baseRun)

		x3 = length[0]*cos(baseAngle) - length[1]*cos(alpha)
		y3 = length[0]*sin(baseAngle) - length[1]*sin(alpha)
		theta3ArccosValue = (x3**2 + y3**2 + (-length[3])**2 - length[5]**2)/2 * (-length[3]) * math.sqrt(x3*x3 + y3*y3)
		theta3ArccosValue = np.float64(theta3ArccosValue)
		theta3pos = math.atan2(y3,x3) + np.arccos(theta3ArccosValue)
		#print(theta3pos)
		theta3neg = math.atan2(y3,x3) - np.arccos(theta3ArccosValue)

		theta3 = theta3pos
		theta5 = math.atan2((y3-(-length[3])*np.sin(theta3))/length[5], (x3-(-length[3])*np.cos(theta3))/length[5])
		point2 = [0,0]
		point2[0] = bases[1][0] + length[3]*np.cos(theta3)
		point2[1] = bases[1][1] + length[3]*np.sin(theta3)

		dyadAngle1 = np.arccos(np.float64(length[5]**2 + length[2]**2 - length[4]**2)/np.float64(2*length[5]*length[2]))
		pointC = [0,0]
		pointC[0] = point1[0]+length[2]*np.cos(theta5 + dyadAngle1)
		pointC[1] = point1[1]+length[2]*np.sin(theta5 + dyadAngle1)
		
		self.mechanismPoints = [bases[0],bases[1],point1,point2,pointC]

		return self.mechanismPoints


	def runFourBarOptimization(self):
		print('Starting Optimization')
		(xExpression,yExpression) = self.fourBarExpression()
		objectiveFunction = 0

		for i in range(0,len(self.xTargets)):
			objectiveExpression = sqrt((self.xTargets[i] - xExpression)**2 + (self.yTargets[i] - yExpression)**2)
			objectiveExpression = objectiveExpression.subs(symbols('alpha'),self.thetaTargets[i])
			objectiveFunction = objectiveFunction + objectiveExpression

		variables = ['l0','l1','l2','l3','l4','l5','a0','b1x','b1y','b2x','b2y']
		rpSymbol = symbols('rp')
		g1 = 'l5 - (l2 + l4) + 0.01'
		g2 = '(l1 + l0 + 0.01) - (l5 + l3)'
		g3 = '-l1'
		g4 = '-l2'
		g5 = '-l3'
		g6 = '-l4'
		g7 = '-l5'
		inequalityConstraints = [g1,g2,g3,g4,g5,g6,g7]
		equalityConstraints = []
		startingPoint = [self.mechanismLengths[0],
			self.mechanismLengths[1],
			self.mechanismLengths[2],
			self.mechanismLengths[3],
			self.mechanismLengths[4],
			self.mechanismLengths[5],
			self.mechanismStartAngle,
			self.mechanismBases[0][0],
			self.mechanismBases[0][1],
			self.mechanismBases[1][0],
			self.mechanismBases[1][1]]
		shouldContinue = True
		rp = 1
		testSteps = [0,0.01,0.02]
		expression = objectiveFunction

		print('Set up initial constants and expressions')
		self.iterations = 0
		position = startingPoint
		objectiveValue = opt.evaluateLinearExtendedPenalty(expression, 
			inequalityConstraints=inequalityConstraints, 
			equalityConstraints=equalityConstraints, 
			variables = variables, 
			values = position, 
			rp = rp,  
			evaluate=True)

		print('Calculated initial objective function:')
		print(objectiveValue)
		while shouldContinue == True:
			# perform 1 iteration
			# check validity
			# plot
			########################
			self.iterations = self.iterations + 1
			
			print('Calculating expression for this rp')
			expressionHere = opt.evaluateLinearExtendedPenalty(expression, 
				inequalityConstraints=inequalityConstraints, 
				equalityConstraints=equalityConstraints, 
				variables = variables, 
				values = position, 
				rp = rp, 
				evaluate=False)
			expressionHere = expressionHere.subs(rpSymbol,float(rp))
			print('Calculating gradients:')
			slopeList = opt.getNumGradient(expressionHere,variables,position,normalize=False)
			print(slopeList)
			#print(slopeList)
			# Get three points in that direction at intervals of 0.5,1,2
			#
			#Constant step toward goal
			"""
			functionValues = [objectiveValue]
			for alphaValue in testSteps:
				if alphaValue != testSteps[0]:
					testLocation = []
					for oldPosition, slope in zip(position,slopeList):
						testLocation.append(oldPosition-slope*alphaValue)
					#print('Location to test:')
					#print(testLocation)
					functionValues.append(opt.evaluateLinearExtendedPenalty(expression, 
						inequalityConstraints=inequalityConstraints, 
						equalityConstraints=equalityConstraints, 
						variables = variables, 
						values = testLocation, 
						rp = rp))
			print('Calculated test positions')
			# Fit parabola to curve
			print(functionValues)
			C = approx.threePointQuadraticApprox(testSteps, functionValues)
			print(C)
			# Check parabola is concave up
			# Calculate alpha that gives minimum
			alphaStar = 0.0
			if C[2] < 0.00001:
				print("Fitted parabola is concave down. Minimum alpha value is not bounded.")
				alphaStar = 0.5
			else:
				(alphaStar,bestY) = opt.minimizeParabola(C)
			"""
			alphaStar = 0.1
			print('Calculating New Position')
			# Move to position of calculated alpha
			newPosition = []
			for oldPosition, slope in zip(position,slopeList):
				newPosition.append(oldPosition-slope*self.damping*alphaStar)
			lastPosition = position
			position = newPosition
			objectiveValueLast = objectiveValue
			print('Finding New error')
			objectiveValue = opt.evaluateLinearExtendedPenalty(expression, 
				inequalityConstraints=inequalityConstraints, 
				equalityConstraints=equalityConstraints, 
				variables = variables, 
				values = position, 
				rp = rp, 
				evaluate=True)
			print('New error is:')
			print(objectiveValue)
			print('New Design is:')
			print(position)
			print('Updating Mechanism parameters')
			self.mechanismLengths = [position[0],position[1],position[2],position[3],position[4],position[5]]
			self.mechanismAngle0 = position[6]
			self.mechanismBases = [[position[7],position[8]],[position[9],position[10]]]

			print('Redrawing')
			# Redraw results labels
			self.redrawResultsLabels()
			# Redraw graph
			self.calculateActualPath(self.plotAngles)
			path = [self.pathX,self.pathY]
			mechanismPoints = self.calculateFourBarPoint(self.initialAlpha)
			targets = [self.thetaTargets,self.xTargets,self.yTargets]
			path = [self.pathX,self.pathY]        

			self.graph_canvas.plotFourBar(targets, mechanismPoints,path)

			print('Checking Convergence')
			# Check convergence
			deltaObjective = abs(float(objectiveValueLast - objectiveValue))
			deltaVariables = [abs(new - old) for new, old in zip(position,lastPosition)]
			maxDelta = max(deltaVariables)
			if max(deltaObjective,maxDelta) < epsilon:
				shouldContinue = False
				print("Local Optimium found")

			if self.iterations > self.maxIterations:
				print("Function timed out. Returning final result")
				shouldContinue = False

			#######################

		path = [self.pathX,self.pathY]
		mechanismPoints = self.calculateFourBarPoint(self.initialAlpha)
		targets = [self.thetaTargets,self.xTargets,self.yTargets]
		self.graph_canvas.plotFourBar(targets, self.mechanismPoints,path)

	def vLen(self,point1,point2):
		# takes in two points in the format [x,y] and returns the float of vector length
		dx = point1[0] - point2[0]
		dy = point1[1] - point2[1]

		length = np.sqrt(dx*dx + dy*dy)
		return length


def main():
	app = QApplication(sys.argv)
	form = MainWindow()
	form.show()
	app.exec_()

if __name__ == '__main__':
	main()
