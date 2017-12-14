#!/bin/env/python
# David Wright
# Copyright 2017
# Written for Python 3.5.2

# TODO: Reset view to input parameters
# TODO: Add cross sections to design variables
# TODO: Calculate euler buckling bounds
# TODO: Export to OpenSCAD?
# TODO: Add beam weights as forces to truss nodes (Can it hold itself up?)
# TODO: Add progress bar corellated to Rp value
# TODO: Show total weight
# TODO: Convergence plot?

#Import 
import sys
from PyQt5.QtWidgets import (QWidget, QTreeView, QMessageBox, QHBoxLayout, 
							 QFileDialog, QLabel, QSlider, QCheckBox, 
							 QLineEdit, QVBoxLayout, QApplication, QPushButton,
							 QTableWidget, QTableWidgetItem,QSizePolicy,
							 QGridLayout,QGroupBox, QMainWindow,QAction)
from PyQt5.QtCore import Qt, QTimer, QCoreApplication, QThread
from PyQt5 import QtCore, QtGui
import TrussUI
import TrussOptimize
from scipy import stats
from sympy import *
import numpy as np
import os
import cmath
import math
import FunctionApproximation as approx
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import rcParams
import matplotlib.mlab as mlab
import csv


class MainWindow(QMainWindow, TrussUI.Ui_MainWindow):
	def __init__(self):
		super(self.__class__, self).__init__()
		self.setupUi(self)
		self.startButton.clicked.connect(self.startOptimization)
		# Initialize variables
		# [[X,Y,Fix X, Fix Y, Rx, Ry,Applied Force, Force Angle],[X,Y,Fix X, Fix Y, Rx, Ry,Applied Force, Force Angle]]
		self.initialNodeArray = [[0,0,1,1,0,0,1,270],[5,0,1,0,1,1,0,0],[5,3,1,0,1,0,0,0]] 
		self.initialBeamArray = [[0,1,1,0,0],[1,2,1,0,0],[2,0,1,0,0]] # [[From, To, Dim1, Dim2,stress], [From, To, Dim1, Dim2,stress]]

		# 1 = "Rectangular - Equal Thickness":
		# 2 = "Rectangular":
		# 3 = "Rectangular - Hollow":
		# 4 = "Square":
		# 5 = "Square - Hollow":
		# 6 = "Round":
		# 7 = "Round - Hollow":
		self.crossSection = 1 

		self.currentNodeArray = self.initialNodeArray
		self.currentBeamArray = self.initialBeamArray

		self.formerForceArray = []
		for i in range(0,len(self.initialNodeArray)):
			if self.initialNodeArray[i][6] != 0:
				self.formerForceArray.append([i,self.initialNodeArray[i][6],self.initialNodeArray[i][7]])

		self.damping = 0.1
		self.iterations = 0
		self.maxIterations = 100
		self.maxStress = 10.0
		self.density = 1.0
		#self.mechanismStartAngle = 0.017
		self.programLoaded = False
		self.dampingSlider.valueChanged.connect(self.dampingChanged)
		self.nodesTable.itemSelectionChanged.connect(self.selectedNodesTable)
		self.nodesTable.cellChanged.connect(self.nodeCellChanged)
		self.beamTable.itemSelectionChanged.connect(self.selectedBeamTable)
		self.beamTable.cellChanged.connect(self.beamCellChanged)
		self.forceTable.itemSelectionChanged.connect(self.selectedForcesTable)
		self.forceTable.cellChanged.connect(self.forceCellChanged)
		self.maxStressTextBox.textChanged.connect(self.maxStressChanged)	
		self.maxIterationsTextBox.textChanged.connect(self.maxIterationsChanged)	
		self.redrawInputTables()
		self.graph_canvas.plotTruss(self.initialNodeArray,self.initialBeamArray)
		self.programLoaded = True
		self.userEdited = True

	def startOptimization(self):
		# Create calculation thread
		self.design = self.packageDesign()
		controls = [self.damping,self.maxIterations,self.maxStress,self.crossSection,self.density]
		self.optimizeThread = TrussOptimize.OptimizeThread(self.design,controls)
		# Connect to emitted signals
		self.optimizeThread.iterationDone.connect(self.iterationDone)
		self.optimizeThread.designOptimized.connect(self.designOptimized)
		# Start the thread
		self.optimizeThread.start()

		self.stopButton.setEnabled(True)
		self.stopButton.clicked.connect(self.endOptimization)
		self.startButton.setEnabled(False)

	def load_data(self):
		#Write this function
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		fileName, _ = QFileDialog.getOpenFileName(self,
			"QFileDialog.getOpenFileName()", 
			"","All Files (*);;Text Files (*.txt)", 
			options=options)
		if fileName:
			self.initialNodeArray = []
			self.initialBeamArray = []
			with open(fileName, 'r') as fin:
				reader = csv.reader(fin, delimiter=',')
				for line in reader:
					#print(line)
					if len(line) == 8:
						self.initialNodeArray.append([float(line[0]),float(line[1]),int(line[2]),int(line[3]),int(line[4]),int(line[5]),float(line[6]),float(line[7])])
					elif len(line) == 5:
						self.initialBeamArray.append([int(line[0]),int(line[1]),float(line[2]),float(line[3]),float(line[4])])
					else:
						print('Unexpected line length')
			self.redrawInputTables()
			self.graph_canvas.plotTruss(self.initialNodeArray,self.initialBeamArray)

	def saveInputData(self):
		design = [self.initialNodeArray,self.initialBeamArray]
		self.saveDesign(design)

	def saveOptimizedData(self):
		design = [self.currentNodeArray,self.currentBeamArray]
		self.saveDesign(design)

	def saveDesign(self,design):
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		fileName, _ = QFileDialog.getSaveFileName(self,
			"QFileDialog.getSaveFileName()","",
			"All Files (*);;Text Files (*.txt)", 
			options=options)
		# Check to make sure ends in .txt
		suffix = fileName[-4:-1] + fileName[-1]
		if suffix != '.txt':
			fileName = fileName + '.txt'

		if fileName:
			f= open(fileName,"w+")
			for node in design[0]:
				f.write("%f, %f, %i, %i, %i, %i, %f, %f \n" % (node[0],node[1],node[2],node[3],node[4],node[5],node[6],node[7]))
			for beam in design[1]:
				f.write("%i, %i, %f, %f, %f \n" % (beam[0],beam[1],beam[2],beam[3],beam[4]))			
			f.close() 


	def endOptimization(self):
		self.optimizeThread.terminate
		self.stopButton.setEnabled(False)
		self.startButton.setEnabled(True)

	def dampingChanged(self,value):
		self.damping = float(value)/1000
		self.dampingLabel.setText("Damping = %1.2f"%(self.damping))
	
	def maxStressChanged(self,newText):
		if newText != "":
			self.maxStress = float(newText)
		else:
			self.maxStress = 1.0
		#print(self.maxStress)

	def maxIterationsChanged(self,newText):
		if newText != "":
			self.maxIterations = int(newText)
		else:
			self.maxIterations = 1

	def densityChanged(self,newText):
		if newText != "":
			self.density = float(newText)
		else:
			self.density = 1

	def crossSectionChanged(self, text):
		if text == "Rectangular - Equal Thickness":
			self.crossSection = 1
		elif text == "Rectangular":
			self.crossSection = 2
		elif text == "Rectangular - Hollow":
			self.crossSection = 3
		elif text == "Square":
			self.crossSection = 4
		elif text == "Square - Hollow":
			self.crossSection = 5
		elif text == "Round":
			self.crossSection = 6
		elif text == "Round - Hollow":
			self.crossSection = 7

	def nodeCellChanged(self,row,column):
		if self.programLoaded == True:
			if self.userEdited == True:
				self.userEdited = False
			
				# Add a new entry to the nodes list if last line is selected
				if row == len(self.initialNodeArray): 
					cell = self.nodesTable.item(row, column)
					cellText = cell.text()
					if cellText != '':
						cellValue = float(cellText)
					else:
						cellValue = 0

					if column == 1:
						self.initialNodeArray.append([cellValue,0,0,0,0,0,0,0])
					elif column == 2:
						self.initialNodeArray.append([0,cellValue,0,0,0,0,0,0])
				else:
					# Grab float value of text input for columns 0 and 1
					if (column == 0 or column == 1):
						cell = self.nodesTable.item(row, column)
						cellText = cell.text()			
						if cellText != '':
							cellValue = float(cellText)
						else:
							cellValue = 0
						self.initialNodeArray[row][column] = cellValue
					# SelectedNodesTable already took care of checkmark assignment
				self.nodesTable.setRowCount(len(self.initialNodeArray)+1)
				self.graph_canvas.plotTruss(self.initialNodeArray,self.initialBeamArray)
				self.userEdited = True

	def selectedNodesTable(self):
		for currentQTableWidgetItem in self.nodesTable.selectedItems():
			row = currentQTableWidgetItem.row()
			col = currentQTableWidgetItem.column()
			# Invert chekcboxes and save for columns 2, 3, 4, and 5
			if (col == 2 or col == 3 or col == 4 or col == 5):
				cellValue = self.initialNodeArray[row][col]
				self.initialNodeArray[row][col] = not self.initialNodeArray[row][col]

				self.userEdited = False
				if self.initialNodeArray[row][col] == True:
					self.nodesTable.setItem(row,col, QTableWidgetItem('\u2714'))
				else:
					self.nodesTable.setItem(row,col, QTableWidgetItem(' '))
				self.userEdited = True
				self.graph_canvas.plotTruss(self.initialNodeArray,self.initialBeamArray)
			# Edit values and save for columns 0 and 1
			elif (col == 0 or col == 1):
				self.nodesTable.editItem(currentQTableWidgetItem)

		self.nodesTable.clearSelection()
		#self.redrawNodeTable()

	def beamCellChanged(self,row,column):
		# Add a new entry to the nodes list if last line is selected
		if self.programLoaded == True:
			if self.userEdited == True:
				self.userEdited = False
				if row == len(self.initialBeamArray): 
					cell = self.beamTable.item(row, column)
					cellText = cell.text()
					if cellText != '':
						cellValue = int(cellText)
					else:
						cellValue = 0
					#cell.setText('%i'%(cellValue))

					if column == 1:
						self.initialBeamArray.append([cellValue-1,0,1,0,0])
						self.beamTable.setItem(row,1, QTableWidgetItem('%i'%(self.initialBeamArray[row][1]+1)))	

					elif column == 2:
						self.initialBeamArray.append([0,cellValue-1,1,0,0])
						self.beamTable.setItem(row,0, QTableWidgetItem('%i'%(self.initialBeamArray[row][0]+1)))
					
					self.beamTable.setRowCount(len(self.initialBeamArray)+1)
				else:
					# No action for column zero
					# Grab integer value of text input for columns 1 and 2
					if (column == 1 or column == 2):
						cell = self.beamTable.item(row, column)
						cellText = cell.text()			
						if cellText != '':
							cellValue = int(cellText)
						else:
							cellValue = 0
						self.initialBeamArray[row][column-1] = cellValue - 1
				# SelectedNodesTable already took care of checkmark assignment
				self.graph_canvas.plotTruss(self.initialNodeArray,self.initialBeamArray)
				self.userEdited = True

	def selectedBeamTable(self):
		for currentQTableWidgetItem in self.beamTable.selectedItems():
			row = currentQTableWidgetItem.row()
			col = currentQTableWidgetItem.column()
			# Do nothing for column 0
			# Edit values and save for columns 1 and 2
			if (col == 1 or col == 2):
				self.beamTable.editItem(currentQTableWidgetItem)
		self.beamTable.clearSelection()
		#self.redrawTable()

	def forceCellChanged(self,row,column):
		# Add a new entry to the nodes list if last line is selected
		# [X,Y,Fix X, Fix Y, Rx, Ry,Applied Force, Force Angle]
		if self.programLoaded == True:
			if self.userEdited == True:
				self.userEdited = False
				forcesArray = []
				nodeCount = len(self.initialNodeArray)

				for i in range(0,nodeCount):
					thisForce = self.initialNodeArray[i][6]
					thisAngle = self.initialNodeArray[i][7]
					if thisForce != 0:
						forcesArray.append([i,thisForce,thisAngle])

				self.formerForceArray = forcesArray
				if row == len(forcesArray): 
					cell = self.forceTable.item(row, column)
					cellText = cell.text()
					if cellText != '':
						cellValue = float(cellText)
					else:
						cellValue = 0
					#cell.setText('%i'%(cellValue))
					if column == 0:				
						self.forceTable.setItem(row,1, QTableWidgetItem('1'))	
						self.forceTable.setItem(row,2, QTableWidgetItem('0'))	
						self.initialNodeArray[int(cellValue)-1][6] = 1
						self.initialNodeArray[int(cellValue)-1][7] = 0
					if column == 1:
						#self.forceTable.setItem(row,0, QTableWidgetItem('1'))	
						#self.forceTable.setItem(row,2, QTableWidgetItem('0'))	
						self.initialNodeArray[0][6] = cellValue
						self.initialNodeArray[0][7] = 0

					elif column == 2:
						#self.forceTable.setItem(row,0, QTableWidgetItem('1'))	
						#self.forceTable.setItem(row,1, QTableWidgetItem('0'))	
						self.initialNodeArray[0][6] = 0
						self.initialNodeArray[0][7] = cellValue
					
					self.forceTable.setRowCount(len(forcesArray)+1)
				else:
					# Grab integer value of text input for columns 1 and 2
					cell = self.forceTable.item(row, column)
					cellText = cell.text()			
					if cellText != '':
						cellValue = float(cellText)
					else:
						cellValue = 0
					
					if column == 0:
						cellValue = int(cellValue) - 1
						# Copy values to new node
						newNode = int(cellValue)
						oldNode = int(self.formerForceArray[row][0])
						self.initialNodeArray[newNode][6] = self.formerForceArray[row][1]
						self.initialNodeArray[newNode][7] = self.formerForceArray[row][2]
						# Clear values from old node
						self.initialNodeArray[oldNode][6] = 0
						self.initialNodeArray[oldNode][7] = 0
					else:
						forceNode = int(self.forceTable.item(row, 0).text())-1
						self.initialNodeArray[forceNode][column+5] = cellValue

				self.graph_canvas.plotTruss(self.initialNodeArray,self.initialBeamArray)
				self.redrawForceTable()
				self.userEdited = True

	def selectedForcesTable(self):
		for currentQTableWidgetItem in self.forceTable.selectedItems():
			row = currentQTableWidgetItem.row()
			col = currentQTableWidgetItem.column()
			# Do nothing for column 0
			# Edit values and save for columns 1 and 2
			self.forceTable.editItem(currentQTableWidgetItem)
		self.forceTable.clearSelection()
		#self.redrawTable()

	def redrawInputTables(self):
		self.userEdited = False
		self.redrawNodeTable()
		self.redrawBeamTable()
		self.redrawForceTable()
		self.userEdited = True

	def redrawNodeTable(self):
		# Node Table
		self.nodesTable.setRowCount(len(self.initialNodeArray)+1)
		for i in range(0,len(self.initialNodeArray)):
			self.nodesTable.setItem(i,0, QTableWidgetItem('%2.3f'%(self.initialNodeArray[i][0])))
			self.nodesTable.setItem(i,1, QTableWidgetItem('%2.3f'%(self.initialNodeArray[i][1])))	
			# Fix X Checkboxes
			if self.initialNodeArray[i][2] == True:
				self.nodesTable.setItem(i,2, QTableWidgetItem('\u2714'))
			else:
				self.nodesTable.setItem(i,2, QTableWidgetItem(''))
			# Fix Y Checkboxes
			if self.initialNodeArray[i][3] == True:
				self.nodesTable.setItem(i,3, QTableWidgetItem('\u2714'))
			else:
				self.nodesTable.setItem(i,3, QTableWidgetItem(''))
			# Reaction X Checkboxes
			if self.initialNodeArray[i][4] == True:
				self.nodesTable.setItem(i,4, QTableWidgetItem('\u2714'))
			else:
				self.nodesTable.setItem(i,4, QTableWidgetItem(''))
			# Reaction Y Checkboxes
			if self.initialNodeArray[i][5] == True:
				self.nodesTable.setItem(i,5, QTableWidgetItem('\u2714'))
			else:
				self.nodesTable.setItem(i,5, QTableWidgetItem(''))
	
	def redrawBeamTable(self):
		# Beam Table
		self.beamTable.setRowCount(len(self.initialBeamArray)+1)
		for i in range(0,len(self.initialBeamArray)):
			self.beamTable.setItem(i,0, QTableWidgetItem('%i'%(self.initialBeamArray[i][0] + 1)))
			self.beamTable.setItem(i,1, QTableWidgetItem('%i'%(self.initialBeamArray[i][1] + 1)))	

	def redrawForceTable(self):
		# Force Table
		forcesCount = 0
		for i in range(0,len(self.initialNodeArray)):
			if self.initialNodeArray[i][6] != 0:
				forcesCount = forcesCount + 1
		self.forceTable.setRowCount(forcesCount+1)
		currentRow = 0
		for i in range(0,len(self.initialNodeArray)):
			if self.initialNodeArray[i][6] != 0:
				self.forceTable.setItem(i,0, QTableWidgetItem('%i'%(currentRow+1)))
				self.forceTable.setItem(i,1, QTableWidgetItem('%2.3f'%(self.initialNodeArray[i][6])))
				self.forceTable.setItem(i,2, QTableWidgetItem('%2.3f'%(self.initialNodeArray[i][7])))
				currentRow = currentRow + 1

	def redrawResultsTables(self):
		self.redrawNodeResultsTable()
		self.redrawBeamResultsTable()

	def redrawNodeResultsTable(self):
		self.nodesResultsTable.setRowCount(len(self.currentNodeArray))
		for i in range(0,len(self.currentNodeArray)):
			self.nodesResultsTable.setItem(i,0, QTableWidgetItem('%i'%(i+1)))
			self.nodesResultsTable.setItem(i,1, QTableWidgetItem('%2.4f'%(self.currentNodeArray[i][0])))
			self.nodesResultsTable.setItem(i,2, QTableWidgetItem('%2.4f'%(self.currentNodeArray[i][1])))	

	def redrawBeamResultsTable(self):
		self.beamResultsTable.setRowCount(len(self.currentBeamArray))
		for i in range(0,len(self.currentBeamArray)):
			# Calculate length
			fromNode = self.currentBeamArray[i][0]
			toNode = self.currentBeamArray[i][1]
			# X Component
			xFrom = self.currentNodeArray[fromNode][0]
			xTo = self.currentNodeArray[toNode][0]
			# Y Component
			yFrom = self.currentNodeArray[fromNode][1]
			yTo = self.currentNodeArray[toNode][1]
			beamLen = self.vLen([xFrom,yFrom],[xTo,yTo])
			self.beamResultsTable.setItem(i,0, QTableWidgetItem('%2.2f'%(beamLen)))
			# Areas
			self.beamResultsTable.setItem(i,1, QTableWidgetItem('%2.2f'%(self.currentBeamArray[i][2])))	
			self.beamResultsTable.setItem(i,2, QTableWidgetItem('%2.2f'%(self.currentBeamArray[i][3])))	
			# Stresses
			self.beamResultsTable.setItem(i,3, QTableWidgetItem('%2.2f'%(self.currentBeamArray[i][4])))	

	def parseDesign(self,design):
		print('Do we need to parse the design?')

	def packageDesign(self):
		return [self.initialNodeArray,self.initialBeamArray]

	def iterationDone(self,design):
		self.currentNodeArray = list(design[0])
		self.currentBeamArray = list(design[1])
		self.redrawResultsTables()
		self.graph_canvas.plotTruss(self.currentNodeArray,self.currentBeamArray)

	def designOptimized(self,design):
		self.startButton.setEnabled(True)
		self.stopButton.setEnabled(False)
		self.optimizeThread.terminate
		#print(design)

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
