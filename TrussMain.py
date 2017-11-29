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

class MainWindow(QMainWindow, TrussUI.Ui_MainWindow):
	def __init__(self):
		super(self.__class__, self).__init__()
		self.setupUi(self)
		self.startButton.clicked.connect(self.startOptimization)
		# Initialize variables
		# [[X,Y,Fix X, Fix Y, Rx, Ry,Applied Force, Force Angle],[X,Y,Fix X, Fix Y, Rx, Ry,Applied Force, Force Angle]]
		self.initialNodeArray = [[0,0,1,1,0,0,1,270],[5,0,1,0,1,1,0,0],[5,3,1,0,1,0,0,0]] 
		self.initialBeamArray = [[0,1,1,0],[1,2,1,0],[2,0,1,0]] # [[From, To, Dim1, Dim2], [From, To, Dim1, Dim2]]

		self.formerForceArray = []
		for i in range(0,len(self.initialNodeArray)):
			if self.initialNodeArray[i][6] != 0:
				self.formerForceArray.append([i,self.initialNodeArray[i][6],self.initialNodeArray[i][7]])

		self.damping = 1.0
		self.iterations = 0
		self.maxIterations = 100
		self.mechanismStartAngle = 0.017
		self.programLoaded = False
		self.dampingSlider.valueChanged.connect(self.dampingChanged)
		self.nodesTable.itemSelectionChanged.connect(self.selectedNodesTable)
		self.nodesTable.cellChanged.connect(self.nodeCellChanged)
		self.beamTable.itemSelectionChanged.connect(self.selectedBeamTable)
		self.beamTable.cellChanged.connect(self.beamCellChanged)
		self.forceTable.itemSelectionChanged.connect(self.selectedForcesTable)
		self.forceTable.cellChanged.connect(self.forceCellChanged)
		self.redrawInputTables()
		self.graph_canvas.plotTruss(self.initialNodeArray,self.initialBeamArray)
		self.programLoaded = True
		self.userEdited = True

	def startOptimization(self):
		# Create calculation thread
		self.design = self.packageDesign()
		self.optimizeThread = TrussOptimize.OptimizeThread(self.design)
		# Connect to emitted signals
		self.optimizeThread.iterationDone.connect(self.iterationDone)
		self.optimizeThread.designOptimized.connect(self.designOptimized)
		# Start the thread
		self.optimizeThread.start()

		self.stopButton.setEnabled(True)
		self.stopButton.clicked.connect(self.optimizeThread.terminate)
		self.startButton.setEnabled(False)

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

	def nodeCellChanged(self,row,column):
		if self.programLoaded == True:
			
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
				
				self.nodesTable.setItem(row,0, QTableWidgetItem('%i'%(row+1)))
			else:
				# No action for column zero
				# Grab float value of text input for columns 1 and 2
				if (column == 1 or column == 2):
					cell = self.nodesTable.item(row, column)
					cellText = cell.text()			
					if cellText != '':
						cellValue = float(cellText)
					else:
						cellValue = 0
					self.initialNodeArray[row][column-1] = cellValue
				# SelectedNodesTable already took care of checkmark assignment
			self.nodesTable.setRowCount(len(self.initialNodeArray)+1)
			self.graph_canvas.plotTruss(self.initialNodeArray,self.initialBeamArray)

	def selectedNodesTable(self):
		for currentQTableWidgetItem in self.nodesTable.selectedItems():
			row = currentQTableWidgetItem.row()
			col = currentQTableWidgetItem.column()
			# Do nothing for column 0
			# Edit values and save for columns 1 and 2

			# Invert chekcboxes and save for columns 3, 4, 5, and 6
			if (col == 3 or col == 4 or col == 5 or col == 6):
				cellValue = self.initialNodeArray[row][col-1]
				self.initialNodeArray[row][col-1] = not self.initialNodeArray[row][col-1]

				if self.initialNodeArray[row][col-1] == True:
					self.nodesTable.setItem(row,col, QTableWidgetItem('\u2714'))
				else:
					self.nodesTable.setItem(row,col, QTableWidgetItem(''))
			elif (col == 1 or col == 2):
				self.nodesTable.editItem(currentQTableWidgetItem)

		self.nodesTable.clearSelection()
		#self.redrawTable()

	def beamCellChanged(self,row,column):
		# Add a new entry to the nodes list if last line is selected
		if self.programLoaded == True:
			if row == len(self.initialBeamArray): 
				cell = self.beamTable.item(row, column)
				cellText = cell.text()
				if cellText != '':
					cellValue = int(cellText)
				else:
					cellValue = 0
				#cell.setText('%i'%(cellValue))

				if column == 1:
					self.initialBeamArray.append([cellValue-1,0,1,1])
					self.beamTable.setItem(row,2, QTableWidgetItem('%i'%(self.initialBeamArray[row][1]+1)))	

				elif column == 2:
					self.initialBeamArray.append([0,cellValue-1,1,1])
					self.beamTable.setItem(row,1, QTableWidgetItem('%i'%(self.initialBeamArray[row][0]+1)))
				
				self.beamTable.setItem(row,0, QTableWidgetItem('%i'%(row+1)))
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
		self.redrawNodeTable()
		self.redrawBeamTable()
		self.redrawForceTable()

	def redrawNodeTable(self):
		# Node Table
		self.nodesTable.setRowCount(len(self.initialNodeArray)+1)
		for i in range(0,len(self.initialNodeArray)):
			self.nodesTable.setItem(i,0, QTableWidgetItem('%i'%(i+1)))
			self.nodesTable.setItem(i,1, QTableWidgetItem('%2.3f'%(self.initialNodeArray[i][0])))
			self.nodesTable.setItem(i,2, QTableWidgetItem('%2.3f'%(self.initialNodeArray[i][1])))	
			# Fix X Checkboxes
			if self.initialNodeArray[i][2] == True:
				self.nodesTable.setItem(i,3, QTableWidgetItem('\u2714'))
			else:
				self.nodesTable.setItem(i,3, QTableWidgetItem(''))
			# Fix Y Checkboxes
			if self.initialNodeArray[i][3] == True:
				self.nodesTable.setItem(i,4, QTableWidgetItem('\u2714'))
			else:
				self.nodesTable.setItem(i,4, QTableWidgetItem(''))
			# Reaction X Checkboxes
			if self.initialNodeArray[i][4] == True:
				self.nodesTable.setItem(i,5, QTableWidgetItem('\u2714'))
			else:
				self.nodesTable.setItem(i,5, QTableWidgetItem(''))
			# Reaction Y Checkboxes
			if self.initialNodeArray[i][5] == True:
				self.nodesTable.setItem(i,6, QTableWidgetItem('\u2714'))
			else:
				self.nodesTable.setItem(i,6, QTableWidgetItem(''))
	
	def redrawBeamTable(self):
		# Beam Table
		self.beamTable.setRowCount(len(self.initialBeamArray)+1)
		for i in range(0,len(self.initialBeamArray)):
			self.beamTable.setItem(i,0, QTableWidgetItem('%i'%(i+1)))
			self.beamTable.setItem(i,1, QTableWidgetItem('%i'%(self.initialBeamArray[i][0] + 1)))
			self.beamTable.setItem(i,2, QTableWidgetItem('%i'%(self.initialBeamArray[i][1] + 1)))	

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

	def redrawResultsLabels(self):
		#self.iterationLabel.setText("Iteration %i"%(self.damping))
		#self.length1Label.setText("Length 1 = %2.4f"%(self.mechanismLengths[1]))
		#self.length2Label.setText("Length 2 = %2.4f"%(self.mechanismLengths[2]))
		#self.length3Label.setText("Length 3 = %2.4f"%(self.mechanismLengths[3]))
		#self.length4Label.setText("Length 4 = %2.4f"%(self.mechanismLengths[4]))
		#self.length5Label.setText("Length 5 = %2.4f"%(self.mechanismLengths[5]))
		#self.base1Label = QLabel("Base 1 Location  = (%2.4f, %2.4f)"%(self.mechanismBases[0][0],self.mechanismBases[0][1]),self)
		#self.base2Label = QLabel("Base 2 Location  = (%2.4f, %2.4f)"%(self.mechanismBases[1][0],self.mechanismBases[1][1]),self)
		print('Redrawing Results')
	def parseDesign(self,design):
		print('going to parse design')

	def packageDesign(self):
		print('going to create design list')

	def iterationDone(self,design):
		print(design)

	def designOptimized(self,design):
		self.startButton.setEnabled(True)
		self.stopButton.setEnabled(False)
		self.optimizeThread.terminate
		print(design)

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
