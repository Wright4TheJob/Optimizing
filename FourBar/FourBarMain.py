#!/bin/env/python
# David Wright
# Copyright 2017
# Written for Python 3.5.2

# TODO: Reset view to input parameters
# TODO: Add cross sections to design variables
# TODO: Calculate euler buckling bounds
# TODO: Export to OpenSCAD?
# TODO: Add beam weights as forces to truss nodes (Can it hold itself up?)
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
import FourBarUI
import FourBarOptimize as fb
from scipy import stats
from sympy import *
import numpy as np
import os
import cmath
import math
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import rcParams
import matplotlib.mlab as mlab
import csv


class MainWindow(QMainWindow, FourBarUI.Ui_MainWindow):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        self.startButton.clicked.connect(self.startOptimization)
        # Initialize variables
        self.damping = 0.1
        self.progress = 0
        self.maxIterations = 1
        self.mechanismStartAngle = 0.017

        # Targets = [[theta,x,y,exct],[theta,x,y,exact]]
        self.targets = [[0,.2,.2,1],[90,0.2,2,1],[180,2,2,1],[270,2,.5,0]]
        self.exactSelected = 0
        self.countExactPoints()

        # Design:
        self.betas = [10,20]
        self.gammas = [15,25]
        self.design = [self.betas,self.gammas]

        # Mechanism encoding: [lengths,base1,base2]
        self.initialMechanismBases = [[0.06,0.04],[1.5,-0.02]]
        baselength = fb.vLen(self.initialMechanismBases[0],self.initialMechanismBases[1])
        #self.initialMechanismLengths = [baselength,0.3,0.75,0.6,1,1.5]
        self.initialMechanismLengths = [.8, 1.2,1.1,1,1]
        self.initialAlpha = 1.5

        # Four Bar Properties
        self.mechanismBases = self.initialMechanismBases
        self.mechanismLengths = self.initialMechanismLengths
        self.mechanism = [self.mechanismBases, self.mechanismLengths]
        self.mechanismPoints = fb.calculateFourBarPoint(self.mechanism, self.initialAlpha)

        # Four bar plotting
        self.plotAngles = np.linspace(0.0,2*np.pi,num=360).tolist()
        self.plotAngles.append(self.plotAngles[0])
        self.pathX = [0]*len(self.plotAngles)
        self.pathY = [0]*len(self.plotAngles)

        # UI Connections
        self.programLoaded = False
        self.dampingSlider.valueChanged.connect(self.dampingChanged)
        self.angleSlider.valueChanged.connect(self.angleChanged)
        self.inputTable.itemSelectionChanged.connect(self.selectedInputTable)
        self.inputTable.cellChanged.connect(self.inputCellChanged)
        self.maxIterationsTextBox.textChanged.connect(self.maxIterationsChanged)
        self.redrawInputTables()

        # Display mechanims plot
        self.redraw_mechanism()
        self.redrawResultsLabels()
        self.programLoaded = True
        self.userEdited = True

    def startOptimization(self):
        # Create calculation thread
        controls = [self.damping,self.maxIterations]
        targets_to_send = self.targets
        for i in range(len(targets_to_send)):
            target = targets_to_send[i]
            target[0] = target[0]*np.pi/18
            targets_to_send[i] = target
        self.optimizeThread = fb.OptimizeThread(self.mechanism,targets_to_send,controls)
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

    def maxIterationsChanged(self,newText):
        if newText != "":
            self.maxIterations = int(newText)
        else:
            self.maxIterations = 1

    def selectedInputTable(self):
        for currentQTableWidgetItem in self.inputTable.selectedItems():
            row = currentQTableWidgetItem.row()
            col = currentQTableWidgetItem.column()
            # Invert chekcboxes and save for column 3
            if (col == 3):
                cellValue = self.targets[row][col]
                if cellValue == True:
                    self.targets[row][col] = not self.targets[row][col]
                else:
                    if self.exactSelected < 3:
                        self.targets[row][col] = not self.targets[row][col]

                self.userEdited = False
                if self.targets[row][col] == True:
                    self.inputTable.setItem(row,col, QTableWidgetItem('\u2714'))
                else:
                    self.inputTable.setItem(row,col, QTableWidgetItem(' '))
                self.userEdited = True
                #self.graph_canvas.plotTruss(self.initialNodeArray,self.initialBeamArray)
                self.countExactPoints()
                if self.exactSelected < 3:
                    self.startButton.setEnabled(False)
                if self.exactSelected == 3:
                    self.startButton.setEnabled(True)

            # Edit values and save for columns 0 and 1
            else:
                self.inputTable.editItem(currentQTableWidgetItem)

        self.inputTable.clearSelection()
        #self.redrawNodeTable()

    def inputCellChanged(self,row,column):
        if self.programLoaded == True:
            if self.userEdited == True:
                self.userEdited = False

                # Add a new entry to the nodes list if last line is selected
                if row == len(self.targets):
                    cell = self.inputTable.item(row, column)
                    cellText = cell.text()
                    if cellText != '':
                        cellValue = float(cellText)
                    else:
                        cellValue = 0

                    if column == 1:
                        self.targets.append([cellValue,0,0,0])
                    elif column == 2:
                        self.targets.append([0,cellValue,0,0])
                else:
                    # Grab float value of text input for columns 0 and 1
                    if (column == 0 or column == 1 or column == 2):
                        cell = self.inputTable.item(row, column)
                        cellText = cell.text()
                        if cellText != '':
                            cellValue = float(cellText)
                        else:
                            cellValue = 0
                        self.targets[row][column] = cellValue
                    # SelectedNodesTable already took care of checkmark assignment
                self.inputTable.setRowCount(len(self.targets)+1)
                self.userEdited = True
                path = [self.pathX,self.pathY]
                self.graph_canvas.plotFourBar(self.targets, self.mechanismPoints, path)

    def redrawInputTables(self):
        self.userEdited = False
        self.redrawInputTable()
        self.userEdited = True

    def redrawInputTable(self):
        # Node Table
        self.inputTable.setRowCount(len(self.targets)+1)
        for i in range(0,len(self.targets)):
            self.inputTable.setItem(i,0, QTableWidgetItem('%2.3f'%(self.targets[i][0])))
            self.inputTable.setItem(i,1, QTableWidgetItem('%2.3f'%(self.targets[i][1])))
            self.inputTable.setItem(i,2, QTableWidgetItem('%2.3f'%(self.targets[i][2])))
            # Exact Checkboxes
            if self.targets[i][3] == True:
                self.inputTable.setItem(i,3, QTableWidgetItem('\u2714'))
            else:
                self.inputTable.setItem(i,3, QTableWidgetItem(''))

    def redrawResults(self):
        self.redrawResultsLabels()
        self.redraw_mechanism()

    def redrawResultsLabels(self):
        #self.iterationLabel.setText("Iteration %i"%(self.damping))
        self.length1Label.setText("Length 1 = %2.4f"%(self.mechanismLengths[0]))
        self.length2Label.setText("Length 2 = %2.4f"%(self.mechanismLengths[1]))
        self.length3Label.setText("Length 3 = %2.4f"%(self.mechanismLengths[2]))
        self.length4Label.setText("Length 4 = %2.4f"%(self.mechanismLengths[3]))
        self.length5Label.setText("Length 5 = %2.4f"%(self.mechanismLengths[4]))
        self.base1Label.setText("Base 1 Location  = (%2.4f, %2.4f)"%(self.mechanismBases[0][0],self.mechanismBases[0][1]))
        self.base2Label.setText("Base 2 Location  = (%2.4f, %2.4f)"%(self.mechanismBases[1][0],self.mechanismBases[1][1]))

    def redraw_mechanism(self):
        # Display mechanims plot
        print("Redrawing Mechanism")
        print("Bases:")
        print(self.mechanismBases)
        print("Lengths:")
        print(self.mechanismLengths)
        self.calculateActualPath(self.plotAngles)
        path = [self.pathX,self.pathY]
        try:
            pose = fb.calculateFourBarPose([self.mechanismBases, self.mechanismLengths],self.initialAlpha)
        except ValueError:
            pose = None
            print("Invalid pose")
        self.graph_canvas.plotFourBar(self.targets, pose,path)

    def iterationDone(self,design):
        print("Iteration Done")
        self.mechanismBases = design[0]
        self.mechanismLengths = design[1]
        self.mechanism = [self.mechanismBases, self.mechanismLengths]

        #self.progress = int(design[2]*100)
        #print(design)
        self.redrawResults()
        self.resultsBar.setValue(self.progress)

    def designOptimized(self,design):
        self.optimizeThread.terminate
        self.startButton.setEnabled(True)
        self.stopButton.setEnabled(False)
        #print(design)

    def calculateActualPath(self,angleList):
        xActual = [0]*len(angleList)
        yActual = [0]*len(angleList)
        for i in range(0,len(angleList)):
            try:
                endpoint = fb.calculateFourBarPoint([self.mechanismBases, self.mechanismLengths],angleList[i])
                xActual[i] = endpoint[0]
                yActual[i] = endpoint[1]
            except ValueError:
                print("Mechanism not solvable")
        self.pathX = xActual
        self.pathY = yActual

    def countExactPoints(self):
        self.exactSelected = 0
        for target in self.targets:
            if target[3] == True:
                self.exactSelected = self.exactSelected + 1

    def angleChanged(self, angle):
        self.initialAlpha = angle * np.pi/ 180
        self.angleLabel.setText("Angle = %i"%(self.initialAlpha*180 / np.pi))
        try:
            pose = fb.calculateFourBarPose([self.mechanismBases, self.mechanismLengths],self.initialAlpha)
        except ValueError:
            pose = fb.dummy_pose([self.mechanismBases, self.mechanismLengths], self.initialAlpha)
            print("Invalid pose")
        self.graph_canvas.plotFourBar(self.targets, pose, [self.pathX,self.pathY])

def main():
    app = QApplication(sys.argv)
    form = MainWindow()
    form.show()
    app.exec_()

if __name__ == '__main__':
    main()
