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
import FourBarOptimize
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
        self.targets = [[0,5,10,1],[30,10,10,1],[60,10,5,1],[90,5,5,0]]
        self.exactSelected = 0
        self.countExactPoints()

        # Design:
        self.betas = [10,20]
        self.gammas = [15,25]
        self.design = [self.betas,self.gammas]

        # Synthesize Mechanism
        self.synthesizeMechanism()
        # Mechanism encoding: [lengths,base1,base2]
        self.initialMechanismBases = [[0.06,0.04],[1.5,-0.02]]
        baselength = self.vLen(self.initialMechanismBases[0],self.initialMechanismBases[1])
        #self.initialMechanismLengths = [baselength,0.3,0.75,0.6,1,1.5]
        self.initialMechanismLengths = [0.25, 1.1,0.83,1.25,1.46]
        self.initialAlpha = 0.75

        # Four Bar Properties
        self.mechanismBases = self.initialMechanismBases
        self.mechanismLengths = self.initialMechanismLengths
        self.mechanismPoints = self.calculateFourBarPoint(self.initialAlpha)

        # Four bar plotting
        self.plotAngles = np.linspace(0.0,6.28318530718,num=15).tolist()
        self.plotAngles.append(self.plotAngles[0])
        self.pathX = [0]*len(self.plotAngles)
        self.pathY = [0]*len(self.plotAngles)

        # UI Connections
        self.programLoaded = False
        self.dampingSlider.valueChanged.connect(self.dampingChanged)
        self.inputTable.itemSelectionChanged.connect(self.selectedInputTable)
        self.inputTable.cellChanged.connect(self.inputCellChanged)
        self.maxIterationsTextBox.textChanged.connect(self.maxIterationsChanged)
        self.redrawInputTables()

        # Display mechanims plot
        self.redraw_mechanism()

        self.programLoaded = True
        self.userEdited = True

    def startOptimization(self):
        # Create calculation thread
        self.design = self.packageDesign()
        controls = [self.damping,self.maxIterations]
        targets_to_send = self.targets
        for i in range(len(targets_to_send)):
            target = targets_to_send[i]
            target[0] = target[0]*2*np.pi/360
            targets_to_send[i] = target
        self.optimizeThread = FourBarOptimize.OptimizeThread(self.design,self.targets,controls)
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
                    self.synthesizeMechanism()
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
        self.calculateActualPath(self.plotAngles)
        path = [self.pathX,self.pathY]
        pose = self.calculateFourBarPoint(self.initialAlpha)
        self.graph_canvas.plotFourBar(self.targets, pose,path)


    def parseDesign(self,design):
        print('Do we need to parse the design?')

    def packageDesign(self):
        # Mechanism encoding: [lengths,base1,base2]
        return self.mechanismLengths + self.mechanismBases[0] + self.mechanismBases[1]

    def iterationDone(self,design):
        self.mechanismBases = [
            [design[5],design[6]],
            [design[7],design[8]]
        ]
        self.mechanismLengths = list(design[0:5])
        #self.progress = int(design[2]*100)
        print(design)
        self.redrawResults()
        self.resultsBar.setValue(self.progress)

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

    def synthesizeMechanism(self):
        precisionPoints = []
        for target in self.targets:
            if target[3] == True:
                precisionPoints.append(target)
        if len(precisionPoints) != 3:
            print("Missing precision points, synthesis terminating")
            return
        delta = [0,0]
        delta[0] = [precisionPoints[1][1] - precisionPoints[0][1],precisionPoints[1][2] - precisionPoints[0][2]]
        delta[1] = [precisionPoints[2][1] - precisionPoints[0][1],precisionPoints[2][2] - precisionPoints[0][2]]


    def calculateFourBarPoint(self,alphaInput):
        # lengths is in format [base, z1,z2,z3,z4,z5]
        # bases are in format [[base1X,base1Y],[base2X,base2Y]]
        # alpha is a single float for crank angle from horizontal

        bases = list(self.mechanismBases)
        length = list(self.mechanismLengths)
        length.insert(0,self.vLen(self.mechanismBases[0],self.mechanismBases[1]))
        alpha = alphaInput + self.mechanismStartAngle

        point1 = [0,0]
        point1[0] = bases[0][0] + length[1]*np.cos(alpha)
        point1[1] = bases[0][1] + length[1]*np.sin(alpha)

        baseRun = np.float64(bases[1][0] - bases[0][0])
        baseRise = np.float64(bases[1][1]-bases[0][1])
        #print(baseRun)
        #print(baseRise)
        baseAngle = 0
        if baseRise == 0:
            if baseRun > 0:
                baseAngle = 0
            elif baseRun < 0:
                baseAngle = np.pi
        elif baseRun == 0:
            if baseRise > 0:
                baseAngle = np.pi/2
            elif baseRise < 0:
                baseAngle = 3*np.pi/2
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

    def calculateActualPath(self,angleList):
        xActual = [0]*len(angleList)
        yActual = [0]*len(angleList)
        for i in range(0,len(angleList)):
            points = self.calculateFourBarPoint(angleList[i])
            xActual[i] = points[4][0]
            yActual[i] = points[4][1]
        self.pathX = xActual
        self.pathY = yActual

    def countExactPoints(self):
        self.exactSelected = 0
        for target in self.targets:
            if target[3] == True:
                self.exactSelected = self.exactSelected + 1

def main():
    app = QApplication(sys.argv)
    form = MainWindow()
    form.show()
    app.exec_()

if __name__ == '__main__':
    main()
