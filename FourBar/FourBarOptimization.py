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

rcParams.update({'figure.autolayout': True})

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        menuBar = self.menuBar()

        file_menu = menuBar.addMenu('&File')

        open_file = QAction(None, '&Open', self)
        open_file.setShortcut('Ctrl+O')
        open_file.setStatusTip('Load Counter Data')
        open_file.triggered.connect(self.load_data)
        file_menu.addAction(open_file)

        save_file = QAction(None, '&Save',self)
        save_file.setShortcut('Ctrl+S')
        save_file.setStatusTip('Save Counter Data')
        save_file.triggered.connect(self.save_data)
        file_menu.addAction(save_file)

        exit_action = QAction(None, '&Exit', self)
        exit_action.setShortcut('Ctrl+Q')
        exit_action.setStatusTip('Exit application')
        exit_action.triggered.connect(self.close) #This is built in
        file_menu.addAction(exit_action)

        self.setCentralWidget(self.mainWidget)

class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        FigureCanvas.mpl_connect(self,'button_press_event', self.double_click)

    def export(self,event):
        filename = "ExportedGraph.pdf"
        self.fig.savefig(filename)
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setText("Saved a copy of the graphics window to {}".format(filename))
        #msg.setInformativeText("This is additional information")
        msg.setWindowTitle("Saved PDF File")
        msg.setDetailedText("The full path of the file is \n{}".format(os.path.abspath(os.getcwd())))
        msg.setStandardButtons(QMessageBox.Ok)
        msg.setWindowModality(Qt.ApplicationModal)
        msg.exec_()
        print("Exported PDF file")

    def double_click(self, event):
        FigureCanvas.mpl_connect(self,'button_press_event', self.export)

class MyDynamicMplCanvas(MyMplCanvas):
    """A canvas that updates itself frequently with a new plot."""
    def __init__(self, *args, **kwargs):
        MyMplCanvas.__init__(self, *args, **kwargs)
        self.axes.set_xlabel("X")
        self.axes.set_ylabel("Y")
        self.axes.set_title("Four Bar Linkage")

    def plotFourBar(self,targets,mechPoints,path):
        self.axes.cla() #Clear axes
        # Plot mechanism bars
        # (crank) Base1 to point 3
        self.axes.plot([mechPoints[0][0],mechPoints[2][0]],[mechPoints[0][1],mechPoints[2][1]],'k')
        # (Rocker) Base 2 to point 4
        self.axes.plot([mechPoints[1][0],mechPoints[3][0]],[mechPoints[1][1],mechPoints[3][1]],'k')
        # (DyadBase) Point 3 to point 4
        self.axes.plot([mechPoints[2][0],mechPoints[3][0]],[mechPoints[2][1],mechPoints[3][1]],'k')
        # (Dyad Side 1) Point 3 to point 5
        self.axes.plot([mechPoints[2][0],mechPoints[4][0]],[mechPoints[2][1],mechPoints[4][1]],'k')
        # (Dyad Side 2) Point 4 to point 5
        self.axes.plot([mechPoints[3][0],mechPoints[4][0]],[mechPoints[3][1],mechPoints[4][1]],'k')

        # Plot target positons
        xTargets = targets[1]
        yTargets = targets[2]
        self.axes.plot(xTargets, yTargets, 'rx')
        # Plot actual path
        #print(path[0])
        self.axes.plot(path[0],path[1],'0.45')

        #self.axes.set_xlabel(data_label)
        #self.axes.set_ylabel("Estimated Prob. Density Funct.")
        #self.axes.set_title(title)
        #self.axes.legend(shadow=True)
        self.draw()
        #print("Finished Drawing Normalized Histogram.")


class FourBarOptimizer(QWidget):

    def __init__(self):
        super().__init__()

        # Upon startup, run a user interface routine
        self.init_ui()

    def init_ui(self):
        self.damping = 0.1
        self.iterations = 0
        self.maxIterations = 1
        self.mechanismStartAngle = 0.017

        #self.thetaTargets = [0,90,180,270]
        #self.xTargets = [1,4,6,3]
        #self.yTargets = [1,0,1,2]
        self.thetaTargets = [0,90,180,270]
        self.xTargets = [1,3,1,2]
        self.yTargets = [1,2,2,1]
        self.betas = [0,10,20]
        self.gammas = [5,15,25]
        self.targets = [self.thetaTargets,self.xTargets,self.yTargets]
        self.exact = [True,True,False,True]
        self.exactSelected = sum(self.exact)
        self.fixedPointCount = min(len(self.thetaTargets),len(self.xTargets),len(self.yTargets))
        # Mechanism encoding: [[lengths],[base1],[base2],]
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
        self.setGeometry(200,30,1000,700)

        self.load_button = QPushButton('Load Data',self)
        self.load_button.clicked.connect(self.load_data)

        self.runButton = QPushButton('Run Optimization',self)
        self.runButton.clicked.connect(self.runFourBarOptimization)

        self.dampingSlider = QSlider(Qt.Horizontal)
        self.dampingSlider.setMinimum(1)
        self.dampingSlider.setMaximum(1000)
        #self.dampingSlider.setTickInterval(10)
        #self.dampingSlider.setSingleStep(0.01)
        self.dampingSlider.setValue(self.damping*1000)
        self.dampingSlider.valueChanged.connect(self.dampingChanged)
        self.dampingLabel = QLabel("Damping = %1.2f"%(self.damping),self)

        # Results Labels
        self.iterationLabel = QLabel("Iterations Not Started",self)
        self.length1Label = QLabel("Length 1 Not Computed Yet",self)
        self.length2Label = QLabel("Length 2 Not Computed Yet",self)
        self.length3Label = QLabel("Length 3 Not Computed Yet",self)
        self.length4Label = QLabel("Length 4 Not Computed Yet",self)
        self.length5Label = QLabel("Length 5 Not Computed Yet",self)
        self.base1Label = QLabel("Base 1 Location Not Computed Yet",self)
        self.base2Label = QLabel("Base 2 Location Not Computed Yet",self)

        #Set up a Table to display data
        self.data_table = QTableWidget()
        self.data_table.itemSelectionChanged.connect(self.selectedTableItem)
        self.data_table.cellChanged.connect(self.cellChanged)
        self.data_table.setColumnCount(4)
        self.data_table.setRowCount(len(self.thetaTargets)+1)
        self.data_table.setHorizontalHeaderLabels(['Exact Match?','Crank Angle','X','Y'])



        self.main_widget = QWidget(self)
        self.graph_canvas = MyDynamicMplCanvas(self.main_widget, width=5, height=4, dpi=120)

        targets = [self.thetaTargets,self.xTargets,self.yTargets]
        self.calculateActualPath(self.plotAngles)
        path = [self.pathX,self.pathY]
        mechanismPoints = self.calculateFourBarPoint(self.initialAlpha)
        self.graph_canvas.plotFourBar(targets, mechanismPoints,path)

        #Define where the widgets go in the window
        #We start by defining some boxes that we can arrange

        #Create a GUI box to put all the table and data widgets in
        table_box = QGroupBox("Data Table")
        #Create a layout for that box using the vertical
        table_box_layout = QVBoxLayout()
        #Add the widgets into the layout
        table_box_layout.addWidget(self.load_button)
        table_box_layout.addWidget(self.data_table)

        #setup the layout to be displayed in the box
        table_box.setLayout(table_box_layout)

        # Results Label Box
        resultsBox = QGroupBox("Results")
        resultsBoxLayout = QVBoxLayout()
        resultsBoxLayout.addWidget(self.iterationLabel)
        resultsBoxLayout.addWidget(self.length1Label)
        resultsBoxLayout.addWidget(self.length2Label)
        resultsBoxLayout.addWidget(self.length3Label)
        resultsBoxLayout.addWidget(self.length4Label)
        resultsBoxLayout.addWidget(self.length5Label)
        resultsBoxLayout.addWidget(self.base1Label)
        resultsBoxLayout.addWidget(self.base2Label)
        resultsBox.setLayout(resultsBoxLayout)

        # Controls Box
        controlsBox = QGroupBox("Controls")
        controlsBoxLayout = QVBoxLayout()
        controlsBoxLayout.addWidget(self.runButton)
        controlsBoxLayout.addWidget(self.dampingLabel)
        controlsBoxLayout.addWidget(self.dampingSlider)
        controlsBox.setLayout(controlsBoxLayout)

        #Now we can set all the previously defined boxes into the main window
        grid_layout = QGridLayout()
        grid_layout.addWidget(table_box,0,0)
        grid_layout.addWidget(resultsBox,1,0)
        grid_layout.addWidget(controlsBox,1,1)
        grid_layout.addWidget(self.graph_canvas,0,1)
        #grid_layout.addWidget(distribution_box,1,1)

        self.redrawTable()
        self.setLayout(grid_layout)

        self.setWindowTitle('Four Bar Linkage Optimization')
        self.activateWindow()
        self.raise_()
        self.show()

    def synthesizeFourBar(self,betas,gammas):
        exactPoints = []
        for i in range(0,len(thetaTargets)):
            if exact[i] == True:
                exactPoints.append([thetaTargets[i],xTargets[i],yTargets[i]])
        delta2 = (exactPoints[1][1]-exactPoints[0][1]) + (exactPoints[1][2]-exactPoints[0][2])j
        delta3 = (exactPoints[2][1]-exactPoints[0][1]) + (exactPoints[2][2]-exactPoints[0][2])j

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
            f.write("%i" % (self.counter_value))
            f.close()

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
        baseAngle = atan((b2y-b1y)/(b2x-b1x))

        x3 = l0*cos(baseAngle) - l1*cos(alpha + a0)
        y3 = l0*sin(baseAngle) - l1*sin(alpha + a0)
        theta3ArccosValue = (x3**2 + y3**2 + (-l3)**2 - l5**2)/2 * (-l3) * sqrt(x3**2 + y3**2)
        theta3 = atan2(y3,x3) + acos(theta3ArccosValue)
        theta5 = atan2((y3-(-l3)*sin(theta3))/l5, (x3-(-l3)*cos(theta3))/l5)

        dyadAngle1 = acos((l5**2 + l2**2 - l4**2)/(2*l5*l2))
        Cx = b1x+l1*cos(alpha + a0) + l2*cos(theta5 + dyadAngle1)
        Cy = b1y+l1*sin(alpha + a0) + l2*sin(theta5 + dyadAngle1)

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
        #print(baseRun)
        #print(baseRise)
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
        # Calculate betas and gammas for target design - 3 points
        # Evaluate error from target path
        #
        print('Starting Optimization')
        (xExpression,yExpression) = self.fourBarExpression()
        objectiveFunction = 0

        for i in range(0,len(self.xTargets)):
            objectiveExpression = (self.xTargets[i] - xExpression)**2 + (self.yTargets[i] - yExpression)**2
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
        rp = 0.5
        testSteps = [0,0.01,0.02]
        expression = objectiveFunction
        epsilon = 0.001
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
            deltaObjective = abs(float(abs(objectiveValueLast) - abs(objectiveValue)))
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

    def dampingChanged(self,value):
        self.damping = float(value)/1000
        self.dampingLabel.setText("Damping = %1.2f"%(self.damping))

    def cellChanged(self,row,column):
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
       self.data_table.setRowCount(len(self.thetaTargets)+1)
       for i in range(0,len(self.thetaTargets)):
            if self.exact[i] == True:
                self.data_table.setItem(i,0, QTableWidgetItem('\u2714'))
            else:
                self.data_table.setItem(i,0, QTableWidgetItem(''))

            self.data_table.setItem(i,1, QTableWidgetItem('%2.4f'%(self.thetaTargets[i])))
            self.data_table.setItem(i,2, QTableWidgetItem('%2.4f'%(self.xTargets[i])))
            self.data_table.setItem(i,3, QTableWidgetItem('%2.4f'%(self.yTargets[i])))

    def redrawResultsLabels(self):
        self.iterationLabel.setText("Iteration %i"%(self.damping))
        self.length1Label.setText("Length 1 = %2.4f"%(self.mechanismLengths[1]))
        self.length2Label.setText("Length 2 = %2.4f"%(self.mechanismLengths[2]))
        self.length3Label.setText("Length 3 = %2.4f"%(self.mechanismLengths[3]))
        self.length4Label.setText("Length 4 = %2.4f"%(self.mechanismLengths[4]))
        self.length5Label.setText("Length 5 = %2.4f"%(self.mechanismLengths[5]))
        self.base1Label.setText("Base 1 Location  = (%2.4f, %2.4f)"%(self.mechanismBases[0][0],self.mechanismBases[0][1]))
        self.base2Label.setText("Base 2 Location  = (%2.4f, %2.4f)"%(self.mechanismBases[1][0],self.mechanismBases[1][1]))

if __name__ == '__main__':
    #Start the program this way according to https://stackoverflow.com/questions/40094086/python-kernel-dies-for-second-run-of-pyqt5-gui
    app = QCoreApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    execute = FourBarOptimizer()
    sys.exit(app.exec_())
