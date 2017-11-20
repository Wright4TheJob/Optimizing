# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'threading_design.ui'
#
# Created: Thu Aug  6 13:47:18 2015
#      by: PyQt4 UI code generator 4.10.4

from PyQt5 import QtCore, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtWidgets import (QWidget, QTreeView, QMessageBox, QHBoxLayout, 
							 QFileDialog, QLabel, QSlider, QCheckBox, 
							 QLineEdit, QVBoxLayout, QApplication, QPushButton,
							 QTableWidget, QTableWidgetItem,QSizePolicy,
							 QGridLayout,QGroupBox, QMainWindow,QAction,QHeaderView)
from PyQt5.QtCore import Qt, QTimer, QCoreApplication
from matplotlib.figure import Figure
from matplotlib import rcParams
import matplotlib.image as image
import math

try:
	_fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
	def _fromUtf8(s):
		return s

"""
try:
	_encoding = QtGui.QApplication.UnicodeUTF8
	def _translate(context, text, disambig):
		return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
	def _translate(context, text, disambig):
		return QtGui.QApplication.translate(context, text, disambig)
"""

rcParams.update({'figure.autolayout': True})

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
		self.axes.set_title('Truss')
			 
	def plotTruss(self,nodeList,beamList):
		# Nodelist Format: [X,Y,Fix X, Fix Y, Rx, Ry,Applied Force, Force Angle]
		self.axes.cla() #Clear axes
		
		# Plot roller symbol for constraints
		rollerSize = 0.1
		rollerForX = image.imread('image/RollerH.png')
		rollerForY = image.imread('image/RollerV.png')
		constraintLocation = [5,0]
		off = 0.05
		arrowLen = 0.5

		for i in range(0,len(nodeList)):
			if nodeList[i][2] == True: # X is constrained
				self.axes.imshow(rollerForX, extent=(nodeList[i][0]-2*rollerSize, nodeList[i][0], 
					nodeList[i][1] - rollerSize, nodeList[i][1] + rollerSize), zorder=2)
			if nodeList[i][3] == True:
				self.axes.imshow(rollerForY, extent=(nodeList[i][0]-rollerSize, nodeList[i][0] + rollerSize, 
					nodeList[i][1] - 2*rollerSize, nodeList[i][1]), zorder=-1)

			# Plot arrows for applied forces
			if nodeList[i][6] != 0:
				dx = arrowLen*math.cos(math.radians(nodeList[i][7]))
				dy = arrowLen*math.sin(math.radians(nodeList[i][7]))
				self.axes.arrow(nodeList[i][0], nodeList[i][1], dx, dy,color='r',zorder=3,shape='full',head_width=0.075, head_length=0.15)

			# Plot nodes
			self.axes.plot([nodeList[i][0]],[nodeList[i][1]],'ko')
			self.axes.text(nodeList[i][0]+off,nodeList[i][1]+off, '%i'%(i+1), fontsize=10)

			# Plot Reaction Forces
			if nodeList[i][4] == True: # X is constrained
				dx = -arrowLen/1.5
				dy = 0
				self.axes.arrow(nodeList[i][0]-dx, nodeList[i][1]-dy, dx, dy,color='g',
					length_includes_head = True,zorder=3,shape='full',head_width=0.075, head_length=0.15)				
			if nodeList[i][5] == True:
				dx = 0
				dy = arrowLen/1.5
				self.axes.arrow(nodeList[i][0]-dx, nodeList[i][1]-dy, dx, dy,color='g',
					length_includes_head = True,zorder=3,shape='full',head_width=0.075, head_length=0.15)				


		# Plot mechanism bars
		for i in range(0,len(beamList)):
			fromNode = beamList[i][0]
			toNode = beamList[i][1]
			self.axes.plot([nodeList[fromNode][0],nodeList[toNode][0]],[nodeList[fromNode][1],nodeList[toNode][1]],'k')
			midX = (nodeList[fromNode][0]+nodeList[toNode][0])/2
			midY = (nodeList[fromNode][1] + nodeList[toNode][1])/2
			self.axes.text(midX+off,midY+off, '%i'%(i+1), fontsize=10)

		#self.axes.set_xlabel(data_label)
		#self.axes.set_ylabel("Estimated Prob. Density Funct.")
		#self.axes.set_title(title)
		#self.axes.legend(shadow=True)
		self.axes.axis('equal')
		self.axes.margins(0.2, 0.2)
		self.draw()
		#print("Finished Drawing Normalized Histogram.")
		  

class Ui_MainWindow(object):
	def setupUi(self, MainWindow):
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
		MainWindow.setObjectName(_fromUtf8("MainWindow"))
		MainWindow.resize(1000, 800)
		self.centralwidget = QWidget(MainWindow)
		self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
		
	 
		### Controls box ###
		controlsBox = QGroupBox("Controls")
		controlsBoxLayout = QGridLayout()
		# Start Button
		self.startButton = QPushButton('Start',self)
		controlsBoxLayout.addWidget(self.startButton,0,0)
		# Stop Button
		self.stopButton = QPushButton('Stop',self)
		self.stopButton.setEnabled(False)
		controlsBoxLayout.addWidget(self.stopButton,0,1)
		# Damping Label and slider
		self.dampingLabel = QLabel("Damping = 1",self)
		controlsBoxLayout.addWidget(self.dampingLabel,1,0)
		self.dampingSlider = QSlider(Qt.Horizontal)
		self.dampingSlider.setMinimum(1)
		self.dampingSlider.setMaximum(1000)
		self.dampingSlider.setValue(1000)
		controlsBoxLayout.addWidget(self.dampingSlider,1,1)
		# Cross section selection dropdown menu
		# Max Iterations text box
		self.maxIterationsLabel = QLabel("Maximum Iterations",self)
		controlsBoxLayout.addWidget(self.maxIterationsLabel,2,0)
		self.maxIterationsTextBox = QLineEdit(self)
		self.maxIterationsTextBox.setText('100')
		controlsBoxLayout.addWidget(self.maxIterationsTextBox,2,1)
		# Max stress text box
		self.maxStressControlLabel = QLabel("Max Allowable Stress",self)
		controlsBoxLayout.addWidget(self.maxStressControlLabel,3,0)
		self.maxStressTextBox = QLineEdit(self)
		self.maxStressTextBox.setText('1')
		controlsBoxLayout.addWidget(self.maxStressTextBox,3,1)
		# Density optional text box
		self.densityLabel = QLabel("Density",self)
		controlsBoxLayout.addWidget(self.densityLabel,4,0)
		self.densityTextBox = QLineEdit(self)
		self.densityTextBox.setText('1')
		controlsBoxLayout.addWidget(self.densityTextBox,4,1)
		controlsBox.setLayout(controlsBoxLayout)

		# Results Labels
		#self.iterationLabel = QLabel("Iterations Not Started",self)
		#self.length1Label = QLabel("Length 1 Not Computed Yet",self)
		#self.length2Label = QLabel("Length 2 Not Computed Yet",self)
		#self.length3Label = QLabel("Length 3 Not Computed Yet",self)
		#self.length4Label = QLabel("Length 4 Not Computed Yet",self)
		#self.length5Label = QLabel("Length 5 Not Computed Yet",self)
		#self.base1Label = QLabel("Base 1 Location Not Computed Yet",self)
		#self.base2Label = QLabel("Base 2 Location Not Computed Yet",self)

		### Input Tables Box ###
		inputBox = QGroupBox('Input')
		inputBoxLayout = QGridLayout()
		# Node Table
		self.nodeTableLabel = QLabel("Enter Node Positions",self)
		self.nodeTableLabel.setAlignment(Qt.AlignCenter)
		inputBoxLayout.addWidget(self.nodeTableLabel,0,0,1,2)		
		self.nodesTable = QTableWidget()
		self.nodesTable.setColumnCount(7)
		self.nodesTable.setRowCount(1) # Make 1 longer than number of elements for manual addition of elements
		self.nodesTable.setHorizontalHeaderLabels(['Node','X','Y','Fix X','Fix Y','Reaction X','Reaction Y'])
		nodeHeader = self.nodesTable.horizontalHeader()
		nodeHeader.setSectionResizeMode(0, QHeaderView.ResizeToContents)
		nodeHeader.setSectionResizeMode(1, QHeaderView.Stretch)
		nodeHeader.setSectionResizeMode(2, QHeaderView.Stretch)
		nodeHeader.setSectionResizeMode(3, QHeaderView.ResizeToContents)
		nodeHeader.setSectionResizeMode(4, QHeaderView.ResizeToContents)
		nodeHeader.setSectionResizeMode(5, QHeaderView.ResizeToContents)
		nodeHeader.setSectionResizeMode(6, QHeaderView.ResizeToContents)
		inputBoxLayout.addWidget(self.nodesTable,1,0,1,2)
		# Beam Table
		self.beamTableLabel = QLabel("Enter Beam Connections",self)
		self.beamTableLabel.setAlignment(Qt.AlignCenter)
		inputBoxLayout.addWidget(self.beamTableLabel,2,0)		
		self.beamTable = QTableWidget()
		self.beamTable.setColumnCount(3)
		self.beamTable.setRowCount(1) # Make 1 longer than number of elements for manual addition of elements
		self.beamTable.setHorizontalHeaderLabels(['Beam','From Node','To Node'])
		beamHeader = self.beamTable.horizontalHeader()
		beamHeader.setSectionResizeMode(0, QHeaderView.ResizeToContents)
		beamHeader.setSectionResizeMode(1, QHeaderView.Stretch)
		beamHeader.setSectionResizeMode(2, QHeaderView.Stretch)
		inputBoxLayout.addWidget(self.beamTable,3,0)
		# External Force Table
		self.forceTableLabel = QLabel("Enter Applied Forces",self)
		self.forceTableLabel.setAlignment(Qt.AlignCenter)
		inputBoxLayout.addWidget(self.forceTableLabel,2,1)		
		self.forceTable = QTableWidget()
		self.forceTable.setColumnCount(3)
		self.forceTable.setRowCount(1) # Make 1 longer than number of elements for manual addition of elements
		self.forceTable.setHorizontalHeaderLabels(['Node','Force','Angle'])
		forceTableHeader = self.forceTable.horizontalHeader()
		forceTableHeader.setSectionResizeMode(0, QHeaderView.ResizeToContents)
		forceTableHeader.setSectionResizeMode(1, QHeaderView.Stretch)
		forceTableHeader.setSectionResizeMode(2, QHeaderView.Stretch)
		inputBoxLayout.addWidget(self.forceTable,3,1)
		inputBox.setLayout(inputBoxLayout)

		# Plot
		self.graph_canvas = MyDynamicMplCanvas(self.centralwidget, width=5, height=4, dpi=120)
				
		# Results Label Box 
		resultsBox = QGroupBox("Results")
		resultsBoxLayout = QVBoxLayout()
		#resultsBoxLayout.addWidget(self.iterationLabel)       
		#resultsBoxLayout.addWidget(self.length1Label)
		#resultsBoxLayout.addWidget(self.length2Label)
		#resultsBoxLayout.addWidget(self.length3Label)
		#resultsBoxLayout.addWidget(self.length4Label)
		#resultsBoxLayout.addWidget(self.length5Label)
		#resultsBoxLayout.addWidget(self.base1Label)
		#resultsBoxLayout.addWidget(self.base2Label)
		resultsBox.setLayout(resultsBoxLayout)
				
		#Now we can set all the previously defined boxes into the main window
		master_layout = QGridLayout()
		master_layout.addWidget(inputBox,0,0) 
		master_layout.addWidget(resultsBox,1,1)
		master_layout.addWidget(controlsBox,1,0)
		master_layout.addWidget(self.graph_canvas,0,1) 
		#master_layout.addWidget(distribution_box,1,1)
		
		#self.centralwidget.addWidget(master_layout)
		self.centralwidget.setLayout(master_layout)
		
		self.setWindowTitle('Four Bar Linkage Optimization')
		self.activateWindow()
		self.raise_()
		self.show()
		MainWindow.setCentralWidget(self.centralwidget)

		menuBar = self.menuBar()

		file_menu = menuBar.addMenu('&File')

		open_file = QAction('&Open', self)
		open_file.setShortcut('Ctrl+O')
		open_file.setStatusTip('Load Truss Design')
		open_file.triggered.connect(self.load_data)
		file_menu.addAction(open_file)

		save_file = QAction('&Save',self)
		save_file.setShortcut('Ctrl+S')
		save_file.setStatusTip('Save Optimized Design')
		save_file.triggered.connect(self.save_data)
		file_menu.addAction(save_file)

		exit_action = QAction('&Exit', self)        
		exit_action.setShortcut('Ctrl+Q')
		exit_action.setStatusTip('Exit application')
		exit_action.triggered.connect(self.close) #This is built in
		file_menu.addAction(exit_action)

		QtCore.QMetaObject.connectSlotsByName(MainWindow)
