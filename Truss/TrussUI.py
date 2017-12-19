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
							 QGridLayout,QGroupBox, QMainWindow,QAction,QHeaderView,QComboBox,QProgressBar)
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
				self.axes.imshow(rollerForX, extent=(nodeList[i][0]-2*rollerSize, nodeList[i][0], nodeList[i][1] - rollerSize, nodeList[i][1] + rollerSize), zorder=2)
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
			if (fromNode != -1 and toNode != -1):
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
		self.dampingLabel = QLabel("Damping = 0.1",self)
		controlsBoxLayout.addWidget(self.dampingLabel,1,0)
		self.dampingSlider = QSlider(Qt.Horizontal)
		self.dampingSlider.setMinimum(1)
		self.dampingSlider.setMaximum(1000)
		self.dampingSlider.setValue(100)
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
		self.maxStressTextBox.setText('10')
		controlsBoxLayout.addWidget(self.maxStressTextBox,3,1)
		# Density optional text box
		self.densityLabel = QLabel("Density",self)
		controlsBoxLayout.addWidget(self.densityLabel,4,0)
		self.densityTextBox = QLineEdit(self)
		self.densityTextBox.setText('1')
		controlsBoxLayout.addWidget(self.densityTextBox,4,1)
		self.crossSectionLabel = QLabel("Cross Section",self)
		controlsBoxLayout.addWidget(self.crossSectionLabel,5,0)
		self.crossSectionBox = QComboBox(self)
		self.crossSectionBox.addItem("Rectangular - Equal Thickness")
		self.crossSectionBox.addItem("Rectangular")
		self.crossSectionBox.addItem("Rectangular - Hollow")
		self.crossSectionBox.addItem("Square")
		self.crossSectionBox.addItem("Square - Hollow")
		self.crossSectionBox.addItem("Round")
		self.crossSectionBox.addItem("Round - Hollow")
		self.crossSectionBox.activated[str].connect(self.crossSectionChanged)

		controlsBoxLayout.addWidget(self.crossSectionBox,5,1)

		controlsBox.setLayout(controlsBoxLayout)

		### Input Tables Box ###
		inputBox = QGroupBox('Input')
		inputBoxLayout = QGridLayout()
		# Node Table
		self.nodeTableLabel = QLabel("Enter Node Positions",self)
		self.nodeTableLabel.setAlignment(Qt.AlignCenter)
		inputBoxLayout.addWidget(self.nodeTableLabel,0,0,1,2)		
		self.nodesTable = QTableWidget()
		self.nodesTable.setColumnCount(6)
		self.nodesTable.setRowCount(1) # Make 1 longer than number of elements for manual addition of elements
		self.nodesTable.setHorizontalHeaderLabels(['X','Y','Fix X','Fix Y','Reaction X','Reaction Y'])
		nodeHeader = self.nodesTable.horizontalHeader()
		nodeHeader.setSectionResizeMode(0, QHeaderView.Stretch)
		nodeHeader.setSectionResizeMode(1, QHeaderView.Stretch)
		nodeHeader.setSectionResizeMode(2, QHeaderView.ResizeToContents)
		nodeHeader.setSectionResizeMode(3, QHeaderView.ResizeToContents)
		nodeHeader.setSectionResizeMode(4, QHeaderView.ResizeToContents)
		nodeHeader.setSectionResizeMode(5, QHeaderView.ResizeToContents)
		inputBoxLayout.addWidget(self.nodesTable,1,0,1,2)
		# Beam Table
		self.beamTableLabel = QLabel("Enter Beam Connections",self)
		self.beamTableLabel.setAlignment(Qt.AlignCenter)
		inputBoxLayout.addWidget(self.beamTableLabel,2,0)		
		self.beamTable = QTableWidget()
		self.beamTable.setColumnCount(2)
		self.beamTable.setRowCount(1) # Make 1 longer than number of elements for manual addition of elements
		self.beamTable.setHorizontalHeaderLabels(['From Node','To Node'])
		beamHeader = self.beamTable.horizontalHeader()
		beamHeader.setSectionResizeMode(0, QHeaderView.Stretch)
		beamHeader.setSectionResizeMode(1, QHeaderView.Stretch)
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
				
		### Results Tables Box ###
		resultsBox = QGroupBox("Results")
		resultsBoxLayout = QGridLayout()
		self.resultsBarLabel = QLabel("Optimization Progress: ",self)
		resultsBoxLayout.addWidget(self.resultsBarLabel,0,0)		

		self.resultsBar = QProgressBar(self)
		resultsBoxLayout.addWidget(self.resultsBar,0,1)		
		# Node Table
		self.nodeResultsTableLabel = QLabel("Optimized Node Positions",self)
		self.nodeResultsTableLabel.setAlignment(Qt.AlignCenter)
		resultsBoxLayout.addWidget(self.nodeResultsTableLabel,1,0)		
		self.nodesResultsTable = QTableWidget()
		self.nodesResultsTable.setColumnCount(3)
		self.nodesResultsTable.setRowCount(1) # Make 1 longer than number of elements for manual addition of elements
		self.nodesResultsTable.setHorizontalHeaderLabels(['Node','X','Y'])
		nodeResultsHeader = self.nodesResultsTable.horizontalHeader()
		nodeResultsHeader.setSectionResizeMode(0, QHeaderView.ResizeToContents)
		nodeResultsHeader.setSectionResizeMode(1, QHeaderView.Stretch)
		nodeResultsHeader.setSectionResizeMode(2, QHeaderView.Stretch)
		resultsBoxLayout.addWidget(self.nodesResultsTable,2,0)
		# Beam Table
		self.beamResultsTableLabel = QLabel("Optimized Beam Properties",self)
		self.beamResultsTableLabel.setAlignment(Qt.AlignCenter)
		resultsBoxLayout.addWidget(self.beamResultsTableLabel,1,1)		
		self.beamResultsTable = QTableWidget()
		self.beamResultsTable.setColumnCount(4)
		self.beamResultsTable.setRowCount(1) # Make 1 longer than number of elements for manual addition of elements
		self.beamResultsTable.setHorizontalHeaderLabels(['Length','OD', 'ID', 'Stress'])
		beamResultsHeader = self.beamResultsTable.horizontalHeader()
		beamResultsHeader.setSectionResizeMode(0, QHeaderView.Stretch)
		beamResultsHeader.setSectionResizeMode(1, QHeaderView.Stretch)
		beamResultsHeader.setSectionResizeMode(2, QHeaderView.Stretch)
		beamResultsHeader.setSectionResizeMode(3, QHeaderView.Stretch)
		resultsBoxLayout.addWidget(self.beamResultsTable,2,1)		
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

		saveInput_file = QAction('&Save Input Design',self)
		saveInput_file.setStatusTip('Save Optimized Design')
		saveInput_file.triggered.connect(self.saveInputData)
		file_menu.addAction(saveInput_file)

		saveOptimized_file = QAction('&Save Optimized Design',self)
		saveOptimized_file.setShortcut('Ctrl+S')
		saveOptimized_file.setStatusTip('Save Optimized Design')
		saveOptimized_file.triggered.connect(self.saveOptimizedData)
		file_menu.addAction(saveOptimized_file)

		exit_action = QAction('&Exit', self)        
		exit_action.setShortcut('Ctrl+Q')
		exit_action.setStatusTip('Exit application')
		exit_action.triggered.connect(self.close) #This is built in
		file_menu.addAction(exit_action)

		QtCore.QMetaObject.connectSlotsByName(MainWindow)
