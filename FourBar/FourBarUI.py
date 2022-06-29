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
		self.axes.set_title("Four Bar Linkage")

	def plotFourBar(self,targets,mechPoints,path):
		self.axes.cla() #Clear axes
		# Plot mechanism bars
		# (crank) Base1 to point 3
		if mechPoints:
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
		for target in targets:
			self.axes.plot(target[1], target[2], 'rx')
		# Plot actual path
		#print(path[0])
		self.axes.plot(path[0],path[1],'0.45')

		#self.axes.set_xlabel(data_label)
		#self.axes.set_ylabel("Estimated Prob. Density Funct.")
		#self.axes.set_title(title)
		#self.axes.legend(shadow=True)
		self.draw()
		#print("Finished Drawing Normalized Histogram.")

	def plot_line(self, a, b):
		self.axes.plot([a[0],b[0]],[a[1],b[1]],'k')

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
		controlsBox.setLayout(controlsBoxLayout)


		### Input Tables Box ###
		inputBox = QGroupBox('Input')
		inputBoxLayout = QGridLayout()
		# Node Table
		self.inputTableLabel = QLabel("Enter Target Path Points",self)
		self.inputTableLabel.setAlignment(Qt.AlignCenter)
		inputBoxLayout.addWidget(self.inputTableLabel,0,0)
		self.inputTable = QTableWidget()
		self.inputTable.setColumnCount(4)
		self.inputTable.setRowCount(1) # Make 1 longer than number of elements for manual addition of elements
		self.inputTable.setHorizontalHeaderLabels(['Crank Angle','X','Y','Exact'])
		inputHeader = self.inputTable.horizontalHeader()
		inputHeader.setSectionResizeMode(0, QHeaderView.Stretch)
		inputHeader.setSectionResizeMode(1, QHeaderView.Stretch)
		inputHeader.setSectionResizeMode(2, QHeaderView.Stretch)
		inputHeader.setSectionResizeMode(3, QHeaderView.ResizeToContents)
		inputBoxLayout.addWidget(self.inputTable,1,0)
		inputBox.setLayout(inputBoxLayout)


		self.main_widget = QWidget(self)
		self.graph_canvas = MyDynamicMplCanvas(self.main_widget, width=5, height=4, dpi=120)
		#Define where the widgets go in the window
		#We start by defining some boxes that we can arrange

		# Results Label Box
		resultsBox = QGroupBox("Results")
		resultsBoxLayout = QVBoxLayout()
		self.resultsBarLabel = QLabel("Optimization Progress",self)
		resultsBoxLayout.addWidget(self.resultsBarLabel)
		self.resultsBar = QProgressBar(self)
		resultsBoxLayout.addWidget(self.resultsBar)
		self.angleLabel = QLabel("Angle = 90",self)
		resultsBoxLayout.addWidget(self.angleLabel)
		self.angleSlider = QSlider(Qt.Horizontal)
		self.angleSlider.setMinimum(0)
		self.angleSlider.setMaximum(360)
		self.angleSlider.setValue(90)
		resultsBoxLayout.addWidget(self.angleSlider)
		self.length1Label = QLabel("Length 1 Not Computed Yet",self)
		resultsBoxLayout.addWidget(self.length1Label)
		self.length2Label = QLabel("Length 2 Not Computed Yet",self)
		resultsBoxLayout.addWidget(self.length2Label)
		self.length3Label = QLabel("Length 3 Not Computed Yet",self)
		resultsBoxLayout.addWidget(self.length3Label)
		self.length4Label = QLabel("Length 4 Not Computed Yet",self)
		resultsBoxLayout.addWidget(self.length4Label)
		self.length5Label = QLabel("Length 5 Not Computed Yet",self)
		resultsBoxLayout.addWidget(self.length5Label)
		self.base1Label = QLabel("Base 1 Location Not Computed Yet",self)
		resultsBoxLayout.addWidget(self.base1Label)
		self.base2Label = QLabel("Base 2 Location Not Computed Yet",self)
		resultsBoxLayout.addWidget(self.base2Label)
		resultsBox.setLayout(resultsBoxLayout)

		#Now we can set all the previously defined boxes into the main window
		grid_layout = QGridLayout()
		grid_layout.addWidget(inputBox,0,0)
		grid_layout.addWidget(resultsBox,1,1)
		grid_layout.addWidget(controlsBox,1,0)
		grid_layout.addWidget(self.graph_canvas,0,1)
		#grid_layout.addWidget(distribution_box,1,1)

		self.centralwidget.setLayout(grid_layout)

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
