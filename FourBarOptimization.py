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
                             QGridLayout,QGroupBox, QMainWindow)
from PyQt5.QtCore import Qt, QTimer, QCoreApplication

from scipy import stats

import statistics
import numpy as np
import os

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
             
    def plotFourBar(self,targets,mechanismData,data_label="Mechanism",title="Mechanism"):
        # Plot mechanism bars

        xTargets = targets[1]
        yTargets = targets[2]

        self.axes.cla() #Clear axes
        
        # Plot target positons
        self.axes.plot(xTargets, yTargets, 'rx')
        # Plot actual path

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
        self.damping = 1.0
        self.thetaTargets = [0,90,180,270]
        self.xTargets = [1,4,6,3]
        self.yTargets = [1,0,1,2]
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
        self.dampingSlider.setValue(1000)
        self.dampingSlider.valueChanged.connect(self.dampingChanged)
        self.dampingLabel = QLabel("Damping = %1.2f"%(self.damping),self)

        # Let's make a simple label widget to keep track of a count
        self.slider_label = QLabel("LED Blink Spacing")
        self.iterationLabel = QLabel("Iterations Not Started",self)
        self.length1Label = QLabel("Length 1 Not Computed Yet",self)
        self.length2Label = QLabel("Length 2 Not Computed Yet",self)
        self.length3Label = QLabel("Length 3 Not Computed Yet",self)
        self.length4Label = QLabel("Length 4 Not Computed Yet",self)
        self.length5Label = QLabel("Length 5 Not Computed Yet",self)
        
        #Set up a Table to display data
        self.data_table = QTableWidget()
        self.data_table.itemSelectionChanged.connect(self.compute_stats)
        
        self.main_widget = QWidget(self)
        self.graph_canvas = MyDynamicMplCanvas(self.main_widget, width=5, height=4, dpi=100)
        
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
        resultsBox.setLayout(resultsBoxLayout)

        # Controls Box
        controlsBox = QGroupBox("Controls")
        controlsBoxLayout = QVBoxLayout()
        controlsBoxLayout.addWidget(self.runButton)
        controlsBoxLayout.addWidget(self.dampingLabel)
        controlsBoxLayout.addWidget(self.dampingSlider)
        controlsBox.setLayout(controlsBoxLayout)
        
        #Create some distribution options
        #Start with creating a check button.
        self.normal_checkbox = QCheckBox('Normal Distribution',self)
        # We want to run the plotting routine for the distribution, but 
        # we need to know the statistical values, so we'll calculate statistics
        # first.
        self.normal_checkbox.stateChanged.connect(self.compute_stats)
        
        
        distribution_box = QGroupBox("Distribution Functions")
        distribution_box_layout= QVBoxLayout()
        distribution_box_layout.addWidget(self.normal_checkbox)
        distribution_box.setLayout(distribution_box_layout)
        #Now we can set all the previously defined boxes into the main window
        grid_layout = QGridLayout()
        grid_layout.addWidget(table_box,0,0) 
        grid_layout.addWidget(resultsBox,1,0)
        grid_layout.addWidget(controlsBox,1,1)
        grid_layout.addWidget(self.graph_canvas,0,1) 
        #grid_layout.addWidget(distribution_box,1,1)
        
        self.setLayout(grid_layout)
        
        self.setWindowTitle('Four Bar Linkage Optimization')
        self.activateWindow()
        self.raise_()
        self.show()
    
    def load_data(self):
        #Write this function
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)
        if fileName:
            f= open(fileName,"r")
            if f.mode == 'r':
                contents =f.read()
                # Do stuff with contents
    
    def save_data(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()","","All Files (*);;Text Files (*.txt)", options=options)
        # Check to make sure ends in .txt
        if fileName:
            f= open(fileName,"w+")
            f.write("%i" % (self.counter_value))
            f.close() 

    def runFourBarOptimization(self):
        print('Analysis should go here')

        """
        #print("Mean = {0:5f} ".format(mean_value))
        self.mean_label.setText("Mean = {0:3f} ".format(mean_value))
        #print("Standard Error = {0:5f} ".format(stdError_value))
        self.length1Label.setText("Length 1 = {0:3f} ".format(stdError_value))
        #print("Median = {0:5f} ".format(median_value))
        self.median_label.setText("Median = {0:3f} ".format(median_value))
        #print("Mode = {0:5f} ".format(mode_value))
        self.mode_label.setText("Mode = {0:3f} ".format(mode_value))
        #print("Standard Deviation = {0:5f} ".format(stdDev_value))
        self.stdDev_label.setText("Standard Deviation = {0:3f} ".format(stdDev_value))
        #print("Stanard Variance = {0:5f} ".format(var_value))
        self.var_label.setText("Standard Variance = {0:3f} ".format(var_value))
        #print("Kurtosis = {0:5f} ".format(kurt_value))
        self.kurt_label.setText("Kurtosis = {0:3f} ".format(kurt_value))
        #print("Skew = {0:5f} ".format(skew_value))
        self.skew_label.setText("Skewness = {0:3f} ".format(skew_value))
        #print("Max = {0:5f} ".format(max_value))
        self.max_label.setText("Max = {0:3f} ".format(max_value))
        #print("Min = {0:5f} ".format(min_value))
        self.min_label.setText("Min = {0:3f} ".format(min_value))
        #print("Range  = {0:5f} ".format(range_value))
        self.range_label.setText("Range = {0:3f} ".format(range_value))
        #print("Sum = {0:5f} ".format(sum_value))
        self.sum_label.setText("Sum = {0:3f} ".format(sum_value))
        #print("Count = {0:5f} ".format(count_value))
        self.count_label.setText("Number of Entries = {0:3f} ".format(count_value))
        """
        mechanismData = 0
        targets = [self.thetaTargets,self.xTargets,self.yTargets]
        self.graph_canvas.plotFourBar(targets, mechanismData)
        #add more distributions here

    def dampingChanged(self,value):
        self.damping = float(value)/1000
        self.dampingLabel.setText("Damping = %1.2f"%(self.damping))
        
    def compute_stats(self):
        
        #setup array
        item_list=[]
        items = self.data_table.selectedItems()
        for item in items:
            try:
                item_list.append(float(item.text()))
            except:
                pass
        
        if len(item_list) > 1: #Check to see if there are 2 or more samples
            data_array = np.asarray(item_list)
            mean_value = np.mean(data_array)            
            stdDev_value = np.std(data_array)
            stdError_value = stats.sem(data_array)
            median_value = np.median(data_array)
            mode_value = statistics.mode(data_array)
            var_value = np.var(data_array)
            kurt_value = stats.kurtosis(data_array)
            skew_value = stats.skew(data_array)
            max_value = max(data_array)
            min_value = min(data_array) 
            range_value = max_value-min_value
            sum_value = sum(data_array)
            count_value = len(data_array)
        
            
            
            #print("Mean = {0:5f} ".format(mean_value))
            self.mean_label.setText("Mean = {0:3f} ".format(mean_value))
            #print("Standard Error = {0:5f} ".format(stdError_value))
            self.stdError_label.setText("Standard Error = {0:3f} ".format(stdError_value))
            #print("Median = {0:5f} ".format(median_value))
            self.median_label.setText("Median = {0:3f} ".format(median_value))
            #print("Mode = {0:5f} ".format(mode_value))
            self.mode_label.setText("Mode = {0:3f} ".format(mode_value))
            #print("Standard Deviation = {0:5f} ".format(stdDev_value))
            self.stdDev_label.setText("Standard Deviation = {0:3f} ".format(stdDev_value))
            #print("Stanard Variance = {0:5f} ".format(var_value))
            self.var_label.setText("Standard Variance = {0:3f} ".format(var_value))
            #print("Kurtosis = {0:5f} ".format(kurt_value))
            self.kurt_label.setText("Kurtosis = {0:3f} ".format(kurt_value))
            #print("Skew = {0:5f} ".format(skew_value))
            self.skew_label.setText("Skewness = {0:3f} ".format(skew_value))
            #print("Max = {0:5f} ".format(max_value))
            self.max_label.setText("Max = {0:3f} ".format(max_value))
            #print("Min = {0:5f} ".format(min_value))
            self.min_label.setText("Min = {0:3f} ".format(min_value))
            #print("Range  = {0:5f} ".format(range_value))
            self.range_label.setText("Range = {0:3f} ".format(range_value))
            #print("Sum = {0:5f} ".format(sum_value))
            self.sum_label.setText("Sum = {0:3f} ".format(sum_value))
            #print("Count = {0:5f} ".format(count_value))
            self.count_label.setText("Number of Entries = {0:3f} ".format(count_value))
            
            
            self.graph_canvas.plotFourBar(mechanismData)
            #add more distributions here
        
'''       
  1. Add the ability to plot a normalized Histogram of the selected data in the table.*
  2. Add a menu option to open a CSV data file.*
  3. Add a checkbox for at least 5 distribution functions to plot over the top of the Histogram. 
    a. Include a legend and appropriate labels on your graph.
    b. Include axes labels. (Challenge: make the labels editable in your program).
  4. Use a grid style layout for your GUI*
  5. Save the plot to a PDF when you double click on it.*
  6. Try to find one of the most obscure distributions as one of your 5. Please try to be different than everyone else. 
  7. Print and turn in a screenshot of your GUI on one page. Be sure your name in in the window title.
  8. Print and turn in the PDF file of the properly labeled Histogram with 2 distributions shown.
'''

if __name__ == '__main__':
    #Start the program this way according to https://stackoverflow.com/questions/40094086/python-kernel-dies-for-second-run-of-pyqt5-gui
    app = QCoreApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    execute = FourBarOptimizer()
    sys.exit(app.exec_())
