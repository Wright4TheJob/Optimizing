from PyQt5.QtCore import QThread
import PyQt5.QtCore as QtCore

class OptimizeThread(QThread):

	iterationDone = QtCore.Signal(object)
	designOptimized = QtCore.Signal(object)
	def __init__(self, design):
		QThread.__init__(self)
		self.design = list(design)
		#self.variables = list(variables)

	def __del__(self):
		self.wait()

	def runOptimization(self, design):
		result = []
		result.append(sum(design))
		self.iterationDone.emit(result)
		return design

	def run(self):
		self.design = self.runOptimization(self.design)
		self.designOptimized.emit(self.design)
		
		self.sleep(2)
		