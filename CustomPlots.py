import sympy as sy
import numpy as np
import random
from sympy import *
import FunctionApproximation as approx
import matplotlib.pyplot as plt

def plotLines(x,y,saveName,xLabel='x',yLabel='y',shouldShow=False):

	saveName = saveName + '.png'

	fig = plt.figure()
	fig.set_size_inches(6,4)
	for ySet in y:
		line, = plt.plot(x,ySet)
	plt.xlabel(xLabel)
	plt.ylabel(yLabel)
	plt.savefig(saveName, dpi=120)
	if shouldShow == True:
		plt.show()
	plt.close()

def plotConstrained1D(x,y,constraints,saveName,xLabel='x',yLabel='y',shouldShow = False):
	saveName = saveName + '.png'

	fig, ax = plt.subplots(1)	
	fig.set_size_inches(6,4)
	line, = plt.plot(x,y)
	gMin = 999999
	gMax = -999999
	for g in constraints:
		ax.fill_between(x, -100, g,facecolor='red',alpha=0.5,interpolate=True)
		thisGMin = min(g)
		thisGMax = max(g)
		if thisGMin < gMin:
			gMin = thisGMin
		if thisGMax > gMax:
			gMax = thisGMax


	yMax = max(max(y),gMax)
	yMin = min(min(y),gMin)
	plt.ylim(yMin,yMax)
	plt.xlim(min(x),max(x))
	plt.xlabel(xLabel)
	plt.ylabel(yLabel)
	ax.axhline(y=0, color='k')
	ax.axvline(x=0, color='k')
	plt.savefig(saveName, dpi=120)
	if shouldShow == True:
		plt.show()
	plt.close()

def plotConstrainedContour(x,y,z,saveName,constraints=[],lineArray = [], xRange=[],yRange=[],shouldShow=False,levels=[]):
	# lineArray is a 2D numpy array 
	saveName = saveName + '.png'

	if xRange == []:
		xRange = [min(x),max(x)]

	if yRange == []:
		yRange = [min(y),max(y)]

	XPlot, YPlot = np.meshgrid(x, y)
	plt.figure()
	if len(levels) > 0:
		CS = plt.contour(XPlot, YPlot, z,levels = levels)
	else:
		CS = plt.contour(XPlot, YPlot, z)

	for i in range(0,lineArray[:,0]-1):
		plt.plot([i, 0], [i+1, 1], '-k')

	for g in constraints:
		plt.contourf(XPlot,YPlot,g,[0,10000],colors=('r'),alpha=0.5)
	plt.xlim(xRange[0],xRange[1])
	plt.ylim(yRange[0],yRange[1])
	plt.clabel(CS, inline=1, fontsize=10)
	#plt.title('Simplest default with labels')
	plt.savefig(saveName, dpi=120)
	if shouldShow == True:
		plt.show()
	plt.close()
