import sympy as sy
import numpy as np
import random
from sympy import *

def threePointCubicApprox(x,y,xSlopePoint,yPrime):
	C3 = (y[2] - y[0])/((x[2] - x[1])*(x[2] - x[0])**2) - (y[1] - y[0])/((x[2] - x[1])*(x[1] - x[0])**2) + yPrime/((x[1]-x[0])*(x[2] - x[0]))
	C2 = (((y[1] - y[0])/(x[1] - x[0])) - yPrime)/(x[1] - x[0]) - C3*(2*x[0] + x[1])
	C1 = yPrime - 2*C2*x[0] - 3*C3*x[0]**2
	C0 = y[0] - C1*x[0] - C2*x[0]**2 - C3*x[0]**3
	
	return [C0,C1,C2,C3]

def threePointQuadraticApprox(x, y):
	#Inputs: Vector of x values and y values. These vectors must be equal length
	#Outputs: Coefficients for polynomial equation according to the form C0 + C1*x + C2*x^2...
	C2 = (((y[2]-y[0])/(x[2]-x[0])) - ((y[1]-y[0])/(x[1]-x[0])))/(x[2]-x[1])
	C1 = (y[1] - y[0])/(x[1]-x[0]) - C2*(x[0]+x[1])
	C0 = y[0] - C1*x[0] - C2*x[0]**2

	return [C0,C1,C2]

def twoPointLinearApprox(x, y):
	#Inputs: Vector of x values and y values. These vectors must be equal length
	#Outputs: Coefficients for polynomial equation according to the form C0 + C1*x + C2*x^2...
	C1 = (y[1] - y[0])/(x[1]-x[0])
	C0 = y[0] - C1*x[0]

	return [C0,C1]

def getValueOfPoly(c,x):
	#Inputs: Coefficients for polynomial equation according to the form C0 + C1*x + C2*x^2...
	#Inputs: x - value to get value at
	constantQuantity = len(c)

	if constantQuantity == 1:
		# Flat line
		y = c[0]
	elif constantQuantity == 2:
		# Linear
		y = c[0] + c[1] * x
	elif constantQuantity == 3:
		# Quadratic
		y = c[0] + c[1]*x + c[2]*x**2
	elif constantQuantity == 4:
		# Cubic
		y = c[0] + c[1]*x + c[2]*x**2 + c[3]*x**3
	else:
		print("Polynomial could not be calculated. Check getValueOfPoly function.")
		y = 99999999

	return y