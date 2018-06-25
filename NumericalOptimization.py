import sympy as sy
import numpy as np
import random
from sympy import *
import FunctionApproximation as approx
import scipy.optimize
import scipy
import CustomPlots
import math
#from sympy.mpmath import *

def vlen(inputs):
	totalSquared = 0
	for input in inputs:
		totalSquared += input*input
	result = math.sqrt(totalSquared)
	return result

def degToRad(deg):
	rad = deg/360.0
	return rad

def make2dList(rows, cols):
	a=[]
	for row in xrange(rows): a += [[0]*cols]
	return a

def variableSymbols(variables):
	if variables:
		variableSymbols = []
		if isinstance(variables[0],str) == True:
			for variable in variables:
				variableSymbols.append(symbols(variable))
		else:
			variableSymbols = variables

	return variableSymbols

def expressionSymbols(expression):

	if isinstance(expression, str):
		expression = sy.sympify(expression)
	return expression

def gradient(function,inputs,delta=0.0001,normalize=False):
	'''returns a list of partial gradients of the function around the input point'''
	# Inputs: function is a python function that accepts only a list of inputs as arguments
	# Inputs is a list representing the point at which to evaluate the function.
	# Optional: delta is the numerical step size of the gradient approximation
	# Normalize returns the slope of each partial of the gradient divided by the total slope

	slopeValues = []
	for i in range(0,len(inputs)):
		negativeInputs = list(inputs)
		negativeInputs[i] = float(negativeInputs[i]) - float(delta)
		negativePoint = function(negativeInputs)

		positiveInputs = list(inputs)
		positiveInputs[i] = float(positiveInputs[i]) + float(delta)
		positivePoint = function(positiveInputs)

		slope = (positivePoint - negativePoint)/(2*delta)
		slopeValues.append(slope)

	if normalize == True:
		totalSlope = vlen(slopeValues)
		for i in range(0,len(slopeValues)):
			slopeValues[i] = slopeValues[i]/totalSlope
	return slopeValues

def hessian(expression,variables):
	n = len(variables)
	H = make2dList(n, n)

	# Core iteration function
	for i in range(0,n):
		firstPartial = diff(expression,variableSymbol[i])
		for j in range(0,n):
			if i > j:
				H[i][j] = H[j][i]
			else:
				H[i][j] = diff(firstPartial,variableSymbol[j])

	return H

def steepestDescentMinimum(function,startingPoint,epsilon=0.0001,nMax=100,damping=1,echo=False,parabolaFitStepSize = 0.1,constantStepSize = 0.1,**kwargs):
	'''minimizes output of function using steepest descent method'''
	# Inputs: python function which returns a single value and takes an input of a list of values
	# Variables is a list of text inputs for each input variable
	# StartingPoint is a vector of intial points for each input variable
	# Convergence and timeout parameters are optional
	alpha = [-parabolaFitStepSize,0,parabolaFitStepSize]
	i = 0

	# Loop
	shouldContinue = True
	position = startingPoint
	objectiveValue = function(position)
	# print("starting loop...")
	# Print current iteration results
	if echo == True:
		headerString = "Iteration\tPosition\t"
		headerString += "Gradient\t"
		headerString += "F(x)"
		print(headerString)

	while shouldContinue == True:
		i = i+1
		# Get gradient at position
		# print("About to get gradient")
		slopeList = gradient(function,position)
		# print("fitting polynomial...")
		# Get three points in that direction at positions of alpha
		functionValues = []
		for alphaValue in alpha:
			testLocation = []
			for oldPosition, slope in zip(position,slopeList):
				testLocation.append(oldPosition-slope*alphaValue)
			functionValues.append(function(testLocation))
		# Fit parabola to curve
		C = approx.threePointQuadraticApprox(alpha, functionValues)
		# Check parabola is concave up
		# Calculate alpha that gives minimum
		alphaStar = 0.0
		if C[2] < 0:
			print("Fitted parabola is concave down. Minimum alpha value is not bounded.")
			alphaStar = constantStepSize
		elif abs(C[2]) < 0.001:
			print("Shallow gradient, using constant step size")
			alphaStar = constantStepSize
		else:
			(alphaStar,bestY) = minimizeParabola(C)
		# Move to position of calculated alpha
		newPosition = []
		for oldPosition, slope in zip(position,slopeList):
			newPosition.append(oldPosition-slope*damping*alphaStar)
		lastPosition = position
		position = newPosition
		objectiveValueLast = objectiveValue
		objectiveValue = function(position)

		# Print current iteration results
		if echo == True:
			resultsString = "%i        \t" %(i)
			resultsString += "{}\t".format(position)
			resultsString += "{}\t".format(slopeList)
			resultsString += "%2.6f" % (objectiveValue)
			print(resultsString)

		# Check convergence
		deltaObjective = objectiveValueLast - objectiveValue
		#print("Delta Objective = %2.4f" % (float(deltaObjective)))
		if abs(deltaObjective) <= epsilon:
			shouldContinue = False
			print("Local Optimium found")

		#print("About to check iteration maximum")
		if i > nMax:
			print("Function timed out. Returning final result")
			shouldContinue = False

	print("#### - Results - ####")
	print("Position is:")
	print(position)
	print("F = %2.6f" % (objectiveValue))
	return (objectiveValue, position)

def SLP(expression,variables,startingPoint,inequalityConstraints=[],epsilon=0.0001,nMax=100,stepMax = 0.5,saveSequence=False,echo=False):
	# Outputs: optimum position
	if isinstance(expression, str):
		expression = sy.sympify(expression)

	shouldContinue = True
	position = startingPoint
	objectiveValue = evaluateExpression(expression, variables = variables, values = position)
	if echo == True:
		headerString = "Iteration\t"
		for variable in variables:
			headerString += "%s\t" % (variable)
		headerString += "F(x)"
		print(headerString)
	designSequence = [position]

	n = 0
	taylorCoeffs = [0]*len(inequalityConstraints)
	b = [0]*len(inequalityConstraints)

	while shouldContinue == True:
		n = n + 1
		oldPosition = list(position)

		# Linearize objective function and constraints
		(expressionCoeffs,intercept) = approx.taylorLinearize(expression,variables = variables, values = position)
		expressionCoeffs = np.array(expressionCoeffs)
		for i in range(0,len(inequalityConstraints)):
			(taylorCoeffs[i],b[i]) = approx.taylorLinearize(inequalityConstraints[i],variables = variables, values = position)

		taylorArray = np.array(taylorCoeffs)
		# Solve linear problem
		res = scipy.optimize.linprog(expressionCoeffs,A_ub=taylorArray,b_ub=b)
		#print(res)
		# Exctract optimum from result
		newOptimum = res.get("fun", -9999)
		objectiveValueLast = objectiveValue
		objectiveValue = newOptimum
		# Extract optimized design from result
		newPosition  = res.get("x",[-9999]*len(variables))
		newPosition = newPosition.tolist()

		# Check movement delta for each variable
		for i in range(0,len(variables)):
			delta =  newPosition[i] - position[i]
			if abs(delta)> stepMax:
				#print("New position is a large move")
				#print("Former position: " + str(position[i]))
				position[i] = position[i] + delta/abs(delta)*stepMax
				#print("Optimum position: " + str(newPosition[i]))
				#print("Chosen move position: " + str(position[i]))
			else:
				position[i] = newPosition[i]
				#print("New position is a valid move")
		designSequence.append(list(position))

		# Check convergence
		# Print current iteration results
		if echo == True:
			resultsString = "%i        \t" %(n)
			for value in position:
				resultsString += "%2.4f\t" % (value)
			resultsString += "%2.6f" % (objectiveValue)
			print(resultsString)

		# Check convergence
		deltaObjective = objectiveValueLast - objectiveValue
		#print("Last position: " + str(oldPosition))
		#print("Current position: " + str(position))
		variableDeltas = [abs(old - new) for old, new in zip(oldPosition,position)]
		#print(variableDeltas)
		deltaVar = max(variableDeltas)
		#print("Delta Objective = %2.4f" % (float(deltaObjective)))
		if (abs(deltaObjective) <= epsilon and deltaVar <= epsilon):
			shouldContinue = False
			print("Local Optimium found")

		#print("About to check iteration maximum")
		if n > nMax:
			print("Function timed out. Returning final result")
			shouldContinue = False

	print("#### - Results - ####")
	for variable, variableValue in zip(variables,position):
		print(variable + " = %2.6f" % (variableValue))
	print("F = %2.6f" % (objectiveValue))

	if saveSequence == True:
		print(designSequence)
		designSequence = np.array(designSequence)
		x = np.arange(0,3,0.1)
		y = np.arange(0,3,0.1)
		z = make2dList(len(y),len(x))
		constraintValues = []
		for i in range(0,len(inequalityConstraints)):
			constraintValues.append(make2dList(len(y),len(x)))
		for i in range(0,len(x)):
			for j in range(0,len(y)):
				z[j][i] = evaluateExpression(expression,variables = variables,values = [x[i],y[j]])
				for n in range(0,len(inequalityConstraints)):
					constraintValues[n][j][i] = evaluateExpression(inequalityConstraints[n],variables = variables,values = [x[i],y[j]])
		CustomPlots.plotConstrainedContour(x,y,z,"DesignSequence",constraints=constraintValues,lineArray = designSequence)

	return (objectiveValue, position)

def augmentedLagrange(expression,variables,equalityConstraints = [], x0 = [],l0 = 0,epsilon=0.0001,nMax=100,damping=1.0,rp=1.0,echo=False,**kwargs):
	# Inputs: expression is a text string or sympy expression for the objective function
	# Variables is a list of text inputs for each input variable
	# StartingPoint is a vector of intial points for each input variable
	# Convergence and timeout parameters are optional
	alpha = [0,1,2]
	i = 0
	l = [l0]*len(equalityConstraints)
	constraintValues = [0]*len(equalityConstraints)

	if isinstance(expression, str):
		expression = sy.sympify(expression)
	objectiveExpression = expression

	equalityConstraints = [sy.sympify(constraint) for constraint in equalityConstraints]

	#if len(x0) == 0:
		# Calculate starting point from initial lagrange function

	# Loop
	shouldContinue = True
	position = x0
	objectiveValue = evaluateExpression(expression, variables = variables, values = position)
	#print("F = %2.6f" % (objectiveValue))
	# print("About to start loop")
	# Print current iteration results
	if echo == True:
		headerString = "Iteration\t"
		for variable in variables:
			headerString += "%s\t" % (variable)
		for i in range(0,len(l)):
			headerString += "L%i\t" % (i+1)
		for i in range(0,len(equalityConstraints)):
			headerString += "h%i\t" % (i+1)
		headerString += "Phi(x)"
		print(headerString)

	n = 0
	while shouldContinue == True:
		n = n+1

		# Construct expression for augmented lagrange function

		expression = objectiveExpression
		for i in range(0,len(equalityConstraints)):
			expression = expression + rp*equalityConstraints[i]*equalityConstraints[i] + l[i]*equalityConstraints[i]
		#print(expression)

		slopeList = getGradient(expression,variables,position,normalize=False)

		# Get three points in that direction at intervals of 0,0.1,0.2
		functionValues = [objectiveValue]
		for alphaValue in alpha:
			if alphaValue != alpha[0]:
				testLocation = []
				for oldPosition, slope in zip(position,slopeList):
					testLocation.append(oldPosition+slope*alphaValue)
				functionValues.append(evaluateExpression(expression, variables = variables, values = testLocation))
		# Fit parabola to curve
		C = approx.threePointQuadraticApprox(alpha, functionValues)
		# Check parabola is concave up
		# Calculate alpha that gives minimum
		alphaStar = 0.0
		if C[2] < 0:
			print("Fitted parabola is concave down. Minimum alpha value is not bounded.")
			alphaStar = 0.1
		else:
			(alphaStar,bestY) = minimizeParabola(C)
		# Move to position of calculated alpha
		newPosition = []
		for oldPosition, slope in zip(position,slopeList):
			newPosition.append(oldPosition+slope*damping*alphaStar)
		lastPosition = position
		position = newPosition
		objectiveValueLast = objectiveValue
		objectiveValue = evaluateExpression(expression, variables = variables, values = position)

		# Update lagrange multipliers
		for i in range(0,len(equalityConstraints)):
			constraintValues[i] = evaluateExpression(equalityConstraints[i], variables = variables, values = position)
			l[i] = l[i] + 2*rp*constraintValues[i]

		# Print current iteration results
		if echo == True:
			resultsString = "%i        \t" %(n)
			for value in position:
				resultsString += "%2.4f\t" % (value)
			for value in l:
				resultsString += "%2.4f\t" % (value)
			for value in constraintValues:
				resultsString += "%2.4f\t" % (value)
			resultsString += "%2.6f" % (objectiveValue)
			print(resultsString)

		# Check convergence
		deltaObjective = objectiveValueLast - objectiveValue
		#print("Delta Objective = %2.4f" % (float(deltaObjective)))
		if abs(deltaObjective) <= epsilon:
			shouldContinue = False
			print("Local Optimium found")

		#print("About to check iteration maximum")
		if n > nMax:
			print("Function timed out. Returning final result")
			shouldContinue = False

	print("#### - Results - ####")
	for variable, variableValue in zip(variables,position):
		print(variable + " = %2.6f" % (variableValue))
	print("F = %2.6f" % (objectiveValue))
	return (objectiveValue, position)

def quasiNewtonMinimization(expression,variables, startingPoint,epsilon=0.0001,nMax=100,method='bfgs',echo=False):
	x = startingPoint
	i = 0
	shouldContinue = True
	n = len(variables)
	A = np.identity(n)
	alphaTestPoints = [0,0.1,0.2]
	fTestPoints = [0,0,0]
	f = 99999
	delFOld = np.asarray([0,0])
	expression = expressionSymbols(expression)
	variables = variableSymbols(variables)
	xNew = [0,0]
	xOld = [0,0]
	if echo == True:
		headerString = "Iteration \t"
		for variable in variables:
			headerString += str(variable) + "      \t"
		headerString += "F(x)"
		print(headerString)

	while shouldContinue == True:
		(slope, delF) = getGradient(expression,variables,x)
		delF = [-delElement for delElement in delF] # Look downhill rather than uphill
		delF = np.asarray(delF)

		# Calculate values for alpha Star
		j = 0
		for alphaTest in alphaTestPoints:
				xTestPoint = x + alphaTest*(np.dot(A,delF))
				fTestPoints[j] = evaluateExpression(expression,variables = variables,values = xTestPoint)
				j = j + 1

		C = approx.threePointQuadraticApprox(alphaTestPoints,fTestPoints)

		# Check parabola is concave up
		# Calculate alpha that gives minimum
		alphaStar = 0.0
		if C[2] < 0:
			print("Fitted parabola is concave down. Minimum alpha value is not bounded.")
			alphaStar = 0.1
		else:
			(alphaStar,bestY) = minimizeParabola(C)
		xNew = x + alphaStar*(np.dot(A,delF))
		xOld = x
		x = xNew

		# Calculate new A matrix
		if method == 'bfgs':
			p = [xElement - xOldElement for xElement, xOldElement in zip(x, xOld)] # Nx1 vector
			y = delF - delFOld # Nx1 vector
			sigma = np.dot(np.transpose(p),y) # Scalar
			tau = np.dot(np.dot(np.transpose(y),A),y) # Scalar
			D = (sigma+tau)/(sigma*sigma)*np.dot(p,np.transpose(p)) - 1/sigma*(np.dot(np.dot(A,y),np.transpose(p)) + np.dot(p,np.transpose(np.dot(A,y))))
			A = A + D
		elif method == 'DFP':
			print("Implementation of DFP still needed")
		else:
			print("No method selected in quasiNewtonMinimization")

		fNew = evaluateExpression(expression,variables = variables,values = x)
		fDelta = f - fNew
		f = fNew
		delFOld = delF
		i = i + 1

		# Print current iteration results
		if echo == True:
			resultsString = "%i       \t" %(i)
			for variable, value in zip(variables,x):
				resultsString += "%2.4f  \t" % (value)
			resultsString += "%2.6f\t" % (f)
			print(resultsString)

		# Check convergence
		if abs(fDelta) < epsilon:
			shouldContinue = False
			print("Local Optimium found")
		if i > nMax:
			print("Function timed out. Returning final result")
			shouldContinue = False

	print("#### - Results - ####")
	for variable, variableValue in zip(variables,x):
			print(str(variable) + " = %2.6f" % (variableValue))
	print("F = %2.6f" % (f))

	return (f, x)

def NewtonRaphson1DFindZeroUnconstrained(functionString,xStart,tolerance=0.0001,maxIterations=100,echo=False):
	xSymbolic = symbols('x')
	objectiveExpression = sy.sympify(functionString)
	#print(objectiveExpression)
	objectivePrime = diff(objectiveExpression, xSymbolic)
	#print(objectivePrime)

	x = xStart
	shouldContinue = True
	epislon = 1000
	i = 0
	if echo==True:
		print("Iter \t X     \t F      \tF'")

	while shouldContinue == True:
		i = i + 1
		f = evaluateExpression(objectiveExpression, [xSymbolic],[x])
		fPrime = evaluateExpression(objectivePrime, [xSymbolic],[x])
		xNew = x - f/fPrime
		#print(f)
		#print(fPrime)

		epsilon = abs(xNew - x)

		x = xNew

		if epsilon <= tolerance or i >= maxIterations:
			shouldContinue = False
		if echo==True:
			print("%i \t %2.4f \t %2.4f \t %2.4f" % (i, x,f,fPrime))

	return x

def NewtonRaphson1DFindMinUnconstrained(functionString,xStart,tolerance=0.0001,maxIterations=100):
	xSymbolic = symbols('x')
	objectiveExpression = sy.sympify(functionString)
	#print(objectiveExpression)
	objectivePrime = diff(objectiveExpression, xSymbolic)
	#print(objectivePrime)
	objectiveDoublePrime = diff(objectivePrime, xSymbolic)
	#print(objectiveDoublePrime)
	x = xStart
	shouldContinue = True
	epislon = 1000
	i = 0

	while shouldContinue == True:
		i = i + 1
		f = evaluateExpression(expression = objectiveExpression, x = x)
		fPrime = evaluateExpression(expression = objectivePrime, x = x)
		fDoublePrime = evaluateExpression(expression = objectiveDoublePrime, x = x)
		xNew = x - fPrime/fDoublePrime
		#print(f)
		#print(fPrime)

		epsilon = abs(xNew - x)

		x = xNew

		if epsilon <= tolerance or i >= maxIterations:
			shouldContinue = False
		print("Iteration = %i, X = %2.4f, F = %2.4f, F' = %2.4f" % (i, x,f,fPrime))

	error = evaluateExpression(expression = objectiveExpression, x = x)

	return x

def evaluateExteriorPenalty(function, position,
	inequalityConstraints=[],
	equalityConstraints=[], rp=1):
	"""returns a float at the location selected with constraint penalties"""

	objectiveValue = function(position)

	penalty_total = 0
	for constraint in inequalityConstraints:
		penalty = constraint(position)
		penalty = max(0,penalty)**2
		penalty_total +=  penalty

	for constraint in equalityConstraints:
		penalty = constraint(position)**2
		penalty_total +=  penalty

	totalValue = objectiveValue + rp * penalty_total
	result = totalValue

	return result

def evaluateLinearExtendedPenalty(
	function, position,
	inequalityConstraints=[],
	equalityConstraints=[],
	rp=1.0,
	epsilon = -9999):
	"""returns a float at the location selected with constraint penalties"""

	if epsilon == -9999:
		epsilon = -0.2*np.sqrt(1/rp)

	rpPrime = 1/rp
	objectiveValue = function(position)

	inconstraintValue = 0
	for constraint in inequalityConstraints:
		newConstraintValue = constraint(position)
		if newConstraintValue > epsilon:
			inconstraintValue += - (2*epsilon - newConstraintValue)/epsilon**2
		else:
			inconstraintValue = inconstraintValue - 1/newConstraintValue

	constraintValue = 0
	for constraint in equalityConstraints:
		newConstraintValue = constraint(position)**2
		constraintValue = constraintValue + newConstraintValue

	totalValue = objectiveValue + inconstraintValue/rp + constraintValue*rp
	result = totalValue

	return result

def evaluateInteriorInverseBarrierPenalty(
	function, position,
	inequalityConstraints=[],
	equalityConstraints=[],
	rp=1.0):
	"""returns a float at the location selected with constraint penalties"""

	objectiveValue = function(position)

	ineq_constraint_penalty = 0
	for constraint in inequalityConstraints:
		constraint_value = constraint(position)
		if ineq_constraint_penalty <= 0:
			ineq_constraint_penalty += - 1/constraint_value
		else:
			ineq_constraint_penalty += 100*rp * constraint_value

	eq_constraint_penalty = 0
	for constraint in equalityConstraints:
		constraint_value = constraint(position)**2
		eq_constraint_penalty += constraint_value

	result = objectiveValue + ineq_constraint_penalty/rp + rp * eq_constraint_penalty

	return result

def constrainedMinimum(function,startingPoint,
	inequalityConstraints=[],
	equalityConstraints=[],
	rp=1,
	method='ExteriorPenalty',
	echo=False,
	damping=0.1,
	epsilon=0.0001,
	nMax=100,
	parabolaFitStepSize = 0.1,
	constantStepSize = 0.1,
	printResults=True,
	**kwargs):
	'''minimizes the given function for n variables subject to boundary constraints'''
	# Input: function whose only argument is a list of values and which returns a single value
	# StartingPoint is a list of values corrosponding to the number of variables
	# Constraints are functions which take a list of values and return a single value
	# Inequality constraints return less than 0 when valid, equality equal 0
	# Method options: 'ExteriorPenalty', 'InteriorPenalty', 'InteriorInverseBarrier','InverseLog', 'InteriorLinearExtended', 'QuadraticExtended'
	if method == 'ExteriorPenalty':
		penalizedFunction = lambda position: evaluateExteriorPenalty(function,
			inequalityConstraints=inequalityConstraints,
			equalityConstraints=equalityConstraints,
			position = position,
			rp = rp)
	elif method == 'InteriorLinearExtended':
		penalizedFunction = lambda position: evaluateLinearExtendedPenalty(function,
			inequalityConstraints=inequalityConstraints,
			equalityConstraints=equalityConstraints,
			position = position,
			rp = rp,
			epsilon = -9999)
	elif method == 'InteriorInverseBarrier':
		penalizedFunction = lambda position: evaluateInteriorInverseBarrier(function,
			inequalityConstraints=inequalityConstraints,
			equalityConstraints=equalityConstraints,
			position = position,
			rp = rp)
	else:
		print('The method ' + method + ' is not implemented yet.')
		return
	(optimum, position) = steepestDescentMinimum(penalizedFunction, startingPoint,
	epsilon=epsilon,
	nMax=nMax,
	damping=damping,
	echo=echo,
	parabolaFitStepSize = parabolaFitStepSize,
	constantStepSize = constantStepSize,**kwargs)

	return (optimum, position)

def minimizeCubic(c):
	# Inputs: Coefficients for polynomial equation according to the form C0 + C1*x + C2*x^2 + C3*x^3
	# Outputs: Values of x and y where y is minimized
	a = 3*c[3]
	b = 2*c[2]
	d = c[1]
	insideSqareroot = np.float64(b*b-4*a*d)
	if insideSqareroot < 0:
		print("Minimize Cubic function encountered imaginary square root. Aborting.")
		return
	x1 = (-b+np.sqrt(insideSqareroot))/(2*a)
	x2 = (-b-np.sqrt(insideSqareroot))/(2*a)

	x = 0
	y = 0

	y1 = approx.getValueOfPoly(c,x1)
	y2 = approx.getValueOfPoly(c,x2)
	if y1 < y2:
		x = x1
		y = y1
	elif y1 > y2:
		x = x2
		y = y1
	else:
		x = x1
		y = y1
		print("More than one solution in Minimize Cubic")

	return (x,y)

def minimizeParabola(c):
	# Inputs: Coefficients for polynomial equation according to the form C0 + C1*x + C2*x^2...
	# Outputs: Values of x and y where y is minimized
	minX = -c[1]/(2*c[2])

	minY = approx.getValueOfPoly(c,minX)
	return (minX,minY)

def convertToPenaltyFunction(coreFunction,constraints,R=1):
	constraintsToSum = []
	newObjective = coreFunction + " - %2.4f*(" % (R)
	for i in range(0,len(constraints)):
		constraint = constraints[i]
		if i == 0:
			newObjective = newObjective + "1/(" + constraint + ")"
		else:
			newObjective = newObjective + " + 1/(" + constraint + ")"


	newObjective = newObjective + ")"
	return newObjective

def goldenSectionSearch(expression,xlow,xu,epsilon = 0.001,n=100,echo=False):
	tau = 0.381966

	if isinstance(expression, str):
		expression = sy.sympify(expression)

	fu = evaluateExpression(expression,variables = ['x'], values = [xu])
	flow = evaluateExpression(expression,variables = ['x'], values = [xlow])

	x1 = (1-tau)*xlow + tau*xu
	f1 = evaluateExpression(expression,variables = ['x'], values = [x1])
	x2 = tau*xlow + (1-tau)*xu
	f2 = evaluateExpression(expression,variables = ['x'], values = [x2])

	k = 3
	shouldContinue = True

	while shouldContinue == True:

		if f1 > f2:
			xlow = x1
			flow = f1
			x1 = x2
			f1 = f2
			x2 = tau*xlow + (1-tau)*xu
			f2 = evaluateExpression(expression,variables = ['x'], values = [x2])
		else:
			xu = x2
			fu = f2
			x2 = x1
			f2 = f1
			x1 = (1-tau)*xlow + tau*xu
			f1 = evaluateExpression(expression,variables = ['x'], values = [x1])

		k = k + 1
		if echo == True:
			print("i = %i \t xLow = %2.4f \t F(xLow) = %2.4f \t xHigh = %2.4f \t F(xHigh) = %2.4f" % (k,xlow,flow,xu,fu))

		if k > n:
			shouldContinue = False
		if xu - xlow < epsilon:
			shouldContinue = False
	fs = [f1,f2,flow,fu]
	xs = [x1,x2,xlow,xu]
	fMin = min(fs)
	xMin = xs[fs.index(fMin)]

	return(fMin,xMin)

def randomSeaoptimimrch2D(objectiveFunction,xStart,yStart,constraints,tolerance=0.0001, maxIterations=100):
	objectiveExpression = sy.sympify(objectiveFunction)
	constraintExpressions = []
	for constraint in constraints:
		constraintExpressions.append(sy.sympify(constraint))

	# Variable initializations
	xBest = 999999999.9
	yBest = 999999999.9
	objectiveBest = 9999999.9999
	x = xStart
	y = yStart
	epsilon = 100
	shouldContinue = True
	i = 0
	sinceLastPrint = 0
	printInterval = 100

	# Iteration loop
	while shouldContinue == True:
		i = i+1
		sinceLastPrint = sinceLastPrint + 1
		if sinceLastPrint >= printInterval:
			print("Running Iteration %i" %(i))
			sinceLastPrint = 0

		xNew = x + random.uniform(-0.1, 0.1)
		yNew = y + random.uniform(-0.1,0.1)
		objectiveNew = evaluateExpression(objectiveExpression,x=xNew,y=yNew)
		validPoint = True
		if objectiveNew < objectiveBest:
			for g in constraintExpressions:
				if evaluateExpression(expression=g,x=x,y=y) > 0:
					validPoint = False

			if validPoint == True:
				# Move and store new location
				xLast = x
				yLast = y
				objectiveLast = objectiveBest

				x = xNew
				y = yNew
				objectiveBest = objectiveNew

				# Check convergence
				epsilon =  objectiveLast - objectiveBest
				print("Best solution so far: %2.4f" % (objectiveBest))

		if epsilon <= tolerance or i >= maxIterations:
			shouldContinue = False

	return (x, y, objectiveBest)

def bruteForceMinimum2D(objectiveFunction,xSearchRange,ySearchRange,constraints,resolution):

	objectiveExpression = sy.sympify(objectiveFunction)
	constraintExpressions = []
	for constraint in constraints:
		constraintExpressions.append(sy.sympify(constraint))

	xArray = np.arange(xSearchRange[0], xSearchRange[1], resolution).tolist()
	yArray = np.arange(ySearchRange[0], ySearchRange[1], resolution).tolist()

	xBest = 999999999.9
	yBest = 999999999.9
	zBest = 9999999.9999

	# Iteration loop
	for x in xArray:
		print("x position: %2.4f" % (x))
		for y in yArray:
			z = evaluateExpression(expression=objectiveExpression, x=x,y=y)
			validPoint = True
			if z < zBest:
				for g in constraintExpressions:
					if evaluateExpression(expression=g,x=x,y=y) > 0:
						validPoint = False

				if validPoint == True:
					xBest = x
					yBest = y
					zBest = z

	return (xBest, yBest, zBest)
