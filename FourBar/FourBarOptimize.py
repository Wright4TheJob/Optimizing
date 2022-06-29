from PyQt5.QtCore import QThread
import PyQt5.QtCore as QtCore
from sympy import *
import sympy as sy
import numpy as np
import math
from numpy.linalg import inv,pinv

class OptimizeThread(QThread):

    iterationDone = QtCore.pyqtSignal(object)
    designOptimized = QtCore.pyqtSignal(object)
    def __init__(self, design, targets, controls):
        QThread.__init__(self)
        self.design = design
        self.optimizedDesign = self.design
        self.targets = targets

        self.objectiveSequence = []
        # Controls = [damping, max Iterations]
        self.damping = controls[0]
        self.maxIterations = controls[1]

        self.mechanismStartAngle = 0

    def __del__(self):
        self.wait()

    def runOptimization(self):
        # Constraint functions
        # Dyad validity analysis
        # TODO: Add penalty function for invalid dyad

        #rp = 1
        #rpMax = 1000000

        #n = 34.
        nMax = self.maxIterations
        #while rp < rpMax:
            #n = n+1.0
        (error,optimumParameters) = steepestDescentMinimum(
            self.evaluateObjectiveFunction,
            design2list(self.optimizedDesign),
            nMax = self.maxIterations,
            damping=self.damping)

        # Return to main thread
        #progress = float(n/nMax)
        self.optimizedDesign = list2design(optimumParameters)
        self.iterationDone.emit(self.optimizedDesign)
        #rp = rp*1.5
        self.objectiveSequence.append(error)
        print(rounded_list(optimumParameters,3))
        return self.optimizedDesign

    def run(self):
        self.optimizedDesign = self.runOptimization()
        self.designOptimized.emit(self.optimizedDesign)

        self.sleep(2)

    def evaluateObjectiveFunction(self, design,rp=1):
        targets = self.targets
        objective = 0.0
        for target in targets:
            # Calculate position at crank angle target[0]
            try:
                position = calculateFourBarPoint(design,target[0])
                objective = objective + vLen(target[1:],position)**2
            except ValueError:
                objective += 2

        #TODO: Add lengths of bars as minor addition to objective function
        lengths = design[4:]
        #print(lengths)
        objective = objective + sum(lengths)/1000
        return objective

def design2list(design):
    return design[0][0] + design[0][1] + design[1]

def list2design(param_list):
    bases = [param_list[0:2], param_list[2:4]]
    lengths = param_list[4:]
    return [bases, lengths]

def calculateFourBarPose(design, alpha, link_positive=True, dyad_positive=True):
    # lengths is in format [base, z1,z2,z3,z4,z5]
    # bases are in format [[base1X,base1Y],[base2X,base2Y]]
    # alpha is a single float for crank angle from horizontal

    # for a given crank angle, there are four potential poses to consider:
    # two possible orientations for the dyad and two possible orientations
    # for the trangle composted of the two passive links, the tip of the crank,
    # and the second base
    # All four possibilities are valid when the distance between the tip of the crank
    # and the second base is less than the sum of the two passive links
    # (the base of the dyad and the crank or rocker)

    # Calculate end point of crank
    # Calculate two possible intersections of dyad base and link (could be 2, 1, or 0)
    # Calculate two possible diad orientations
    # Total of 0, 2, or 4 solutions

    # Constraints:
    # all lengths > 0
    # Sum of two Dyad legs >= dyad base

    #         C
    #     3  / \ 4
    #      /    \
    #     1--------2
    # 0  /    1    \
    #   /          \ 2
    #  B1---- 5 --- B2
    bases = design[0]
    length_without_base = design[1]

    lengths = length_without_base + [vLen(bases[0],bases[1])]

    # Crank Endpoint
    point1 = vector_add_length_angle(bases[0], lengths[0], alpha)

    try:
        point2options = link_intersections(point1, bases[1], lengths[1], lengths[2])
        if link_positive:
            point2 = point2options[0]
        else:
            point2 = point2options[1]
    except ValueError:
        raise ValueError("Mechanism base is not solvable in this configuration")

    try:
        pointCoptions = link_intersections(point1, point2, lengths[4], lengths[5])
        if dyad_positive:
            pointC = pointCoptions[0]
        else:
            pointC = pointCoptions[1]
    except ValueError:
        raise ValueError("Mechanism dyad is not solvable in this configuration")

    mechanismPoints = [bases[0],bases[1],point1,point2,pointC]
    return mechanismPoints

def vector_add_length_angle(start, length, angle):
    # Start is a 2-element list of x,y coordinates
    # Length is a value
    # Angle is a value in radians
    end = [0,0]
    end[0] = start[0] + length*np.cos(angle)
    end[1] = start[1] + length*np.sin(angle)
    return end

def dummy_pose(design,alpha):
    bases = design[0]
    length_without_base = design[1]

    lengths = [vLen(bases[0],bases[1])] + length_without_base

    # Crank Endpoint
    point1 = [0,0]
    point1[0] = bases[0][0] + lengths[1]*np.cos(alpha)
    point1[1] = bases[0][1] + lengths[1]*np.sin(alpha)
    mechanismPoints = [bases[0],bases[1],point1,bases[0],bases[0]]
    return mechanismPoints

def calculateFourBarPoint(design,alpha):
    # lengths is in format [base, z1,z2,z3,z4,z5]
    # bases are in format [[base1X,base1Y],[base2X,base2Y]]
    # alpha is a single float for crank angle from horizontal

    # Assumes Nested list design structure
    if len(design) != 2:
        optimizing_mode = True
        #print("Making design a list in calculateFourBarPoint")
        design = list2design(design)
    else:
        design = list(design)

    try:
        pose = calculateFourBarPose(design, alpha)
    except ValueError:
        raise ValueError("No end point for invalid mechanism")
    #print(pointC)
    return pose[-1]

def link_intersections(base1, base2, l1, l2):
    if vLen(base1,base2) > l1 + l2:
        raise ValueError("Links do not intersect")
    elif vLen(base1,base2) == l1 + l2:
        l1_fraction = l1/(l1 + l2)
        x = base1[0] + dx(base1, base2) * l1_fraction
        y = base1[1] + dy(base1, base2) * l1_fraction
        return [[x, y]]
    else:
        alpha = angle_from_triangle_lengths(l2, l1,  vLen(base1,base2))
        base_angle = absolute_vector_angle(base1, base2)
        x1 = base1[0] + l1*cos(base_angle + alpha)
        y1 = base1[1] + l1*sin(base_angle + alpha)
        x2 = base1[0] + l1*cos(base_angle - alpha)
        y2 = base1[1] + l1*sin(base_angle - alpha)
        return [[x1, y1],[x2, y2]]

def absolute_vector_angle(a, b):
    # returns the angle of the vector from a to b relative to positive horizontal
    dx = np.float64(b[0] - a[0])
    dy = np.float64(b[1] - a[1])
    if dx == 0 and dy == 0:
        raise ValueError("Non-zero vector needed to calculate angle")

    if dx == 0:
        # line is vertical
        if dy > 0:
            return np.pi/2
        else:
            return 3*np.pi/2
    elif dy == 0:
        # line is horizontal
        if dx > 0:
            return 0
        else:
            return np.pi

    atan_value = np.arctan(dy/dx)
    if in_quadrent_1(dx, dy):
        return atan_value
    elif in_quadrent_2(dx, dy):
        return np.pi + atan_value
    elif in_quadrent_3(dx, dy):
        return np.pi + atan_value
    elif in_quadrent_4(dx, dy):
        return 2*np.pi + atan_value
    else:
        raise Error("Unexpected geometry condition encountered")


def in_quadrent_1(dx,dy):
    return (dx > 0 and dy > 0)

def in_quadrent_2(dx,dy):
    return (dx < 0 and dy > 0)

def in_quadrent_3(dx,dy):
    return (dx < 0 and dy < 0)

def in_quadrent_4(dx,dy):
    return (dx > 0 and dy < 0)

def dx(v1, v2):
    # distance in x between two 2D points
    return v2[0] - v1[0]

def dy(v1, v2):
    # distance in x between two 2D points
    return v2[1] - v1[1]

def vLen(point1,point2):
    # takes in two points in the format [x,y] and returns the float of vector length
    dx = point1[0] - point2[0]
    dy = point1[1] - point2[1]

    length = math.sqrt(dx*dx + dy*dy)
    return length

def threePointQuadraticApprox(x,y):
    C2 = (((y[2]-y[0])/(x[2]-x[0])) - ((y[1]-y[0])/(x[1]-x[0])))/(x[2]-x[1])
    C1 = (y[1] - y[0])/(x[1]-x[0]) - C2*(x[0]+x[1])
    C0 = y[0] - C1*x[0] - C2*x[0]**2

    return [C0,C1,C2]

def minimizeParabola(c):
    # Inputs: Coefficients for polynomial equation according to the form C0 + C1*x + C2*x^2...
    # Outputs: Values of x and y where y is minimized
    minX = -c[1]/(2*c[2])

    minY = getValueOfPoly(c,minX)
    return (minX,minY)

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

def steepestDescentMinimum(function,
    startingPoint,
    epsilon=0.0001,
    nMax=100,
    damping=0.1,
    echo=False,
    parabolaFitStepSize = 0.01,
    constantStepSize = 0.1,**kwargs):
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
        C = threePointQuadraticApprox(alpha, functionValues)
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

def angle_from_triangle_lengths(a,b,c):
    # Sides defined counter-clockwise from starting
    # Return angle is opposite side a
    acos_input = ((b*b) + (c*c) - (a*a)) / (2*b*c)
    alpha = math.acos(acos_input)
    return alpha

def rounded_list(xs, places):
    return [round(x,places) for x in xs]
