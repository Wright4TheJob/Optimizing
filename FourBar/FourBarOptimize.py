from PyQt5.QtCore import QThread
import PyQt5.QtCore as QtCore
from sympy import *
import sympy as sy
import numpy as np
import math
from numpy.linalg import inv,pinv

class OptimizeThread(QThread):

    iterationDone = QtCore.Signal(object)
    designOptimized = QtCore.Signal(object)
    def __init__(self, design, targets, controls):
        QThread.__init__(self)
        self.design = list(design)
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
        (error,optimumParameters) = self.steepestDescentMinimum(
            self.evaluateObjectiveFunction,
            list(self.optimizedDesign),
            nMax = self.maxIterations,
            damping=self.damping)

        # Return to main thread
        #progress = float(n/nMax)
        self.optimizedDesign = optimumParameters
        self.iterationDone.emit(self.optimizedDesign)
        #rp = rp*1.5
        self.objectiveSequence.append(error)

        return self.optimizedDesign

    def evaluateObjectiveFunction(self,design,rp=1):

        objective = 0.0
        for target in self.targets:
            # Calculate position at crank angle target[0]
            position = self.calculateFourBarPoint(design,target[0])
            if position == [99999,99999]:
                objective += 999
            else:
                objective = objective + self.vLen(target[1:],position)**2
            #print(target)
        # if fails graphoff conditions, add large penalty value
        #penalty = 0

        #epsilon = -0.2*np.sqrt(1/rp)

        #objective = objective + penalty*rp
        #print(objective)
        return objective

    def calculateFourBarPoint(self,design,alphaInput):
        # lengths is in format [base, z1,z2,z3,z4,z5]
        # bases are in format [[base1X,base1Y],[base2X,base2Y]]
        # alpha is a single float for crank angle from horizontal

        bases = [[design[5],design[6]],
                [design[7],design[8]]]
        length = list(design[0:5])
        length.insert(0,self.vLen(bases[0],bases[1]))
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
        theta3ArccosValue = (x3**2 + y3**2 + (-length[3])**2 - length[5]**2)/2 * (-length[3]) * math.sqrt(x3**2 + y3**2)
        if abs(theta3ArccosValue) > 1:
            #print('Invalid theta3ArccosValue of ' + repr(theta3ArccosValue) + ' encountered. Crashing.')
            return [99999,99999]
        theta3ArccosValue = np.float64(theta3ArccosValue)
        theta3pos = math.atan2(y3,x3) + np.arccos(theta3ArccosValue)
        #print(theta3pos)
        theta3neg = math.atan2(y3,x3) - np.arccos(theta3ArccosValue)

        theta3 = theta3pos
        theta5 = math.atan2((y3-(-length[3])*np.sin(theta3))/length[5], (x3-(-length[3])*np.cos(theta3))/length[5])
        point2 = [0,0]
        point2[0] = bases[1][0] + length[3]*np.cos(theta3)
        point2[1] = bases[1][1] + length[3]*np.sin(theta3)

        dyad1_acos_value = np.float64(length[5]**2 + length[2]**2 - length[4]**2)/np.float64(2*length[5]*length[2])
        if abs(dyad1_acos_value) > 1:
            return [99999,99999]
        dyadAngle1 = np.arccos(dyad1_acos_value)
        pointC = [0,0]
        pointC[0] = point1[0]+length[2]*np.cos(theta5 + dyadAngle1)
        pointC[1] = point1[1]+length[2]*np.sin(theta5 + dyadAngle1)
        #print(pointC)
        return pointC

    def vLen(self,point1,point2):
        # takes in two points in the format [x,y] and returns the float of vector length
        dx = point1[0] - point2[0]
        dy = point1[1] - point2[1]

        length = np.sqrt(dx*dx + dy*dy)
        return length

    def threePointQuadraticApprox(self,x,y):
        C2 = (((y[2]-y[0])/(x[2]-x[0])) - ((y[1]-y[0])/(x[1]-x[0])))/(x[2]-x[1])
        C1 = (y[1] - y[0])/(x[1]-x[0]) - C2*(x[0]+x[1])
        C0 = y[0] - C1*x[0] - C2*x[0]**2

        return [C0,C1,C2]

    def minimizeParabola(self,c):
        # Inputs: Coefficients for polynomial equation according to the form C0 + C1*x + C2*x^2...
        # Outputs: Values of x and y where y is minimized
        minX = -c[1]/(2*c[2])

        minY = self.getValueOfPoly(c,minX)
        return (minX,minY)

    def getValueOfPoly(self,c,x):
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

    def steepestDescentMinimum(self,function,
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
            slopeList = self.gradient(function,position)
            # print("fitting polynomial...")
            # Get three points in that direction at positions of alpha
            functionValues = []
            for alphaValue in alpha:
                testLocation = []
                for oldPosition, slope in zip(position,slopeList):
                    testLocation.append(oldPosition-slope*alphaValue)
                functionValues.append(function(testLocation))
            # Fit parabola to curve
            C = self.threePointQuadraticApprox(alpha, functionValues)
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
                (alphaStar,bestY) = self.minimizeParabola(C)
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

    def gradient(self,function,inputs,delta=0.0001,normalize=False):
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


    def threePointCubicApprox(self,x,y,xSlopePoint,yPrime):
        C3 = (y[2] - y[0])/((x[2] - x[1])*(x[2] - x[0])**2) - (y[1] - y[0])/((x[2] - x[1])*(x[1] - x[0])**2) + yPrime/((x[1]-x[0])*(x[2] - x[0]))
        C2 = (((y[1] - y[0])/(x[1] - x[0])) - yPrime)/(x[1] - x[0]) - C3*(2*x[0] + x[1])
        C1 = yPrime - 2*C2*x[0] - 3*C3*x[0]**2
        C0 = y[0] - C1*x[0] - C2*x[0]**2 - C3*x[0]**3

        return [C0,C1,C2,C3]

    def threePointQuadraticApprox(self,x, y):
        #Inputs: Vector of x values and y values. These vectors must be equal length
        #Outputs: Coefficients for polynomial equation according to the form C0 + C1*x + C2*x^2...
        C2 = (((y[2]-y[0])/(x[2]-x[0])) - ((y[1]-y[0])/(x[1]-x[0])))/(x[2]-x[1])
        C1 = (y[1] - y[0])/(x[1]-x[0]) - C2*(x[0]+x[1])
        C0 = y[0] - C1*x[0] - C2*x[0]**2

        return [C0,C1,C2]

    def twoPointLinearApprox(self,x, y):
        #Inputs: Vector of x values and y values. These vectors must be equal length
        #Outputs: Coefficients for polynomial equation according to the form C0 + C1*x + C2*x^2...
        C1 = (y[1] - y[0])/(x[1]-x[0])
        C0 = y[0] - C1*x[0]

        return [C0,C1]

    def getValueOfPoly(self,c,x):
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

    def run(self):
        self.optimizedDesign = self.runOptimization()
        self.designOptimized.emit(self.optimizedDesign)

        self.sleep(2)
