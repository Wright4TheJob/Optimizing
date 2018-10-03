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

    def __del__(self):
        self.wait()

    def runOptimization(self, design):
        # Design parameter list creation
        designInitialPosition = design

        # Constraint functions
        # Dyad validity analysis
        # TODO: Add penalty function for invalid dyad

        rp = 1
        rpMax = 1000000

        optimumParameters = list(designInitialPosition)
        n = 34.
        nMax = self.maxIterations
        while rp < rpMax:
            n = n+1.0
            (error,optimumParameters) = self.numericalMin(optimumParameters,
                rp=rp,
                echo=False,
                damping=self.damping,
                epsilon=0.00001,
                nMax=self.maxIterations,
                alpha = [0,0.01,0.02],
                printResults=False)

            # Return to main thread
            progress = float(n/nMax)
            self.optimizedDesign = optimumParameters
            self.iterationDone.emit(self.optimizedDesign)
            rp = rp*1.5
            self.objectiveSequence.append(error)

        return self.optimizedDesign

    def evaluateObjectiveFunction(self,design,rp=1):

        for target in self.targets:
            # Calculate position at crank angle target[0]
            # error = error + (target[1]-position[0])**2 + (target[2]-position[1])**2
            print(target)
        # if fails graphoff conditions, add large penalty value
        penalty = 0

        epsilon = -0.2*np.sqrt(1/rp)

        objective = 0

        objective = objective + penalty*rp

        return objective

    def calculateFourBarPoint(self,alphaInput):
        # lengths is in format [base, z1,z2,z3,z4,z5]
        # bases are in format [[base1X,base1Y],[base2X,base2Y]]
        # alpha is a single float for crank angle from horizontal

        bases = list(self.mechanismBases)
        length = list(self.mechanismLengths)
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
        theta3ArccosValue = (x3**2 + y3**2 + (-length[3])**2 - length[5]**2)/2 * (-length[3]) * math.sqrt(x3*x3 + y3*y3)
        theta3ArccosValue = np.float64(theta3ArccosValue)
        theta3pos = math.atan2(y3,x3) + np.arccos(theta3ArccosValue)
        #print(theta3pos)
        theta3neg = math.atan2(y3,x3) - np.arccos(theta3ArccosValue)

        theta3 = theta3pos
        theta5 = math.atan2((y3-(-length[3])*np.sin(theta3))/length[5], (x3-(-length[3])*np.cos(theta3))/length[5])
        point2 = [0,0]
        point2[0] = bases[1][0] + length[3]*np.cos(theta3)
        point2[1] = bases[1][1] + length[3]*np.sin(theta3)

        dyadAngle1 = np.arccos(np.float64(length[5]**2 + length[2]**2 - length[4]**2)/np.float64(2*length[5]*length[2]))
        pointC = [0,0]
        pointC[0] = point1[0]+length[2]*np.cos(theta5 + dyadAngle1)
        pointC[1] = point1[1]+length[2]*np.sin(theta5 + dyadAngle1)

        return pointC

    def vLen(self,point1,point2):
        # takes in two points in the format [x,y] and returns the float of vector length
        dx = point1[0] - point2[0]
        dy = point1[1] - point2[1]

        length = np.sqrt(dx*dx + dy*dy)
        return length

    def getNumGradient(self,design,rp = 1,delta = 0.0001):
        testDesign = list(design)
        for v in range(0,len(design)):
            newVariableValue = design[v] - delta/2.0
            testDesign[v] = newVariableValue

        startingValue = self.evaluateObjectiveFunction(testDesign, rp = rp)
        slopeList = []
        for v in range(0,len(design)):
            newVariableValue = design[v] + delta/2.0
            testDesign = list(design)
            testDesign[v] = newVariableValue
            testObjectiveValue = self.evaluateObjectiveFunction(testDesign, rp = rp)
            slopeList.append((testObjectiveValue-startingValue)/delta)

        return slopeList

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

    def numericalMin(self,design,rp=1,echo=False,damping=0.1,epsilon=0.001,nMax=100,alpha = [0,0.01,0.02],printResults=False):

        # Constants
        alphaStarMax = 10
        deltaMax = 0.25 + 100/rp
        i = 0
        oscillationEpsilon = epsilon/100

        # Loop
        shouldContinue = True
        position = list(design)
        objectiveValue = self.evaluateObjectiveFunction(position,rp=rp)
        lastPosition = position
        if echo == True:
            headerString = "Iteration\t"
            headerString += "Design\t\t\t"
            #headerString += "Gradient\t"
            headerString += "F(x)"
            print(headerString)

        while shouldContinue == True:
            i = i+1

            #print("Total Iterations should be  %i" %(nMax))
            #print("Iteration %i" %(i))
            # Get gradient at position
            # print("About to get gradient")

            slopeList = self.getNumGradient(design,rp=rp)
            #print("Slope values from getGradient")
            #print(slopeList)
            # print("About to fit polynomial")
            # Get three points in that direction at intervals of 0.5,1,2
            functionValues = [objectiveValue]
            for alphaValue in alpha:
                if alphaValue != alpha[0]:
                    testLocation = []
                    for oldPosition, slope in zip(position,slopeList):
                        testLocation.append(oldPosition-slope*alphaValue)
                    functionValues.append(self.evaluateObjectiveFunction(testLocation,rp=rp))

            # Fit parabola to curve
            C = self.threePointQuadraticApprox(alpha, functionValues)
            # Check parabola is concave up
            # Calculate alpha that gives minimum
            alphaStar = 0.0
            if C[2] < 0:
                if echo == True:
                    print("Fitted parabola is concave down. Minimum alpha value is not bounded.")
                alphaStar = 0.1
            elif abs(C[2]) < 0.0001:
                alphaStar = 0.1
            else:
                (alphaStar,bestY) = self.minimizeParabola(C)

            if alphaStar > alphaStarMax:
                alphaStar = alphaStarMax
            # Move to position of calculated alpha
            newPosition = []
            for oldPosition, slope in zip(position,slopeList):
                moveDist = slope*damping*alphaStar
                if abs(moveDist) < deltaMax:
                    newPosition.append(oldPosition-moveDist)
                else:
                    newPosition.append(oldPosition-moveDist/abs(moveDist)*deltaMax)

            twoPositionsAgo = lastPosition
            lastPosition = position
            position = newPosition
            objectiveValueLast = objectiveValue
            objectiveValue = self.evaluateObjectiveFunction(position,rp=rp)

            # Print current iteration results
            if echo == True:
                resultsString = "%i        \t" %(i)
                resultsString += "{}\t".format(position)
                resultsString += "%2.6f" % (objectiveValue)
                print(resultsString)

            # Check convergence
            deltaObjective = objectiveValueLast - objectiveValue
            variableChanges = [abs(x - y) for x, y in zip(position,lastPosition)]
            deltaVar = float(max(variableChanges))
            #print("Delta Objective = %2.4f" % (float(deltaObjective)))
            if abs(float(deltaObjective)) < epsilon and i > 1 and deltaVar < epsilon:
                shouldContinue = False
                if printResults == True:
                    print("Local Optimium found")

            # Check oscillation
            oscillationDelta = [abs(x - y) for x, y in zip(position,twoPositionsAgo)]
            oscillationDelta = max(oscillationDelta)

            if oscillationDelta < oscillationEpsilon:
                damping = damping*0.75

            #print("About to check iteration maximum")
            if i > nMax:
                if printResults == True:
                    print("Function timed out. Returning final result")
                shouldContinue = False

        if printResults==True:
            print("#### - Results - ####")
            for variableValue in position:
                print("%2.6f" % (variableValue))
            print("F = %2.6f" % (objectiveValue))

        return (objectiveValue, position)

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
        self.design = self.runOptimization(self.design)
        self.designOptimized.emit(self.design)

        self.sleep(2)
