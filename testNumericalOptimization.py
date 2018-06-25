# Tests numerical optimization functions
import NumericalOptimization as numOpt
import math

def sinFunction(inputs):
    x = inputs[0]
    result = math.sin(x)
    return result

def parabola(inputs):
    x = float(inputs[0])
    result = x*x
    return result

def multiDimFunction(inputs):
    x = inputs[0]
    y = inputs[1]
    return x + y*y

def upperBound(inputs):
    x = inputs[0]
    return x-4

def lowerBound(inputs):
    x = inputs[0]
    return 2-x

#slopeList = numOpt.gradient(parabola,[3])
#print(slopeList)



optimizedLocation = numOpt.constrainedMinimum(sinFunction,[2.5],inequalityConstraints = [upperBound,lowerBound],damping=.1,echo=True,nMax=2000,rp=100,method='InteriorLinearExtended')
print(optimizedLocation)
