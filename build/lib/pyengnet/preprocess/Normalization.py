import math

class Normalization:
    
    def normalizedArray(gen, genNormalized, size):
        maxValue = 0
        if (size > 0):
            minValue = int(math.floor(gen[0]))
            maxValue = int(math.floor(gen[0]))
            iCont = 0
            while (iCont < size):
                iExp = int(math.floor(gen[iCont]))
                genNormalized[iCont] = iExp
                if (iExp < minValue):
                    minValue = iExp
                if (iExp > maxValue):
                    maxValue = iExp
                iCont += 1
            iCont = 0
            while (iCont < size):
                genNormalized[iCont] -= minValue
                iCont += 1
            maxValue = maxValue - minValue + 1
        return maxValue