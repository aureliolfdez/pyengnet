from .Normalization import Normalization

import math

class NMI:
    
    """
        NMI measurement class coded in a parallel ecosystem with CPUs.
    """
    
    def __calculateEntropy(genNormalized, maxValGene):
        LOG_BASE = 2.0
        probMap = [0.0] * (maxValGene)
        iColumn = 0
        size = len(genNormalized)
        
        while (iColumn < size):
            iExpr = genNormalized[iColumn]
            probMap[iExpr] = probMap[iExpr] + 1
            iColumn += 1
        
        iCont = 0        
        while (iCont < maxValGene):
            probMap[iCont] = probMap[iCont] / size
            iCont += 1
        
        dEntropy = 0.0                                      # n = 0.0
        iCont = 0
        while (iCont < maxValGene):
            varAux = probMap[iCont]
            if (varAux > 0.0):                              # if (n2 > 0.0)
                dEntropy -= varAux * math.log(varAux)           # n -= n2 * Math.log(n2)
            iCont += 1
            
        dEntropy /= math.log(LOG_BASE) # n / Math.log(Entropy.LOG_BASE)
        
        return dEntropy
    
    @staticmethod
    def __calculateMutualInformation(gen1, gen2, size, maxValGene1, maxValGene2):
        LOG_BASE = 2.0
        total = maxValGene1 + maxValGene2 + (maxValGene1 * maxValGene2)
        probMap = [0.0] * (total)
        iColumn = 0
        while (iColumn < size):
            valGen1Column = gen1[iColumn]
            valGen2Column = gen2[iColumn]
            
            probMap[valGen1Column] = probMap[valGen1Column] + 1
            probMap[valGen2Column + maxValGene1] = probMap[valGen2Column + maxValGene1] + 1            
            probMap[(valGen1Column + maxValGene1 * valGen2Column) + (maxValGene1 + maxValGene2)] = probMap[(valGen1Column + maxValGene1 * valGen2Column) + (maxValGene1 + maxValGene2)] + 1
            iColumn += 1
                
        iCont = 0
        while (iCont < total):
            if probMap[iCont] > 0.0:
                probMap[iCont] = probMap[iCont] / size
            iCont += 1
        
        nMI = 0.0
        iCont = 0
        subtotal = maxValGene1 * maxValGene2
        while (iCont < subtotal):
            if probMap[iCont + (maxValGene1 + maxValGene2)] > 0.0:
                doubleValue = probMap[iCont + (maxValGene1 + maxValGene2)]
                doubleValue2 = probMap[(iCont % maxValGene1)]
                doubleValue3 = probMap[(int)(((iCont / maxValGene1) + maxValGene1))]
                if (doubleValue > 0.0 and doubleValue2 > 0.0 and doubleValue3 > 0.0):
                    nMI += doubleValue * math.log(doubleValue / doubleValue2 / doubleValue3)
            
            iCont += 1
        nMI /= math.log(LOG_BASE)
        
        return nMI
    
    
    @staticmethod
    def __calculationNMI(gen1, gen2, size):
        value = 0.0
        # Normalized arrays
        gen1Normalized = [0] * (len(gen1))
        gen2Normalized = [0] * (len(gen2))
        maxValGene1 = Normalization.normalizedArray(gen1, gen1Normalized, len(gen1))
        maxValGene2 = Normalization.normalizedArray(gen2, gen2Normalized, len(gen2))
        
        try:
            value = 2.0 * NMI.__calculateMutualInformation(gen1Normalized, gen2Normalized, size, maxValGene1, maxValGene2) / (
                    NMI.__calculateEntropy(gen1Normalized,maxValGene1) + NMI.__calculateEntropy(gen2Normalized,maxValGene2))
        except Exception as e:
            value = 0.0
        return value
    
    @staticmethod
    def process(dataset, arr1, arr2):
        """
        
        Function that performs NMI's measure for two genes.
                
        Parameters
        ----------
        arr1: numpy.array
            Array storing all the values of the dataset for a gene X.
        
        arr2: numpy.array
            Array storing all the values of the dataset for a gene Y.
        
        Return : int, double
        -------------------------
            ans: int
                This value indicates whether the calculated NMI's coefficient is valid according to the threshold offered by the user. If this variable stores the value 1, it is valid, whereas it is invalid if the stored value is 0.
            
            corr: double
                NMI's coefficient calculated.
            
        """
        corr = NMI.__calculationNMI(arr1, arr2, len(arr1))
        if math.isnan(corr):
            corr = 0
        ans = 0
        if corr >= dataset.nmi_threshold:
            ans = 1
        return ans, corr