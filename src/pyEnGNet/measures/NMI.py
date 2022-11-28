from preprocess import *
import math

class NMI:
    
    @staticmethod
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
    def __calculateMutualInformation(gen1, gen2, size, maxVal):
        LOG_BASE = 2.0
        probMap = [0.0] * (8)
        # 2 (gen1) + 2 (gen2) + 4 (joint)
        iColumn = 0
        while (iColumn < size):
            valGen1Column = gen1[iColumn]
            valGen2Column = gen2[iColumn]
            probMap[valGen1Column] = probMap[valGen1Column] + 1
            probMap[valGen2Column + 2] = probMap[valGen2Column + 2] + 1
            probMap[(valGen1Column + maxVal * valGen2Column) +
                    4] = probMap[(valGen1Column + maxVal * valGen2Column) + 4] + 1
            iColumn += 1
        iCont = 0
        while (iCont < 8):
            probMap[iCont] = probMap[iCont] / size
            iCont += 1
        nMI = 0.0
        iCont = 0
        while (iCont < 4):
            if (probMap[iCont + 4] > 0.0 and probMap[iCont % maxVal] > 0.0 and probMap[
                (int(iCont / maxVal)) + 2] > 0.0):
                nMI += probMap[iCont + 4] * math.log(
                    probMap[iCont + 4] / probMap[iCont % maxVal] / probMap[(int(iCont / maxVal)) + 2])
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
        
        
        print("calculateEntropy G1: ",NMI.__calculateEntropy(gen1Normalized, maxValGene1))
        print("calculateEntropy G2: ",NMI.__calculateEntropy(gen2Normalized, maxValGene2))
        try:
            value = 2.0 * NMI.__calculateMutualInformation(gen1Normalized, gen2Normalized, size, maxValGene1) / (
                    NMI.__calculateEntropy(gen1Normalized,maxValGene1) + NMI.__calculateEntropy(gen2Normalized,maxValGene2))
        except Exception as e:
            value = 0.0
        return value
    
    """
        public static strictfp double calculateMutualInformation(final double[] array, final double[] array2) {
        final JointProbabilityState jointProbabilityState = new JointProbabilityState(array, array2);
        final int firstMaxVal = jointProbabilityState.firstMaxVal;
        double n = 0.0;
        for (final Integer key : jointProbabilityState.jointProbMap.keySet()) {
            final double doubleValue = jointProbabilityState.jointProbMap.get(key);
            final double doubleValue2 = jointProbabilityState.firstProbMap.get(key % firstMaxVal);
            final double doubleValue3 = jointProbabilityState.secondProbMap.get(key / firstMaxVal);
            if (doubleValue > 0.0 && doubleValue2 > 0.0 && doubleValue3 > 0.0) {
                n += doubleValue * Math.log(doubleValue / doubleValue2 / doubleValue3);
            }
        }
        return n / Math.log(Entropy.LOG_BASE);
    }
    """
    
    @staticmethod
    def process(dataset, arr1, arr2):
        corr = NMI.__calculationNMI(arr1, arr2, len(arr1))
        if math.isnan(corr):
            corr = 0
        ans = 0
        if corr >= dataset.nmi_threshold:
            ans = 1
        return ans, corr