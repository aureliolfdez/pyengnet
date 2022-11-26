from preprocess import *
import math

class NMI:
    
    @staticmethod
    def __calculateEntropy(gen, genNormalized, size):
        LOG_BASE = 2.0
        probMap = [0.0] * (2)
        iColumn = 0
        while (iColumn < size):
            iExpr = genNormalized[iColumn]
            probMap[iExpr] = probMap[iExpr] + 1
            iColumn += 1
        iCont = 0
        while (iCont < 2):
            probMap[iCont] = probMap[iCont] / size
            iCont += 1
        dEntropy = 0.0
        iCont = 0
        while (iCont < 2):
            varAux = probMap[iCont]
            if (varAux > 0.0):
                dEntropy -= varAux * math.log(varAux)
            iCont += 1
        dEntropy /= math.log(LOG_BASE)
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
        maxVal = Normalization.normalizedArray(gen1, gen1Normalized, len(gen1))
        Normalization.normalizedArray(gen2, gen2Normalized, len(gen2))
        try:
            print("\naqui esta: ")
            value = 2.0 * NMI.__calculateMutualInformation(gen1Normalized, gen2Normalized, size, maxVal) / (
                    NMI.__calculateEntropy(gen1, gen1Normalized, size) + NMI.__calculateEntropy(gen2, gen2Normalized,
                                                                                              size))
        except Exception as e:
            value = 0.0
        return value
    
    @staticmethod
    def process(dataset, arr1, arr2):
        corr = NMI.__calculationNMI(arr1, arr2, len(arr1))
        if math.isnan(corr):
            corr = 0
        ans = 0
        if corr >= dataset.nmi_threshold:
            ans = 1
        return ans, corr