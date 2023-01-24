import scipy.stats as scp
import numpy as np
import warnings

class Spearman:
    
    def process(dataset, arr1, arr2):
        """
        Recibe dos vectores y calcula el coeficiente de correlación de spearman entre ellos, devuelve:
            ans -> 1 si pasa el umbral, 0 en caso contrario
            corr -> El coeficiente de correlación de spearman en valor absoluto
        :param arr1:
        :param arr2:
        :return:
        """
        ans = 0
        corr = 0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            corr, pv = scp.spearmanr(arr1, arr2, nan_policy='omit')
            if corr is np.nan:
                corr = 0
        corr = abs(corr)
        if corr >= dataset.spearman_threshold:
            ans = 1
        return ans, corr