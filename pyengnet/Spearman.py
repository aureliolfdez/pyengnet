import scipy.stats as scp
import numpy as np
import warnings

class Spearman:
    
    """
        Spearman measurement class coded in a parallel ecosystem with CPUs.
    """
    
    def process(dataset, arr1, arr2):
        """
        
        Function that performs Spearman's measure for two genes.
                
        Parameters
        ----------
        arr1: numpy.array
            Array storing all the values of the dataset for a gene X.
        
        arr2: numpy.array
            Array storing all the values of the dataset for a gene Y.
        
        Return : int, double
        -------------------------
            ans: int
                This value indicates whether the calculated Spearman's coefficient is valid according to the threshold offered by the user. If this variable stores the value 1, it is valid, whereas it is invalid if the stored value is 0.
            
            corr: double
                Spearman's coefficient calculated.
            
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