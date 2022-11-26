from .models import *

import pandas as pd
import numpy as np
import multiprocessing as mp
import os

def load(path, separator = "\t", nmi_th=0.7, spearman_th=0.7, kendall_th=0.7, readded_th=0.7, cores=int(mp.cpu_count() / 2)):
    """
    Load data from a txt or csv file. (Reuse function)
    
    Parameters
    ----------
    path : str
        The path where the file is stored.
    
    separator: str
        An attribute indicating how the columns of the file are separated. 
        Default: tab.
    
    Return : pandas.DataFrame
    -------------------------
        A NumPy array from the reading of a file.
        
    """
    
    dataset = None
    if path is not None:
        extensionsCsv = [".txt",".csv"]    
        fileName, fileExtension = os.path.splitext(path)
        if fileExtension in extensionsCsv:
            dataset = np.asarray(pd.read_csv(path, sep=separator))
            dataset = np.delete(dataset, 0, 1)
            dataset = Dataset(dataset, nmi_th, spearman_th, kendall_th, readded_th, cores)        
    return dataset