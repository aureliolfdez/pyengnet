from .models import *

import pandas as pd
import numpy as np
import multiprocessing as mp
import os
import matplotlib.pyplot as plt
import networkx as nx

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
            gene = dataset[:,0]
            dataset = np.delete(dataset, 0, 1)
            dataset = Dataset(data=dataset, gene=gene, nmi_th=nmi_th, spearman_th=spearman_th, kendall_th=kendall_th, readded_th=readded_th, cores=cores)  
    return dataset

def saveFile(path = '', graph = None):
    
    if graph == None or len(graph) > 0:
        with open(path+'graphComplete.csv', 'w') as fp:
            for iLine in graph:
                fp.write("%s\n" % iLine)
        print('INFO: Graph complete exported to file')
    
def showGraph(graph,title=None):
    
    graphComplete = nx.Graph()
    graphComplete.add_edges_from(graph)
    
    pos = nx.spring_layout(graphComplete,k=100)
    nx.draw_networkx_nodes(graphComplete, pos, cmap=plt.get_cmap('jet'), node_size = 500)
    nx.draw_networkx_labels(graphComplete, pos)
    nx.draw_networkx_edges(graphComplete, pos, arrows=False)
    nx.draw_networkx_edge_labels(graphComplete, pos, edge_labels=nx.get_edge_attributes(graphComplete, 'weight'))
    plt.title(title)
    plt.show()