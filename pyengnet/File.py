from .Dataset import Dataset

import pandas as pd
import numpy as np
import multiprocessing as mp
import os
import matplotlib.pyplot as plt
import networkx as nx

class File:
    
    """
        Class used to manage file I/O operations and data visualization.
    """
    
    def load(path, separator = "\t", nmi_th=0.7, spearman_th=0.7, kendall_th=0.7, readded_th=0.7, hub_th=2, cores=int(mp.cpu_count() / 2)):
        """
        Load data from a txt or csv file.
        
        Parameters
        ----------
        path : str
            The path where the file is stored.
        
        separator: str, optional (default is a tab)
            An attribute indicating how the columns of the file are separated. 
        
        nmi_th: float, optional (default=0.7)
            Threshold for the NMI classifier (Used in phase 1 of the EnGNet algorithm).
        
        spearman_th: float, optional (default=0.7)
            Threshold for the Spearman classifier (Used in phase 1 of the EnGNet algorithm).
        
        kendall_th: float, optional (default=0.7)
            Threshold for the Kendall classifier (Used in phase 1 of the EnGNet algorithm).
        
        readdded_th: float, optional (default=0.7)
            Threshold to determine if the edge would be return into the network after the pruning step
        
        hub_th: int, optional (default=2)
            Threshold to determine if the node studie is a hub. Set this threshold to -1 to run the algorithm with standard selection.
        
        cores: int, optional (default = CPU cores / 2)
            Number of CPU cores used for parallelisation.
        
        Return : numpy.Array
        -------------------------
            A NumPy array from the reading of a file.
            
        """
        
        dataset = None
        if path is not None:
            extensionsCsv = [".txt",".csv"]    
            fileName, fileExtension = os.path.splitext(path)
            if fileExtension in extensionsCsv:
                dataset = np.asarray(pd.read_csv(path, sep=separator).dropna())
                gene = dataset[:,0]
                dataset = np.delete(dataset, 0, 1)
                dataset = Dataset(filePath = path, data=dataset, gene=gene, nmi_th=nmi_th, spearman_th=spearman_th, kendall_th=kendall_th, readded_th=readded_th, hub_th=hub_th, cores=cores)    
        return dataset

    def saveFile(path = None, graph = None):
        
        """
        Save a complete and/or pruned network in an external file.
        
        Parameters
        ----------
        path : str
            The path where the file will be stored.
        
        graph: Networkx object (https://pypi.org/project/networkx/)
            Object of type networkx to be stored in a file. This object can be a complete network or a pruned network.
        
        """
        
        if (graph != None or len(graph) > 0) and path != None:
            with open(path, 'w') as fp:
                for iLine in graph:
                    fp.write("%s\n" % iLine)
            print('INFO: Graph complete exported to file')
        
    def showGraph(graph = None, title = ''):
        
        """
        Display a complete and/or pruned network.
        
        Parameters
        ----------
       
        graph: Networkx object (https://pypi.org/project/networkx/)
            Object of type networkx to be stored in a file. This object can be a complete network or a pruned network.
            
        title : str, optional (default='')
            The main title of the network
        
        """
        if graph != None or len(graph) > 0:
            pos = nx.spring_layout(graph,k=100)
            nx.draw_networkx_nodes(graph, pos, cmap=plt.get_cmap('jet'), node_size = 500)
            nx.draw_networkx_labels(graph, pos)
            nx.draw_networkx_edges(graph, pos, arrows=False)
            nx.draw_networkx_edge_labels(graph, pos, edge_labels=nx.get_edge_attributes(graph, 'weight'))
            plt.title(title)
            plt.show()