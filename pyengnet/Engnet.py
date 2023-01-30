from .Normalization import Normalization
from .Dataset import Dataset
from .Kendall import Kendall
from .NMI import NMI
from .Spearman import Spearman

import concurrent.futures
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from tqdm import tqdm
import multiprocessing as mp
import math
import subprocess
import os
import networkx as nx
import pandas as pd


class OutputEngnet:
    
    """
        Class used to save the complete network (phase 1 of the algorithm).        
    """
    
    def __init__(self, accepted=None, testsOutput = None):
       
        self._accepted = accepted
        self._testsOutput = testsOutput  
        
    @property
    def accepted(self):
        """
            Valid correlations after running phase 1 of the algorithm.
        """
        return self._accepted
    
    @accepted.setter
    def accepted(self, accepted):
        self._accepted = accepted
    
    @property
    def testsOutput(self):
        """
            String associated to the valid correlations to be used if the user wants to save the network in an external file.
        """
        return self._testsOutput
    
    @testsOutput.setter
    def testsOutput(self, testsOutput):
        self._testsOutput = testsOutput

class Engnet:
    
    """
        Main class of the EnGNet algorithm.
    """
    
    def __intervals(total, partes):
        lenpartes = int(total / partes)
        lista = []
        for i in range(partes):
            if i == partes - 1:
                end = total
            else:
                end = lenpartes * (i + 1)
            lista.append((lenpartes * i, end))
        return lista

    def __readd_edges(dataset, graph0, graph2):
        """
        Recibe el grafo inicial y el árbol resultante.
        1. Se obtienen cuáles han sido las aristas eliminadas
        2. Se detectan los hubs de la red
            2.a. Si no hay hubs en la red, entonces se devuelve el grafo original
            2.b. Si hay hubs, se volverán a añadir las aristas cuyo peso supere un cierto umbral, determinado por el usuario
        :param graph0:
        :param graph2:
        :return:
        """
        origin_edges = nx.to_edgelist(graph0)
        tree_edges = nx.to_edgelist(graph2)
        
        # Obtain eliminated edges
        eliminated_edges = [edge for edge in origin_edges if edge not in tree_edges]
        
        # Obtain hubs
        detected_hubs = [graph2.degree(node) for node in graph2.nodes if graph2.degree(node) > dataset.hub_threshold]
        
        if len(detected_hubs) != 0:
            mean_degree = round(np.mean([graph2.degree(node) for node in graph2.nodes if graph2.degree(node) > dataset.hub_threshold]))
            hubs = []
            for node in graph2.nodes:
                if graph2.degree(node) >= mean_degree:
                    hubs.append(node)

            readded = []
            for node in hubs:
                for e in eliminated_edges:
                    if (node in e) and (e not in readded) and (e[2]['weight'] > dataset.readded_edges_threshold):
                        readded.append(e)
            G3 = nx.Graph()
            G3.add_edges_from(nx.to_edgelist(graph2))
            G3.add_edges_from(readded)
            return G3
        else:
            return graph2
    
    def __validate_corr(dataset, i, j, accepted_values, testsOutput_values, saveComplete = False):
        major_voting = 0

        # Las dos filas que vamos a utilizar
        v = dataset.data[i]
        w = dataset.data[j]
        
        # Agregamos las respuestas de los tests a una lista que servirá para el calculo de los pesos
        tests = [NMI.process(dataset, v, w), Kendall.process(dataset, v, w), Spearman.process(dataset, v, w)]
                
        for test in tests:
            major_voting += test[0]

        if major_voting >= 2:
            if saveComplete == True:
                testsOutput_values.append(str(dataset.gene[i])+"\t"+str(dataset.gene[j])+"\t"+str(tests[0][1])+"\t"+str(tests[1][1])+"\t"+str(tests[2][1]))
            
            #accepted_values.append(
            #    [dataset.gene[i], dataset.gene[j], {'weight': round(Engnet.__calculate_weight(tests),2)}])
            
            media = (tests[0][1]+tests[1][1]+tests[2][1])/3.0
            accepted_values.append(
                [dataset.gene[i], dataset.gene[j], {'weight': round(media,4)}])
            
    
    def edge_corr_validation(dataset, start, end, saveComplete = False):
        
        """
        
        A function that represents what a CPU core should do.
                
        Parameters
        ----------
        dataset: numpy.array
            The dataset stored in memory.
        
        start: int
            Number of the row of the dataset to be started by a specific CPU core.
            
        end: int
            Number of the row of the dataset to be ended by a specific CPU core.
            
        saveComplete: boolean, optional (default=False)
            Flag indicating whether or not the user wants to store the networks in external files.
            
        """
        
        accepted = []
        testsOutput = []
        
        for i in tqdm(range(start, end)):
            for j in range(i + 1, dataset.row_size):
                Engnet.__validate_corr(dataset, i, j, accepted, testsOutput, saveComplete)
        
        return OutputEngnet(accepted,testsOutput)
    
    def __CPUmethod(dataset, saveComplete = False):
        intervalos = Engnet.__intervals(len(dataset.data), dataset.ncores)
        edges = []
        testsOutput = []
        
        if saveComplete == True:
            testsOutput.append("Source\tDestination\tNMI\tKendall\tSpearman")
            
        with ProcessPoolExecutor(max_workers=dataset.ncores) as executor:            
            results = []
            for rango in intervalos:
                start = rango[0]
                end = rango[1]
                results.append(executor.submit(Engnet.edge_corr_validation, dataset, start, end, saveComplete))
        
        for f in concurrent.futures.as_completed(results):
            for val in f.result().testsOutput:
                testsOutput.append(val)
                
            for val in f.result().accepted:
                edges.append(val)
        
        graphComplete = nx.Graph()
        graphComplete.add_edges_from(edges)

        return graphComplete, testsOutput
    
    def __GPUmethod(dataset, saveComplete = False, numGpus = 1, CC = 61):
        os.system('cd pyengnet && make clean')
        os.system('cd pyengnet && make CC="'+str(CC)+'"')
        os.system('./pyengnet/cuEnGNet --data '+str(dataset.filePath)+' --cors '+str(dataset.spearman_threshold)+' --cork '+str(dataset.kendall_threshold)+' --corn '+str(dataset.nmi_threshold)+' --gpu '+str(numGpus)+' --out')
        
        accepted = []
        testsOutput = []
        
        if saveComplete == True:
            testsOutput.append("Source\tDestination\tNMI\tKendall\tSpearman")
        
        for iGpu in range(numGpus):        
            resultsTest = np.asarray(pd.read_csv("results/results_GPU_"+str(iGpu)+".csv", sep=";"))
            for test in resultsTest:
                media = (test[2]+test[3]+test[4])/3.0
                accepted.append([test[0], test[1], {'weight': round(media,4)}])

                if saveComplete == True: 
                    testsOutput.append(str(test[0])+"\t"+str(test[1])+"\t"+str(test[2])+"\t"+str(test[3])+"\t"+str(test[4]))

        graphComplete = nx.Graph()
        graphComplete.add_edges_from(accepted)
        
        return graphComplete, testsOutput
        
    def process(dataset, saveComplete = False, numGpus = None, computeCapability = None):
        
        """
        
        Function that run EnGNet algorithm.
                
        Parameters
        ----------
        dataset: numpy.array
            The dataset stored in memory.
        
        saveComplete: boolean, optional (default=False)
            Flag indicating whether or not the user wants to store the networks in external files.
        
        numGpus: int, optional (default=None)
            Number of GPU devices to be used.
            
        computeCapability: int
            Compute Capability version of your GPU device (https://developer.nvidia.com/cuda-gpus).
            
        """
        
        # First step: Ensemble
        if(numGpus == None or numGpus < 1):
            graphComplete, infoGraphComplete = Engnet.__CPUmethod(dataset, saveComplete)
        else:
            graphComplete, infoGraphComplete = Engnet.__CPUmethod(dataset, saveComplete)
            #graphComplete, infoGraphComplete = Engnet.__GPUmethod(dataset, saveComplete, numGpus, computeCapability)
        
        # Second step: Pruning
        graphFiltered = nx.minimum_spanning_tree(graphComplete) # Kruskal        
        graphFiltered = Engnet.__readd_edges(dataset, graphComplete, graphFiltered)
        edgesFiltered = nx.to_edgelist(graphFiltered)
        
        # Filtered output
        infoGraphFiltered = []
        
        if saveComplete == True:
            infoGraphFiltered.append("Gene1\tGene2\tWeight")
            for edge in edgesFiltered:
                infoGraphFiltered.append(str(edge[0])+"\t"+str(edge[1])+"\t"+str(edge[2]["weight"]))
        
        return graphFiltered, infoGraphFiltered, graphComplete, infoGraphComplete


    
