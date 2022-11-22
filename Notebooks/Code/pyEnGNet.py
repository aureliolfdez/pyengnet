import concurrent.futures

import pandas as pd
import scipy.stats as scp
import numpy as np
from tqdm import tqdm
import networkx as nx
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
import math

sameval = 0.7


def intervals(total, partes):
    lenpartes = int(total / partes)
    lista = []
    for i in range(partes):
        if i == partes - 1:
            end = total
        else:
            end = lenpartes * (i + 1)
        lista.append((lenpartes * i, end))
    return lista


class PyEnGNet:
    def __init__(self, nparr=None, nmi_th=sameval, spearman_th=sameval, kendall_th=sameval, pearson_th=sameval,
                 readded_th=0.7, cores=int(mp.cpu_count() / 2)):
        if nparr is not None:
            self.maindata = nparr
        else:
            raise Exception("No data given")

        self.row_size = len(self.maindata)
        self.column_size = len(self.maindata[0])

        self.spearman_threshold = spearman_th
        self.kendall_threshold = kendall_th
        self.nmi_threshold = nmi_th
        self.pearson_threshold = pearson_th

        self.readded_edges_threshold = readded_th
        self.ncores = cores

    def calculateEntropy(self, gen, genNormalized, size):
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

    def normalizedArray(self, gen, genNormalized, size):
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

    def calculateMutualInformation(self, gen1, gen2, size, maxVal):
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

    def calculationNMI(self, gen1, gen2, size):
        value = 0.0
        # Normalized arrays
        gen1Normalized = [0] * (len(gen1))
        gen2Normalized = [0] * (len(gen2))
        maxVal = self.normalizedArray(gen1, gen1Normalized, len(gen1))
        self.normalizedArray(gen2, gen2Normalized, len(gen2))
        try:
            value = 2.0 * self.calculateMutualInformation(gen1Normalized, gen2Normalized, size, maxVal) / (
                    self.calculateEntropy(gen1, gen1Normalized, size) + self.calculateEntropy(gen2, gen2Normalized,
                                                                                              size))
        except Exception as e:
            value = 0.0
        return value

    def single_nmi(self, arr1, arr2):
        corr = self.calculationNMI(arr1, arr2, len(arr1))
        if math.isnan(corr):
            corr = 0
        ans = 0
        if corr >= self.nmi_threshold:
            ans = 1
        return ans, corr

    def single_spearman(self, arr1, arr2):
        """
        Recibe dos vectores y calcula el coeficiente de correlación de spearman entre ellos, devuelve:
            ans -> 1 si pasa el umbral, 0 en caso contrario
            corr -> El coeficiente de correlación de spearman en valor absoluto
        :param arr1:
        :param arr2:
        :return:
        """
        ans = 0
        corr, pv = scp.spearmanr(arr1, arr2)
        corr = abs(corr)
        if corr >= self.spearman_threshold:
            ans = 1
        return ans, corr

    def single_pearson(self, arr1, arr2):
        """
        Recibe dos vectores y calcula el coeficiente de correlación de spearman entre ellos, devuelve:
            ans -> 1 si pasa el umbral, 0 en caso contrario
            corr -> El coeficiente de correlación de spearman en valor absoluto
        :param arr1:
        :param arr2:
        :return:
        """
        ans = 0
        corr, pv = scp.pearsonr(arr1, arr2)
        corr = abs(corr)
        if corr >= self.pearson_threshold:
            ans = 1
        return ans, corr

    def single_kendall(self, arr1, arr2):
        """
        Recibe dos vectores y calcula el coeficiente de correlación de Kendall entre ambos
            ans -> 1 si pasa el umbral, 0 en caso contrario
            corr -> El valor resultante, en valor absoluto, del coeficiente de correlación de kendall
        :param arr1:
        :param arr2:
        :return:
        """
        ans = 0
        corr, pv = scp.kendalltau(arr1, arr2)
        corr = abs(corr)
        if corr >= self.kendall_threshold:
            ans = 1
        return ans, corr

    def validate_corr(self, i, j, accepted_values):
        major_voting = 0

        # Las dos filas que vamos a utilizar
        v = self.maindata[i]
        w = self.maindata[j]

        # Agregamos las respuestas de los tests a una lista que servirá para el calculo de los pesos
        tests = [self.single_nmi(v, w), self.single_kendall(v, w), self.single_spearman(v, w)]

        for test in tests:
            major_voting += test[0]

        if major_voting >= 2:
            accepted_values.append(
                [i, j, {'weight': self.calculate_weight(tests)}])

    def calculate_weight(self, tests):
        weight = []
        for test in tests:
            if bool(test[0]):
                weight.append(test[1])
        return np.mean(weight)

    def edge_corr_validation(self, start, end):
        accepted = []
        for i in tqdm(range(start, end)):
            for j in range(i + 1, self.row_size):
                self.validate_corr(i, j, accepted)
        return accepted

    def mainmethod(self):
        intervalos = intervals(len(self.maindata), self.ncores)
        edges = []

        with ProcessPoolExecutor(max_workers=self.ncores) as executor:
            results = []
            for rango in intervalos:
                start = rango[0]
                end = rango[1]
                results.append(executor.submit(self.edge_corr_validation, start, end))

        for f in concurrent.futures.as_completed(results):
            for val in f.result():
                edges.append(val)

        return edges

    def engnet_1_0(self):
        oedges = self.mainmethod()
        G = nx.Graph()
        G.add_edges_from(oedges)
        G2 = nx.maximum_spanning_tree(G, weight='weight', algorithm="kruskal")

        G3 = self.readd_edges(G, G2)
        fedges = nx.to_edgelist(G3)

        return G3, fedges

    def readd_edges(self, graph0, graph2):
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
        detected_hubs = [graph2.degree(node) for node in graph2.nodes if graph2.degree(node) > 2]

        if len(detected_hubs) != 0:
            mean_degree = round(np.mean([graph2.degree(node) for node in graph2.nodes if graph2.degree(node) > 2]))
            hubs = []
            for node in graph2.nodes:
                if graph2.degree(node) >= mean_degree:
                    hubs.append(node)

            readded = []
            for node in hubs:
                for e in eliminated_edges:
                    if (node in e) and (e not in readded) and (e[2]['weight'] > self.readded_edges_threshold):
                        readded.append(e)
            G3 = nx.Graph()
            G3.add_edges_from(nx.to_edgelist(graph2))
            G3.add_edges_from(readded)
            return G3
        else:
            return graph2

    def __str__(self):
        s = f"PyEnGNet Object with shape ({self.row_size} rows,{self.column_size} columns)"
        s += f"\n\n\tSpearman threshold selected: {self.spearman_threshold}"
        s += f"\n\tKendall Tau threshold selected: {self.kendall_threshold}"
        s += f"\n\tNMI threshold selected: {self.nmi_threshold}"
        s += f"\n\n\tWeight Threshold selected for edge readition: {self.readded_edges_threshold}"
        s += f"\nUsing {self.ncores} of your cores to process results"

        return s


if __name__ == "__main__":
    df = pd.read_csv("/home/daiego/PycharmProjects/pyEnGNet/pyEnGNet/Notebooks/Data/113_exp_mat_cond_1.csv")
    df = df.drop(df.columns[[0, 2]], axis=1)
    data = df.to_numpy()
    peg = PyEnGNet(nparr=data)
    print(peg)
    G, aristas = peg.engnet_1_0()
