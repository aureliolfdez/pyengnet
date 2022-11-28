from base import *
from algorithms import *
import numpy as np

if __name__ == "__main__":
    filePath = os.getcwd()+"/src/pyEnGNet/datasets/Spellman_v2.csv"
    dataset = load(path=filePath, separator=",", nmi_th=0.5, spearman_th=0.5, kendall_th=0.5, readded_th=0.5)    
    G, aristas = Engnet.process(dataset)
    #nx.draw_kamada_kawai(G)