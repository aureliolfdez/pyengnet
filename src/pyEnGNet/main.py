from base import *
from algorithms import *

if __name__ == "__main__":
    filePath = os.getcwd()+"/src/pyEnGNet/datasets/113_exp_mat_cond_1.csv"
    dataset = load(path=filePath, separator=",")    
    G, aristas = Algorithm.engnet(dataset)
    #nx.draw_kamada_kawai(G)