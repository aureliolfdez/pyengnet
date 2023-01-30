import os

from pyengnet.File import File
from pyengnet.Engnet import Engnet

if __name__ == "__main__":
    
    dataset = File.load(path=os.getcwd()+"/datasets/Spellman.csv", separator=",", nmi_th=0.6, spearman_th=0.7, kendall_th=0.7, readded_th=0.7, hub_th = 3)    

    # Run pyEnGNet
    graphFiltered, infoGraphFiltered, graphComplete, infoGraphComplete = Engnet.process(dataset, saveComplete = True, numGpus = 1, computeCapability = 61)
    
    # Save files
    File.saveFile(path='/home/principalpc/Escritorio/graphComplete.csv',graph=infoGraphComplete)
    File.saveFile(path='/home/principalpc/Escritorio/graphFiltered.csv',graph=infoGraphFiltered)
    
    # Print graphs examples
    #File.showGraph(graph=graphComplete,title='Complete graph')
    #File.showGraph(graph=graphFiltered,title="Filtered graph")
    
    