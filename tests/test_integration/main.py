import os

from pyengnet.File import File
from pyengnet.Engnet import Engnet

if __name__ == "__main__":
    
    dataset = File.load(path=os.getcwd()+"/datasets/TCGA/expression_good.csv", separator=",", nmi_th=0.6, spearman_th=0.7, kendall_th=0.7, readded_th=0.7, hub_th = 3)    
    #dataset = File.load(path=os.getcwd()+"/datasets/Spellman.csv", separator=",", nmi_th=0.6, spearman_th=0.7, kendall_th=0.7, readded_th=0.7, hub_th = 3)    

    # Execute EnGNet
    graphFiltered, infoGraphFiltered, graphComplete, infoGraphComplete = Engnet.process(dataset, saveComplete = True, numGpus = 2, computeCapability = 61)
    #graphFiltered, infoGraphFiltered, graphComplete, infoGraphComplete = Engnet.process(dataset, saveComplete = True)
    
    # Save files info MacOSX
    #saveFile(path='/Users/aurelio/Desktop/graphComplete.csv',graph=infoGraphComplete)
    #saveFile(path='/Users/aurelio/Desktop/graphFiltered.csv',graph=infoGraphFiltered)
    
    # Save files info Ubuntu
    File.saveFile(path='/home/principalpc/Escritorio/graphComplete.csv',graph=infoGraphComplete)
    File.saveFile(path='/home/principalpc/Escritorio/graphFiltered.csv',graph=infoGraphFiltered)
    
    # Save files info Windows
    #saveWindowsFile = os.getcwd()+"/src/pyEnGNet/graphComplete.csv"
    #saveWindowsFilteredFile = os.getcwd()+"/src/pyEnGNet/graphFiltered.csv"
    #saveFile(path=saveWindowsFile,graph=infoGraphComplete)
    #saveFile(path=saveWindowsFilteredFile,graph=infoGraphFiltered)
    
    # Print graphs examples
    #showGraph(graph=graphComplete,title='Complete graph')
    #showGraph(graph=graphFiltered,title="Filtered graph")
    
    