from base import *
from algorithms import *

if __name__ == "__main__":
    filePath = os.getcwd()+"/src/pyEnGNet/datasets/Spellman_v2.csv"
    dataset = load(path=filePath, separator=",", nmi_th=0.5, spearman_th=0.5, kendall_th=0.5, readded_th=0.5)    
    #graphFiltered, infoGraphFiltered, graphComplete, infoGraphComplete = Engnet.process(dataset, saveComplete = True, numGpus = 1, computeCapability = 61)
    graphFiltered, infoGraphFiltered, graphComplete, infoGraphComplete = Engnet.process(dataset, saveComplete = True)
    
    # Save files info MacOSX
    #saveFile(path='/Users/aurelio/Desktop/graphComplete.csv',graph=infoGraphComplete)
    #saveFile(path='/Users/aurelio/Desktop/graphFiltered.csv',graph=infoGraphFiltered)
    
    # Save files info Ubuntu
    saveFile(path='/home/principalpc/Escritorio/graphComplete.csv',graph=infoGraphComplete)
    saveFile(path='/home/principalpc/Escritorio/graphFiltered.csv',graph=infoGraphFiltered)
    
    # Print graphs examples
    #showGraph(graph=graphComplete,title='Complete graph')
    #showGraph(graph=graphFiltered,title="Filtered graph")
    
    