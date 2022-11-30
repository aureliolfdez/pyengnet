from base import *
from algorithms import *

if __name__ == "__main__":
    filePath = os.getcwd()+"/src/pyEnGNet/datasets/Spellman.csv"
    dataset = load(path=filePath, separator=",", nmi_th=0.5, spearman_th=0.5, kendall_th=0.5, readded_th=0.5)    
    graphFiltered, infoGraphFiltered, graphComplete, infoGraphComplete = Engnet.process(dataset, saveComplete = True)
    
    # Save files info
    saveFile(path='/Users/aurelio/Desktop/graphComplete.csv',graph=infoGraphComplete)
    saveFile(path='/Users/aurelio/Desktop/graphFiltered.csv',graph=infoGraphFiltered)
    
    # Print graphs examples
    showGraph(graph=graphComplete,title='Complete graph')
    showGraph(graph=graphFiltered,title="Filtered graph")
    
    