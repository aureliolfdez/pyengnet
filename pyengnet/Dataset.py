import multiprocessing as mp

class Dataset:
    
    """
    This class represents a Dataset model
    """
    
    def __init__(self, filePath = None, data=None, gene = None, nmi_th=0.7, spearman_th=0.7, kendall_th=0.7, readded_th=0.7, hub_th=2, cores=int(mp.cpu_count() / 2)):
       
        """
        Initialises a Dataset object with the data obtained by parameter.
        
        Parameters
        ----------
        
        filePath : str
            The path where the dataset file is stored.
            
        data : numpy.Array
            The dataset load in memory. (This dataset is obtained by pyengnet.File.load())
        
        gene: numpy.Array
            An array storing the names of the genes in the dataset.
        
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
        
        """
       
        self._data = data
        self._gene = gene
        self._filePath = filePath
        
        if data is not None:
            self._row_size = len(data)
            self._column_size = len(data[0])
        else:
            self._row_size = -1
            self._column_size = -1

        self._spearman_threshold = spearman_th
        self._kendall_threshold = kendall_th
        self._nmi_threshold = nmi_th

        self._readded_edges_threshold = readded_th
        self._hub_threshold = hub_th
        self._ncores = cores
    
    @property
    def filePath(self):
        """
            The path where the dataset file is stored.
        """
        return self._filePath
    
    @filePath.setter
    def filePath(self, filePath):
        self._filePath = filePath
        
    @property
    def data(self):
        """
            The dataset stored in memory.
        """
        return self._data
    
    @data.setter
    def data(self, data):
        self._data = data
    
    @property
    def gene(self):
        """
            Gene names list of the dataset stored in a numpy.array
        """
        return self._gene
    
    @gene.setter
    def gene(self, gene):
        self._gene = gene
    
    @property
    def row_size(self):
        """
            Number of rows of the dataset
        """
        return self._row_size
    
    @row_size.setter
    def row_size(self, row_size):
        self._row_size = row_size
    
    @property
    def column_size(self):
        """
            Number of columns of the dataset
        """
        return self._column_size
    
    @column_size.setter
    def column_size(self, column_size):
        self._column_size = column_size
    
    @property
    def spearman_threshold(self):
        """
            Threshold for the Spearman classifier (Used in phase 1 of the EnGNet algorithm).
        """
        return self._spearman_threshold
    
    @spearman_threshold.setter
    def spearman_threshold(self, spearman_threshold):
        self._spearman_threshold = spearman_threshold
    
    @property
    def kendall_threshold(self):
        """
            Threshold for the Kendall classifier (Used in phase 1 of the EnGNet algorithm).
        """
        return self._kendall_threshold
    
    @kendall_threshold.setter
    def kendall_threshold(self, kendall_threshold):
        self._kendall_threshold = kendall_threshold
    
    @property
    def nmi_threshold(self):
        """
            Threshold for the NMI classifier (Used in phase 1 of the EnGNet algorithm).
        """
        return self._nmi_threshold
    
    @nmi_threshold.setter
    def nmi_threshold(self, nmi_threshold):
        self._nmi_threshold = nmi_threshold
      
    @property
    def readded_edges_threshold(self):
        """
            Threshold to determine if the edge would be return into the network after the pruning step
        """
        return self._readded_edges_threshold
    
    @readded_edges_threshold.setter
    def readded_edges_threshold(self, readded_edges_threshold):
        self._readded_edges_threshold = readded_edges_threshold
    
    @property
    def hub_threshold(self):
        """
            Threshold to determine if the node studie is a hub. Set this threshold to -1 to run the algorithm with standard selection.
        """
        return self._hub_threshold
    
    @hub_threshold.setter
    def hub_threshold(self, hub_threshold):
        self._hub_threshold = hub_threshold
    
    @property
    def ncores(self):
        """
            Number of CPU cores used for parallelisation.
        """
        return self._ncores
    
    @ncores.setter
    def ncores(self, ncores):
        self._ncores = ncores
    
    def __str__(self):
        s = f"PyEnGNet with shape ({self.row_size} rows, {self.column_size} columns)"
        s += f"\n###########################################"
        s += f"\n\n\tSpearman threshold selected: {self.spearman_threshold}"
        s += f"\n\tKendall Tau threshold selected: {self.kendall_threshold}"
        s += f"\n\tNMI threshold selected: {self.nmi_threshold}"
        s += f"\n\tWeight Threshold selected for edge readition: {self.readded_edges_threshold}"
        s += f"\n\tUsing {self.ncores} of your cores to process results\n\n"
        return s
    
    
    