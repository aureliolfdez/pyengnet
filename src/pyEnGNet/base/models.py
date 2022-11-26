import multiprocessing as mp

class Dataset:
    
    """
    This class is a Dataset model

    Parameters
    ----------
    data : Array
        Original dataset
        
    spearman_threshold : float
        Spearman threshold
    
    kendall_threshold : float
        Kendall threshold
    
    nmi_threshold : float
        NMI threshold
    
    readded_edges_threshold : float
        If hubs exist in a network, edges whose weight exceeds the threshold indicated by this parameter shall be added.
    
    ncores : int
        Number of cores used        
    """
    
    def __init__(self, data=None, nmi_th=0.7, spearman_th=0.7, kendall_th=0.7, readded_th=0.7, cores=int(mp.cpu_count() / 2)):
       
        self._data = data
        
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
        self._ncores = cores
    
    @property
    def data(self):
        return self._data
    
    @data.setter
    def data(self, data):
        self._data = data
    
    @property
    def row_size(self):
        return self._row_size
    
    @row_size.setter
    def row_size(self, row_size):
        self._row_size = row_size
    
    @property
    def column_size(self):
        return self._column_size
    
    @column_size.setter
    def column_size(self, column_size):
        self._column_size = column_size
    
    @property
    def spearman_threshold(self):
        return self._spearman_threshold
    
    @spearman_threshold.setter
    def spearman_threshold(self, spearman_threshold):
        self._spearman_threshold = spearman_threshold
    
    @property
    def kendall_threshold(self):
        return self._kendall_threshold
    
    @kendall_threshold.setter
    def kendall_threshold(self, kendall_threshold):
        self._kendall_threshold = kendall_threshold
    
    @property
    def nmi_threshold(self):
        return self._nmi_threshold
    
    @nmi_threshold.setter
    def nmi_threshold(self, nmi_threshold):
        self._nmi_threshold = nmi_threshold
      
    @property
    def readded_edges_threshold(self):
        return self._readded_edges_threshold
    
    @readded_edges_threshold.setter
    def readded_edges_threshold(self, readded_edges_threshold):
        self._readded_edges_threshold = readded_edges_threshold
    
    @property
    def ncores(self):
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
    
    
    