API Reference
=============

I/O Management
^^^^^^^^^^^^^^^^^^^^^^

* :class:`pyengnet.File`: Class used to manage file I/O operations and data visualization.
* :func:`pyengnet.File.load()`: Load dataset from a txt or csv file.
* :func:`pyengnet.File.saveFile()`: Save network to file (can be used to store full and/or pruned networks)
* :func:`pyengnet.File.showGraph()`: Display a specific network

----


Ensemble
^^^^^^^^^^^^^^^^^^^
* :class:`pyengnet.Engnet`: Class in charge of controlling the execution of the EnGNet algorithm.
* :func:`pyengnet.Engnet.process()`: Function that runs the EngNet algorithm. Depending on the parameters of this function, the algorithm will be executed in parallel with CPU processors or GPU devices.
* :class:`pyengnet.Kendall`: Kendall measurement is coded in a parallel ecosystem with CPUs.
* :class:`pyengnet.NMI`: NMI measurement is coded in a parallel ecosystem with CPUs.
* :class:`pyengnet.Spearman`: Spearman measurement is coded in a parallel ecosystem with CPUs.
* :class:`pyengnet.src.correlations`: Execution of Kendall, NMI, and Spearman measures under a parallel multi-GPU ecosystem (CUDA). In addition, it detects those pairs of genes that exceed the threshold for major voting.



----

All Models
^^^^^^^^^^

.. toctree::
   :maxdepth: 4

   pyengnet
