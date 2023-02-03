pyEnGNet: optimized reconstruction of gene co-expression networks using multi-GPU
=================================================================================

**Deployment & Documentation & Stats**

.. image:: https://img.shields.io/badge/pypi-v0.0.3-brightgreen
   :target: https://pypi.org/project/pyengnet/
   :alt: PyPI version


.. image:: https://readthedocs.org/projects/pyengnet/badge/?version=latest
   :target: https://pyengnet.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status


.. image:: https://img.shields.io/github/stars/aureliolfdez/pyengnet.svg
   :target: https://github.com/aureliolfdez/pyEnGNet/stargazers
   :alt: GitHub stars


.. image:: https://img.shields.io/github/forks/aureliolfdez/pyengnet.svg?color=blue
   :target: https://github.com/aureliolfdez/pyEnGNet/network
   :alt: GitHub forks


.. image:: https://img.shields.io/badge/license-GPL--3.0%20license-green
   :target: https://github.com/aureliolfdez/pyEnGNet/blob/main/LICENSE
   :alt: License


----


Abstract here


**pyengnet** is featured for:

* **Unified APIs, detailed documentation, and interactive examples** available to the community.
* **Complete coverage** for reconstruction of massive gene co-expression networks.
* **Optimized models** to generate results in the shortest possible time.
* **Optimization of a High-Performance Computing (HPC) and Big Data ecosystem**, using `cuda <https://developer.nvidia.com/cuda-zone>`_ and `multiprocess <https://github.com/uqfoundation/multiprocess>`_.

**API Demo**\ :

.. code-block:: python


      import os
      from pyengnet.File import File
      from pyengnet.Engnet import Engnet

      if __name__ == "__main__":
         
         # Load dataset
         dataset = File.load(path=os.getcwd()+"/datasets/Spellman.csv", separator=",", nmi_th=0.6, spearman_th=0.7, kendall_th=0.7, readded_th=0.7, hub_th = 3)    

         # Run pyEnGNet on CPUs
         graphFiltered, infoGraphFiltered, graphComplete, infoGraphComplete = Engnet.process(dataset, saveComplete = True)

         # Run pyEnGNet on GPU devices
         # graphFiltered, infoGraphFiltered, graphComplete, infoGraphComplete = Engnet.process(dataset, saveComplete = True, numGpus = 2, computeCapability = 61)
         
         # Save gene co-expression networks and additional information
         File.saveFile(path='/home/principalpc/Escritorio/graphComplete.csv',graph=infoGraphComplete) # Full network
         File.saveFile(path='/home/principalpc/Escritorio/graphFiltered.csv',graph=infoGraphFiltered) # Filtered network
         
         # Print gene co-expression networks
         File.showGraph(graph=graphComplete,title='Complete graph') # Full network
         File.showGraph(graph=graphFiltered,title="Filtered graph") # Filtered network


**Citing pyEnGNet**\ :

`pyEnGNet paper <#>`_ is published in
`(under review) <#>`_.
If you use pyEnGNet in a scientific publication, we would appreciate citations to the following paper::

   Under review

or::

    Under review


**Key Links and Resources**\ :

* `View the latest codes on Github <https://github.com/aureliolfdez/pyEnGNet>`_
* `View the documentation & API <https://pyengnet.readthedocs.io/>`_
* `View all examples <https://github.com/aureliolfdez/pyEnGNet/tree/main/tests/test_integration>`_

----

Installation
============

It is recommended to use **pip** for installation. Please make sure
**the latest version** is installed, as pyengnet is updated frequently:

.. code-block:: bash

   pip install pyengnet            # normal install
   pip install --upgrade pyengnet  # or update if needed
   pip install --pre pyengnet      # or include pre-release version for new features

Alternatively, you could clone and run setup.py file:

.. code-block:: bash

   git clone https://github.com/aureliolfdez/pyEnGNet.git
   pip install .

**Required Dependencies**\ :

* Python>=3.10
* numpy>=1.24.0
* tqdm>=4.64.0
* multiprocess>=0.70.14
* pandas>=1.5.3
* matplotlib>=3.6.3
* networkx>=3.0
* scipy>=1.10.0

API Reference
=============

I/O Management
^^^^^^^^^^^^^^^^^^^^^^

* **pyengnet.File**\: Class used to manage file I/O operations and data visualization.
* **pyengnet.File.load()**\: Load dataset from a txt or csv file.
* **pyengnet.File.saveFile()**\: Save network to file (can be used to store full and/or pruned networks)
* **pyengnet.File.showGraph()**\: Display a specific network

----


Ensemble
^^^^^^^^^^^^^^^^^^^
* **pyengnet.Engnet**\: Class in charge of controlling the execution of the EnGNet algorithm.
* **pyengnet.Engnet.process()**\: Function that runs the EngNet algorithm. Depending on the parameters of this function, the algorithm will be executed in parallel with CPU processors or GPU devices.
* **pyengnet.Kendall**\: Kendall measurement is coded in a parallel ecosystem with CPUs.
* **pyengnet.NMI**\: NMI measurement is coded in a parallel ecosystem with CPUs.
* **pyengnet.Spearman**\: Spearman measurement is coded in a parallel ecosystem with CPUs.
* **pyengnet.src.correlations**\: Execution of Kendall, NMI, and Spearman measures under a parallel multi-GPU ecosystem (CUDA). In addition, it detects those pairs of genes that exceed the threshold for major voting.

Examples by Tasks
=================


**All implemented modes** are associated with examples, check
`"pyEnGNet examples" <https://github.com/aureliolfdez/pyEnGNet/tree/main/tests/test_integration>`_
for more information.


----


Run on CPU
^^^^^^^^^^^^^^^^^^^^^^^^^^^

`"tests/test_integration/test_cpu.py" <https://github.com/aureliolfdez/pyEnGNet/tree/main/tests/test_integration/test_cpu.py>`_
demonstrates the basic API for the generation of co-expression gene networks using CPUs.

#. Load gene co-expression dataset from input file

   .. code-block:: python

      from pyengnet.File import File

      dataset = File.load(path=os.getcwd()+"/datasets/Spellman.csv", separator=",", nmi_th=0.6, spearman_th=0.7, kendall_th=0.7, readded_th=0.7, hub_th = 3)    


#. Run pyEnGNet based on CPUs.

   .. code-block:: python

      from pyengnet.Engnet import Engnet

      graphFiltered, infoGraphFiltered, graphComplete, infoGraphComplete = Engnet.process(dataset, saveComplete = True)

#. Save gene co-expression networks output (optional)

   .. code-block:: python
      
      from pyengnet.File import File
      
      File.saveFile(path='/home/user/Desktop/graphComplete.csv',graph=infoGraphComplete)
      File.saveFile(path='/home/user/Desktop/graphFiltered.csv',graph=infoGraphFiltered)

#. Print gene co-expression networks output  (optional)

   .. code-block:: python
      
      from pyengnet.File import File
      
      File.showGraph(graph=graphComplete,title='Complete graph')
      File.showGraph(graph=graphFiltered,title="Filtered graph")


Run on GPU devices
^^^^^^^^^^^^^^^^^^^^^^^^^^^

`"tests/test_integration/test_gpu.py" <https://github.com/aureliolfdez/pyEnGNet/tree/main/tests/test_integration/test_gpu.py>`_
demonstrates the basic API for the generation of co-expression gene networks using GPU devices.

#. Load gene co-expression dataset from input file

   .. code-block:: python

      from pyengnet.File import File

      dataset = File.load(path=os.getcwd()+"/datasets/Spellman.csv", separator=",", nmi_th=0.6, spearman_th=0.7, kendall_th=0.7, readded_th=0.7, hub_th = 3)    

#. Run pyEnGNet based on CPUs.

   .. code-block:: python

      from pyengnet.Engnet import Engnet

      graphFiltered, infoGraphFiltered, graphComplete, infoGraphComplete = Engnet.process(dataset, saveComplete = True, numGpus = 2, computeCapability = 61)

#. Save gene co-expression networks output (optional)

   .. code-block:: python
      
      from pyengnet.File import File
      
      File.saveFile(path='/home/user/Desktop/graphComplete.csv',graph=infoGraphComplete)
      File.saveFile(path='/home/user/Desktop/graphFiltered.csv',graph=infoGraphFiltered)

#. Print gene co-expression networks output  (optional)

   .. code-block:: python
      
      from pyengnet.File import File
      
      File.showGraph(graph=graphComplete,title='Complete graph')
      File.showGraph(graph=graphFiltered,title="Filtered graph")
