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