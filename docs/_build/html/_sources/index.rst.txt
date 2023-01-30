.. pyEnGNet documentation master file, created by
   sphinx-quickstart on Fri Jan 27 00:37:03 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyEnGNet's documentation!
====================================

**Deployment & Documentation & Stats**

.. image:: https://img.shields.io/badge/pypi-v0.0.0-brightgreen
   :target: https://pypi.org/project/pyengnet/
   :alt: PyPI version


.. image:: https://readthedocs.org/projects/pyengnet/badge/?version=latest
   :target: https://pyengnet.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status


.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/aureliolfdez/pyengnet/master
   :alt: Binder


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

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   install
   examples

.. toctree::
   :maxdepth: 2
   :caption: Documentation

   api

.. toctree::
   :maxdepth: 2
   :caption: Additional information

   about
   faq
   release


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
