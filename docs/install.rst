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