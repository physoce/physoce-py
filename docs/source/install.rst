Installation
============

**Method 1: Install using pip**

To install using pip, run the following in the terminal/command prompt:

.. code-block:: none

   pip install git+https://github.com/physoce/physoce-py

**Method 2: Install using git**

To install the package for the first time, first create an empty directory where you plan to store the code, then change into that directory.

Run the following in the terminal/command prompt to initialize an empty git repository:

.. code-block:: none

   git init .

Add the remote Github url.

.. code-block:: none

   git remote add origin https://github.com/physoce/physoce-py.git

Get code and install.

.. code-block:: none

   git pull origin master
   python setup.py install

After the initial installation, you can update to the latest version on Github by changing to the directory where the code is located running the following in the terminal/command prompt:

.. code-block:: none

   git pull origin master
   python setup.py install
