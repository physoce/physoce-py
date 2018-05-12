Installation
============

**Initialize git repository**

To install the package for the first time, first create an empty directory where you plan to store the code, then change into that directory.

Run the following in the terminal/command prompt to initialize an empty git repository:

    git init .

**Add remote url**

    git remote add origin https://github.com/physoce/physoce-py.git

**Get code and install**

    git pull origin master
    python setup.py install

Updating
========

After the initial installation, you can update to the latest version on Gitlab by changing to the directory where the code is located running the following in the terminal/command prompt:

    git pull origin master
    python setup.py install
