# physoce 
Python tools for physical oceanography at Moss Landing Marine Labs

## Installation

To install the package for the first time, first create an empty directory where
you plan to store the code, then go to that directory and run the following in 
the terminal/command prompt:
```
git init .
git remote add origin https://gitlab.com/physoce/physoce-py.git
git pull origin master
python setup.py install
```

## Updating

After the initial installation, you can update to the latest version on Gitlab
by changing to the directory where the code is located running the following in
the terminal/command prompt:
```
git pull origin master
python setup.py install
```

## Importing modules

After installing, the package can be accessed in Python. Here is code to import 
the main package and see the available modules and sub-packages:
```python
import(phyoce)
help(physoce)
```

Python code to import a module from the main physoce package and see the 
functions available in that module:
```python
from physoce import graph
help(graph)
```

Python code to import module from the physoce.obs sub-package and see the 
functions available in that module:
```python
from physoce.obs import nerr
help(nerr)
```

## Contents

* graph.py module         - graphing and plotting
* io.py module         	- input/output 
* stats.py module         - statistics
* surfacewaves.py module  - surface gravity wave calculations
* tseries.py module  		- time series 
* obs/ sub-package        - modules for specific types of oceanographic data
* obs/lobo.py module      - load MBARI LOBO mooring data
* obs/nerr.py module      - load National Estuarine Research Reserve data

## Adding new code to the package

See the CONTRIBUTING.md file