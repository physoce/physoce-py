# physoce

Python tools for physical oceanography

## Contents

#### physoce package

* io.py                   - input/output 
* graph.py                - graphing and plotting
* oceans.py               - oceanography-specific
* stats.py                - statistics
* tseries.py              - time series 
* util.py                 - general-purpose

#### physoce.obs sub-package 
Modules for specific types of publically-available oceanographic data

* obs/lobo.py module      - MBARI LOBO mooring data
* obs/nerr.py module      - National Estuarine Research Reserve data
* obs/elkhorn.py module	  - Elkhorn Slough Reserve GIS data
* obs/mlml.py module	  - Moss Landing Marine Labs public data

## Installation

#### Inititialize git repository

To install the package for the first time, first create an empty directory where
you plan to store the code, then change into that directory.

Run the following in the terminal/command prompt to initialize an empty git repository:

```
git init .
```

#### Add remote url

Use one of the following commands to add a remote repository url, depending on whether you
are pulling from the MLML development version on Gitlab,

```
git remote add origin https://gitlab.com/physoce/physoce-py.git
```

or the public version on Github.

```
git remote add origin https://github.com/tompc35/physoce.git
```

#### Get code and install

```
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

After installing, the package can be accessed in Python. Here is Python code to import 
the main package and see the available modules and sub-packages:

```python
import phyoce
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
<<<<<<< HEAD
```
=======
```

## Contents

* io.py module         	- input/output 
* graph.py module         - graphing and plotting
* oceans.py module  - oceanography-specific
* stats.py module         - statistics
* tseries.py module  		- time series 
* util.py - general-purpose
* obs/ sub-package        - modules for specific types of oceanographic data
* obs/lobo.py module      - MBARI LOBO mooring data
* obs/nerr.py module      - National Estuarine Research Reserve data
* obs/mlml.py module	  - Moss Landing Marine Labs public data

<<<<<<< HEAD
## Adding new code to the package

See the CONTRIBUTING.md file
>>>>>>> Adding a unit test for tseries.py and updated the gitignore to ignore some pycharm files
=======
## Unit Testing
*** Reformat for Push *** <br>
I am going to add a unit test for the tseries.py modules. I have spent enough hacking together a test in my processing code,
that I think it is worth thinking about.
[Here is where I'm starting](https://jeffknupp.com/blog/2013/12/09/improve-your-python-understanding-unit-testing/)



## Adding new code to the package

See the CONTRIBUTING.md file
>>>>>>> unit_test
>>>>>>> 4dc387ba900be3f4190b99ed733fb04159ba205b
