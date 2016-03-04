# physoce 
Python tools for physical oceanography at Moss Landing Marine Labs

## Installation

To install the package for the first time, first run the following in the 
terminal/command prompt:
```
git init .
git remote add origin https://gitlab.com/physoce/physoce-py.git
```

## Importing modules

Python code to import the main package:
```python
import physoce
```

Python code to import sub-packages:
```python
import physoce.obs
```

CONTENTS

io.py module         	- input/output 

stats.py module         - statistics

surfacewaves.py module  - surface gravity wave calculations

tseries.py module  		- time series 

obs/ sub-package        - modules for specific types of oceanographic data

obs/lobo.py module      - load MBARI LOBO mooring data