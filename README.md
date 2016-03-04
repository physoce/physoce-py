physoce 
-------
Python tools for physical oceanography at Moss Landing Marine Labs

INSTALLATION

To install the package, run the following in the terminal/command prompt:
python setup.py install

IMPORTING

Python code to import the main package:
```import physoce```

Python code to import sub-packages:
```import physoce.obs```

CONTENTS

io.py module         	- input/output 

stats.py module         - statistics

surfacewaves.py module  - surface gravity wave calculations

tseries.py module  		- time series 

obs/ sub-package        - modules for specific types of oceanographic data

obs/lobo.py module      - load MBARI LOBO mooring data