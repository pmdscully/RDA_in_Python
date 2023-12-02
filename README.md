# RDA in Python

This is an implementation of the 'Receptor Density Algorithm' (RDA) Python 2.7. RDA is an Artificial Immune System (AIS) anomaly detection algorithm modelled upon how T-cell receptors respond to antigen, originally modelled by Owens et al in 2009 at University of York. This was part of a sub-project to investigate its applicability as an AIS anomaly classifier for our self-healing architecture.


Details of the experimentation and implementation to be added. See /images/ directory for outputs of our trial run with k=10 cross fold validation on some sampled two-stage normal data.

<img src="https://raw.githubusercontent.com/pmdscully/RDA_in_Python/master/images/out.gif" width="600" />

### Usage:
```
git clone https://github.com/pmdscully/RDA_in_Python.git
cd src/
python rda.py
```

### Adapting the Code:
- Essentially it performs binary-class classification on a given dataset.
- Feel free to modify and adapt.
- *Note: The module probably should be reworked for improved reusability, broader test evaluations, etc. (future work!)*

### Versions:

Version 0.1 implementation of RDA in python 2.7.

- Uses one numerical sensor input from a generated model.
- Includes signature classification of 1 class.
- Suitable for post-processing style classification. The anomaly matching approach is not suitable for real-time responding.
- Requires python 2.7 and libraries: scipy, numpy, matplotlib

### Tested on: Python 2.7.5+ (default, Feb 27 2014, 19:37:08)

In [1]: import matplotlib
In [2]: matplotlib.__version__
Out[2]: '1.2.1'
In [3]: import numpy
In [4]: numpy.__version__
Out[4]: '1.7.1'
In [5]: import scipy
In [6]: scipy.__version__
Out[6]: '0.12.0'


### References:
```
Hilder, J. A., Owens, N. D., Hickey, P. J., Cairns, S. N., Kilgour, D. P., Timmis, J., & Tyrrell, A. (2011). Parameter optimisation in the receptor density algorithm. In Artificial Immune Systems (pp. 226-239). Springer Berlin Heidelberg.

Owens, N. D., Greensted, A., Timmis, J., & Tyrrell, A. (2009). T cell receptor signalling inspired kernel density estimation and anomaly detection. In Artificial Immune Systems (pp. 122-135). Springer Berlin Heidelberg.
```
