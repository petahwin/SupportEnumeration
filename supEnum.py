#!/usr/bin/python

try:
    import numpy
except:
    print "numpy is not installed"

print "NumPy is installed: verion " + numpy.__version__

numpy.test()

