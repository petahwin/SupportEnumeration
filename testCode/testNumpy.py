#!/usr/bin/python
import numpy

a = numpy.array([[4,3,1], [3,7,1], [1,1,0]])
b = numpy.array([0,0,1])

print numpy.linalg.lstsq(a,b)[0]



