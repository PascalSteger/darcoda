#!/usr/bin/env python2

## \file
# class for reading fortran unformatted files

# (c) 2014 ETHZ, Pascal Steger, pascal@steger.aero

import numpy as np

class FortranFile(file):
	def __init__(self, fname, mode='rb', buf=0):
		file.__init__(self, fname, mode, buf)

	def read_fortran_record(self, dtype, endian="="):
		mytype = np.dtype(dtype).newbyteorder(endian)
		mint32 = np.dtype('i4').newbyteorder(endian)

		nbytes = np.array(0,dtype=mytype).itemsize

		n1 = np.fromfile(self, dtype=mint32, count=1)[0]/nbytes
		data = np.fromfile(self, dtype=mytype, count=n1)
		n2 = np.fromfile(self, dtype=mint32, count=1)[0]/nbytes

		if(n1 != n2):
			raise IOError("Error reading FORTRAN binary record")

		return data
    # \fn read_fortran_record(self, dtype, endian="=")
    # Read a FORTRAN record of given numpy dtype
    # @param dtype data type to be read
    # @param endian = "="


	
# \class FortranFile(file)
# Class for reading Fortran unformatted files.

