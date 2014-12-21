#!/usr/bin/env ipython3

## start program with
# LD_PRELOAD=/usr/lib64/libmpi.so mpirun -n 1 ./reprun.py

import numpy as np
import ctypes as ct



#_libfunctions = np.ctypeslib.load_library('libfunctions', '.')
#_libfunctions.hello.argtypes = [ct.c_int]
#_libfunctions.hello.restype  =  ct.c_void_p
#_libfunctions.dprod.argtypes = [np.ctypeslib.ndpointer(dtype=np.float), ct.c_int]
#_libfunctions.dprod.restype  = ct.c_double
#_libfunctions.dcumsum.argtypes = [np.ctypeslib.ndpointer(dtype=np.float), np.ctypeslib.ndpointer(dtype=np.float), ct.c_int]
#_libfunctions.dcumsum.restype  = ct.c_void_p


#it = ct.CDLL('./it.so')
itlib = np.ctypeslib.load_library('it', '.')
itlib.__itimes_MOD_test.argtypes = [ct.c_int]
itlib.__itimes_MOD_test.restype = ct.c_float

nm = 200
itlib.__itimes_MOD_test(ct.byref(ct.c_int(nm)))
