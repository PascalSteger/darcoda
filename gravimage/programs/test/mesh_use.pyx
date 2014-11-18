from numpy cimport ndarray
from numpy import empty, double

cdef extern:
    void c_mesh_exp(double *r_min, double *r_max, int *N, double *mesh)

def mesh_exp(double r_min, double r_max, int N):
    cdef ndarray[double, mode="c"] mesh = empty(N, dtype=double)
    c_mesh_exp(&r_min, &r_max,  &N, &mesh[0])
    return mesh

mm = mesh_exp(0., 1., 6)
print(mm)
