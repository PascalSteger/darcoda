from numpy cimport ndarray
from numpy import empty, double


cdef extern:
    void c_mesh_exp(double *r_min, double *r_max, int *N, double *mesh, bint *pp, char msg[]);

cdef extern:
    void pyx_mdef_4grav_3samplefunc();

def samplefunc():
    pyx_mdef_4grav_3samplefunc()

def mesh_exp(double r_min, double r_max, int N, bint pp, bytes msg):
    cdef ndarray[double, mode="c"] mesh = empty(N, dtype=double)
    c_mesh_exp(&r_min, &r_max,  &N, &mesh[0], &pp, msg)
    return mesh


my_str = "hello world"
print('the string message in python is: ', my_str)

byt = str.encode(my_str)
print('the bytes message in python is: ', byt)
print('types: ', type(byt))


mm = mesh_exp(0., 1., 6, True, byt)
print('.. and in python, we can access it as .. ')
print(mm)
