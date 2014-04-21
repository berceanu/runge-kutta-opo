from numpy import linspace, empty
from numpy cimport ndarray as ar

cdef extern from "pyrk_adaptive.h":
    void c_setg(int* Nx, int* Ny, double* Lx, double* Ly, double* kinetic)

# void c_gfunc(double* a, int* n, int* m, double* a, double* b, double* c)


def py_setg(int nx, int ny, double lx, double ly):
    cdef ar[double,ndim=2] kinetic = empty((nx, ny), order='F')
    c_setg(&nx, &ny, &lx, &ly, <double*> kinetic.data)
    return kinetic
