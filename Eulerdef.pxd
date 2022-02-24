#cython:language level = 3
cdef extern from "Euler1D.cpp":
    pass

from libcpp.string cimport string

cdef extern from "Euler1D.h":
    cdef cppclass Euler1D:
        Euler1D(string bctype, int n) except + 
        void initialize() except + 
        void set_boundary() except + 
        void reconstruction_0() except + 
        void reconstruction_TVD() except +
        void avg_roe() except + 
        void cal_roe_flux() except +
        void writefile(string filename) except +
        void timeadvancement(double dt, int oldset) except +