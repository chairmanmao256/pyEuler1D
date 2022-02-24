#distutils:language=c++
#cython:language_level=3

from Eulerdef cimport Euler1D as _Euler
from libcpp.string cimport string

cdef class Euler1D:

    cdef _Euler* pEuler

    def __cinit__(self):
        self.pEuler = NULL

    def __dealloc__(self):
        if self.pEuler != NULL:
            del self.pEuler

    def __init__(self, string & bc, int N):
        self.pEuler = new _Euler(bc, N)

    def initialize(self):
        self.pEuler.initialize()
    
    def set_boundary(self):
        self.pEuler.set_boundary()

    def reconstruction_0(self):
        self.pEuler.reconstruction_0()

    def reconstruction_TVD(self):
        self.pEuler.reconstruction_TVD()

    def avg_roe(self):
        self.pEuler.avg_roe()

    def cal_roe_flux(self):
        self.pEuler.cal_roe_flux()

    def writefile(self, string & filename):
        self.pEuler.writefile(filename)

    def timeadvancement(self,double dt, int oldset):
        self.pEuler.timeadvancement(dt, oldset)