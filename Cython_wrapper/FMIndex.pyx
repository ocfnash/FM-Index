# distutils: language = c++
# distutils: include_dirs = ../FM-Index ../openbwt-v1.5
# distutils: sources = ../FM-Index/FMIndex.cpp ../FM-Index/WaveletTree.cpp ../FM-Index/BitVector.cpp ../openbwt-v1.5/BWT.c

from libcpp.list cimport list
from libcpp.string cimport string

cdef extern from "FMIndex.h":
    cdef cppclass FMIndex:
        FMIndex(string) except +
        int findn(string)
        list[string] find_lines(string)
        void serialize_to_file(string)
        size_t size()
    #cdef FMIndex * new_from_serialized_file "FMIndex::new_from_serialized_file"(string)

cdef extern from "FMIndex.h" namespace "FMIndex": # static member function hack
    FMIndex * new_from_serialized_file(string)

cdef class PyFMIndex:
    cdef FMIndex * thisptr
    def __cinit__(self, s):
        self.thisptr = new FMIndex(s)
    def __dealloc__(self):
        del self.thisptr
    def findn(self, pattern):
        return self.thisptr.findn(pattern)
    def find_lines(self, pattern):
        return self.thisptr.find_lines(pattern)
    def new_from_serialized_file(self, filename):
        del self.thisptr
        self.thisptr = new_from_serialized_file(filename)
    def serialize_to_file(self, filename):
        self.thisptr.serialize_to_file(filename)
    def size(self):
        return self.thisptr.size()
