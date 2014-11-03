FM-Index
========

This project contains a C++ implementation of the FM Index data structure together with the necessary Cython code to compile it as a Python module.

The FM Index data structure stores text (or even just a blob of binary data) in such a way that a substring query can be answered in optimal asymptotic time, i.e., O(L) time for a substring of length L. Note the perhaps surprising fact that this time complexity does not depend on the length of the text.

The FM Index data structure is built by calculating the Burrows-Wheeler transform of the text which takes some time and is not easy to update if the text changes. The data structure is thus best suited to the problem of repeatedly searching the same text for many substrings.

An example application and a detailed description of the data structure can be found at my blog:
[link](http://ocfnash.wordpress.com/2014/01/03/dna-of-a-password-disaster/)

A simple shell script is provided to create the Python module which should work on any system on which Cython is installed:
Cython_wrapper/make.sh

GoogleTest (gtest) must be installed in order to run the unit tests:
UnitTests/unit_tests.cpp

I have also included various papers which introduce the key ideas of the FM Index data structure in the docs directory.

Finally note that I have used Yuta Mori's OpenBWT code to compute the Burrows-Wheeler transform.

## Building and using

% pip install cython  
% cd Cython_wrapper  
% ./make.sh  
% python  
`>>>` from FMIndex import PyFMIndex as fmi  
`>>>` x = fmi('hello-there')  
`>>>` x.findn('hello')  
1  
`>>>` x.findn('l')  
2
