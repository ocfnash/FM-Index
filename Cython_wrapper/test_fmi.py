#!/usr/bin/python

from FMIndex import PyFMIndex as fmi

print("---------------: Creating fmi object\n\'hello-there\'")
x = fmi(b'hello-there')

print("---------------: find \'hello\'")
find1 = x.findn(b'hello')
print(find1)

print("---------------: find \'l\'")
find2 = x.findn(b'l')
print(find2)