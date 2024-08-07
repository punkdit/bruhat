#!/bin/sh

cython3 -X language_level=3 _element.pyx 

gcc -DNDEBUG -g -O3 -Wstrict-prototypes -fPIC -I/usr/include/python3.10 \
    -c _element.c -o _element.o

gcc -shared _element.o -o _element.so


./test_element.py

