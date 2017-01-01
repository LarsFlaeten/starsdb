# starsdb
A minimal app for converting NAIF/SPICE's mast star cataloges (binary E-kernels) to a more easy to use internal format.

Work in progress..

starsdb now has the following minimal features:
- verbose mode
- Convert E-kernel maser star catalogues to IFS-format (see file table description in IFS.txt)
- List star properties
- List stars brighter than given apparent magnitude


To build:
```
mkdir build
cd build
cmake ..
make
```

Executable will be placed in bin/

For info on usage:
'''
starsdb -h
'''
