# starsdb
A minimal app for converting NAIF/SPICE's mast star cataloges (bdb-files, binary E-kernels) to a more easy to use internal format.

Work in progress..

starsdb now has the following minimal features:
- verbose mode
- Convert E-kernel master star catalogues to IFS-format (see file table description in IFS.txt)
- List star properties
- List stars brighter than given apparent magnitude

How to get the star catalogues from NIAF:
* navigate or use wget from 
```
http://naif.jpl.nasa.gov/pub/naif/generic_kernels/stars/
```
* The Hipparcos catalogue has about 117.000 stars and should be sufficient for most uses
* unzip the transfer file (xdb):
```
gzip -d hipparcos.xdb.Z
```
* Convert to binary E-kernel format with Spice's "tobin":
```
tobin hipparcos.xdb hipparcos.bdb
```

The bdb files can then be read by this app.

To build:
```
mkdir build
cd build
cmake ..
make
```

Executable will be placed in bin/

For info on usage:
```
starsdb -h
```

