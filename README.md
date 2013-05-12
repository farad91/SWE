SWE
===

The Shallow Water Equations teaching code.

Documentation
-------------

The documentation is available in the [Wiki](https://github.com/TUM-I5/SWE/wiki)

License
-------

SWE is release unter GPLv3 (see [gpl.txt](gpl.txt))

Group 2 additional informations and scons changes
-------------------------------------------------

Added the options: block, SWE_dimsplit, method and scenario to scons options.
block=dimsplit compiles swe_simple with our block implementatio.
SWE_dimsplit=true compiles with our replacement of swe_simple.
method=runTimestep just effects swe_simple and improves performance of our block implementation in this content. 
for more information see: scons --help 
