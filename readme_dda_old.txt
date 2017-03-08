Old DDA README
--------------
Potvin and Shapiro dual Doppler analysis (DDA).
cd to src, edit Makefile. Go:
make
make install
or
make install PREFIX=my_binary_dir

Requires C99 and FORTRAN 77 compilers, and NetCDF.

To run, edit a copy of src/params.dda. Then go "DDA myparams.dda"
