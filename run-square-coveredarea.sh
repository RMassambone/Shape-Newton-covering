gfortran -c -O3 square-coveredarea.f90
gfortran -c -O3 vorintpols.f90
gfortran -c -O3 geometry.f90
gfortran -c -O3 geompack2.f90

gfortran square-coveredarea.o vorintpols.o geometry.o geompack2.o \
         -o square-coveredarea

./square-coveredarea > square-coveredarea.out
