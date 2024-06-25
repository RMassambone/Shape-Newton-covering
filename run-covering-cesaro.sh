export ALGENCAN=$HOME/algencan-4.0.0

gfortran -c -O3 modcesaro.f90
gfortran -c -O3 vorintpols.f90
gfortran -c -O3 drand.f90
gfortran -c -O3 geometry.f90
gfortran -c -O3 geompack2.f90
gfortran -c -O3 -I$ALGENCAN/sources/algencan/inc/ covering-cesaro.f90

gfortran covering-cesaro.o modcesaro.o vorintpols.o drand.o geometry.o geompack2.o \
         -L$ALGENCAN/sources/algencan/lib -L$ALGENCAN/sources/hsl/lib \
	 -L$ALGENCAN/sources/blas/lib -lalgencan -lhsl -lblas -o covering-cesaro

for nballs in `seq 1 1 100`
do
  export bestrad=1000000.0d0
  for itrial in `seq 1 1 10000`
  do
    echo $nballs $itrial $bestrad | timeout --preserve-status 30 ./covering-cesaro
    export statdtris=$?
    if [ -e tabline.txt ] ; then
      sed -e 's/$/  S/' -i tabline.txt
      cat tabline.txt >> table-covering-cesaro.txt ;
      bestrad=`cat bestrad.txt`
    else
      if [ $statdtris -eq 0 ] && [ -e alsolver-interrupted-tabline.txt ]; then
        sed -e 's/$/  S/' -i alsolver-interrupted-tabline.txt
        cat alsolver-interrupted-tabline.txt >> table-covering-cesaro.txt ;
      elif [ $statdtris != 0 ] && [ -e alsolver-interrupted-tabline.txt ]; then
        sed -e 's/$/  F/' -i alsolver-interrupted-tabline.txt
        cat alsolver-interrupted-tabline.txt >> table-covering-cesaro.txt ;
      else
        printf " no information\n" >> table-covering-cesaro.txt ;
      fi
    fi
    rm -f tabline.txt alsolver-interrupted-tabline.txt \
          solver-interrupted-tabline.txt newtkkt-interrupted-tabline.txt
  done
done