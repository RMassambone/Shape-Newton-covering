# Shape-Newton-covering
Source codes related to the paper: E. G. Birgin, A. Laurain, R. Massambone, and A. G. Santana, A Shape-Newton approach to the problem of covering with identical balls, published in SIAM Journal on Scientific Computing.

[SIAM](https://epubs.siam.org/doi/10.1137/21M1426067)

[PDF](https://www.ime.usp.br/~egbirgin/publications/coveringsecond.pdf)


All codes needed to reproduce the numerical experiments in the paper
are included. The third party codes Algencan 4.0.0 and GEOMPACK are
required.

ALGENCAN
========

To install Algencan 4.0.0, visit:

https://www.ime.usp.br/~egbirgin/sources/bmcomper/

(Note that Algencan 4.0.0 requires routines from BLAS and HSL as
well.)

It is assumed that Algencan 4.0.0 was successfully installed, that the
environment variable ALGENCAN points to the folder were Algencan 4.0.0
was installed, and that within folders
$ALGENCAN/sources/algencan/lib/, $ALGENCAN/sources/blas/lib/, and
$ALGENCAN/sources/hsl/lib/ there are files named libalgencan.a,
libblas.a, and libhsl.a, respectively.

GEOMETRY AND GEOMPACK
=====================

The two third party packages GEOMETRY and GEOMPACK are also
required. They correspond to files geometry.f90 and geompack2.f90,
taken from

https://people.sc.fsu.edu/~jburkardt/f_src/geometry/geometry.html

and

https://people.math.sc.edu/Burkardt/f_src/geompack2/geompack2.html

that are included in the current distribution. They are compiled
together with the other codes and no specific installation process is
required.

COVERING CODE
=============

This distribution includes:

(a) Five different files named covering-stoyan.f90,
covering-america.f90, covering-star.f90, covering-cesaro.f90, and
covering-minkowski.f90 associated with the five problems considered in
the paper. An additional file modamerica.f90 borrowed from the
distribution of Algencan 3.1.1 is also included.

(b) Each of the five problems has its own script to compile, run, and
reproduce the numerical experiments in the paper. The scripts are
run-stoyan, run-america, run-covering-star, run-cesaro, and
run-minkowski.

(c) The Scharge random number generator is incldued in file drand.f90

(d) There is also a file named vorintpols.f90 that contains
subroutines to deal with Voronoi cells and convex polygons. It
includes the Sutherland-Hodgman algorithm for 2D clipping polygons and
an algorithm to compute the intersection of a convex polygon with a
circle.

HOW TO
======

To solve any of the five problems considered in the manuscript, simple
run the corresponding script and check the output in the screen.
