#/bin/bash
RAMASTER=/projects/frra1220/rayleigh_libs/Rayleigh
PATCHDIR=$PWD
cp patch/* patch_fredy/.
cd $RAMASTER
make distclean
./configure -devel -mkl --with-custom=$PATCHDIR/patch_fredy
make
