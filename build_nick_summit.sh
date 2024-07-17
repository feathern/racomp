#/bin/bash
RAMASTER=/home/feathern/devel/Rayleigh_copy
PATCHDIR=$PWD
cp patch/* patch_nick/.
cd $RAMASTER
make distclean
./configure -mkl --with-custom=$PATCHDIR/patch_nick2 --FC=mpif90 --prefix=$PATCHDIR
make -j8
make install
