#/bin/bash
RAMASTER=/home/feathern/devel/forks/Rayleigh
cp -r $RAMASTER/src .
PATCHDIR=$PWD
cd $RAMASTER
make distclean
./configure -mkl --with-custom=$PATCHDIR/patch_nick2 --FC=mpifort -conda-mkl --prefix=$PATCHDIR $1
make -j4
make install
