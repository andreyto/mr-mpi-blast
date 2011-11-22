#!/bin/bash 

if [ $MRMPIBLAST_PREFIX == "" ]; then
    echo "The environment variable, MRMPIBLAST_PREFIX is not set."
    exit
fi


mkdir -p $MRMPIBLAST_PREFIX/bin2 &&
cp $MRMPIBLAST_PREFIX/bin/mrblast $MRMPIBLAST_PREFIX/bin/makeblastdb $MRMPIBLAST_PREFIX/bin2/ &&
rm -rf $MRMPIBLAST_PREFIX/bin &&
mv $MRMPIBLAST_PREFIX/bin2 $MRMPIBLAST_PREFIX/bin &&
rm -rf $MRMPIBLAST_PREFIX/include &&
rm -rf $MRMPIBLAST_PREFIX/lib

# EOF
