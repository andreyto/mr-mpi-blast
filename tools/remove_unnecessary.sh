#!/bin/bash 

if [ -z "$MRMPIBLAST_PREFIX" ]; then
    echo "The environment variable, MRMPIBLAST_PREFIX is not set."
    exit
fi

echo -e "\n\n### Remove unnecessary files"
find $MRMPIBLAST_PREFIX/bin/ -maxdepth 1 -type f ! -iname 'mrblast' -and ! -iname 'makeblastdb' -delete &&
rm -rf $MRMPIBLAST_PREFIX/include &&
rm -rf $MRMPIBLAST_PREFIX/lib

echo -e "\n### Done!"
# EOF
