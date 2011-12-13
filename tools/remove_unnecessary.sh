#!/bin/bash 

if [ -z "$MRMPIBLAST_PREFIX" ]; then
    echo "The environment variable, MRMPIBLAST_PREFIX is not set."
    exit
fi

echo -e "\n\n### Remove unnecessary files"
find $MRMPIBLAST_PREFIX/bin/ -maxdepth 1 -type f ! -name 'mrblast' -and ! -name 'makeblastdb' -and ! -name 'load_*.py' -and ! -name 'splitter.py' -and ! -name 'seqindexer.py' -and ! -name 'blastdbcmd' -exec rm -f {} \; &&
rm -rf $MRMPIBLAST_PREFIX/include &&
rm -rf $MRMPIBLAST_PREFIX/lib

echo -e "\n### Done!"
# EOF
