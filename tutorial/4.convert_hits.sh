#!/bin/bash 

if [ $MRMPIBLAST_PREFIX == "" ]; then
    echo "The environment variable, MRMPIBLAST_PREFIX is not set."
    exit
fi

echo -e "\n\n### Convert hit files into CSV file ###"
cd $MRMPIBLAST_PREFIX/tutorial/blastdb &&
mkdir -p hits &&
#mv *log.txt ./logs &&
mv *.bin ./hits &&
#python ../tools/load_hd5.py ./hits/ hits 1 &&
#python ../tools/load_sql.py ./hits/ hits 1 &&
python ../tools/load_csv.py ./hits/ hits &&

echo -e "\n### Done!"

# EOF
