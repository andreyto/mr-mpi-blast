#!/bin/sh

Directory="./"
for direc in $Directory* ; do
    bDone=0
    if [[ -d $direc ]]; then
        #echo $direc
        numFiles=$(ls $direc/*.bin | wc)        
        echo $numFiles
    fi
done
