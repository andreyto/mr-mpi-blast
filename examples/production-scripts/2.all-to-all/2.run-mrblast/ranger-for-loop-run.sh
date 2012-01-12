#!/bin/sh

Directory="./"
for direc in $Directory* ; do
    bDone=0
    if [[ -d $direc ]]; then
        echo $direc
        FileList=$(find $direc -type f)
        for file2 in $FileList ; do
            #echo $file2
            if [[ $file2 == *SUCCESS ]] ; then
                echo "already done in $direc"
                bDone=1
                break
            fi
        done
        
        if [[ bDone -eq 0 ]]; then
            topdir=$(pwd)
            cd $direc            
            rm -rf output-*
            echo "run job in $direc"
            time ibrun ../mrblast -db - -dbsize 12494903041 -evalue 1e-4 -num_threads 1 -window_size 0 -word_size 11 -searchsp 0 -num_descriptions 500 -num_alignments 10000 -penalty -5 -reward 4  -lcase_masking -dust yes -soft_masking true -max_target_seqs 2147483647
            echo "job done in $direc"
	    cp ../*.o* .
	    cp ../*.e* .
            touch SUCCESS
            cd $topdir
        fi
    fi
done

#Directory="./"
#for direc in $Directory* ; do
    #echo $direc
    #if find $direc -name SUCCESS; then
        #pushd $direc
        #echo $(pwd)
        #echo "run job in $direc"
        ##ibrun ...
        ##sleep 10; touch OUTPUT
        #echo "job done in $direc"
        #touch SUCCESS
        #popd
    #fi
#done


