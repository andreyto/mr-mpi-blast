for file in `ls *.fasta`
do
    day=`ls -lt $file|cut -d' ' -f7`
    echo $day
    if [[ $day = "23" ]]
    then
        echo "$file matched"
	mv $file /media/disk/illu_100x_reads_4/
    fi
done
