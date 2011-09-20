for file in `ls *.sam`
do
    day=`ls -lt $file|cut -d' ' -f7`
    echo $day
    if [[ $day = "23" ]]
    then
        echo "$file matched"
	mv $file /media/disk-1/illu_100x_sams_4_5/
    fi
done
