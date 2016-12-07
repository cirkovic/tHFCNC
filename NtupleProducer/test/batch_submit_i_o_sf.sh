INPUT=${1}
OUTPUT=${2}

NUM=0
#TOT=10
TOT=-1
while read p; do
    NUM=$((NUM+1))
    #echo $p
    #echo `sed "${NUM}q;d" $OUTPUT`
    CMD="bsub -q 1nh batch_job_sf.sh $p "`sed "${NUM}q;d" $OUTPUT`
    echo $CMD
    #exit
    eval $CMD
    if [[ $NUM -eq $TOT ]]; then
        break
    fi
done < ${INPUT}
