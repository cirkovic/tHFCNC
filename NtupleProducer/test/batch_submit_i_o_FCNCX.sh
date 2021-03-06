INPUT=${1}
OUTPUT=${2}
TABLE=${3}

NUM=0
#TOT=10
TOT=-1
while read p; do
    NUM=$((NUM+1))
    #echo $p
    #echo `sed "${NUM}q;d" $OUTPUT`
    OUT=`sed "${NUM}q;d" $OUTPUT`
    CMD="bsub -q 1nh batch_job_FCNCX.sh $p $OUT $TABLE"
    if [[ `eos ls $OUT` != "" ]]; then
        continue
    fi
    echo $CMD
    eval $CMD
    if [[ $NUM -eq $TOT ]]; then
        break
    fi
done < ${INPUT}
