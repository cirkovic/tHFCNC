cd /afs/cern.ch/work/c/cirkovic/fcnc_ana/CMSSW_8_0_12/src/
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/c/cirkovic/fcnc_ana/CMSSW_8_0_12/src/tHFCNC/NtupleProducer/test
export LD_LIBRARY_PATH=/afs/cern.ch/work/c/cirkovic/fcnc_ana/CMSSW_8_0_12/src/tHFCNC/NtupleProducer:/afs/cern.ch/work/c/cirkovic/fcnc_ana/CMSSW_8_0_12/src/tHFCNC/NtupleProducer/obj:$LD_LIBRARY_PATH

echo ARGS: ${1} ${2} ${3}

#./NtupleProducer
#--file list.txt // txt file with the list of input FlatTree files to read
#--tree FlatTree/tree // name of FlatTree TTree directory
#--outfile "output" // name of the output ROOT file
#--noe 666 // initial number of events for considered dataset to be used for normalization
#--xsec 6.66 // cross section for considered dataset to be used for normalization
#--nmax -1 // max number of events to process (-1: no limit)
#--isdata 0 // run on data or MC
#--stream -1 // 0: SingleElectron, 1: SingleMuon (only for data)
#--issig 0 // flag to mark signal events (1: signal)

#./NtupleProducer --file lists/list_1.txt --tree FlatTree/tree --outfile "output" --noe 666 --xsec 6.66 --nmax -1 --isdata 0 --stream -1 --issig 0

while read p; do
    echo ${p}
    IFS=' ' read -r -a array <<< "$p"
    if [[ ${1} == *${array[0]}* ]] ; then

        ISDATA=0
        if [[ ${array[0]} == "/Single"*"/" ]] ; then
            ISDATA=1
        fi

        STREAM=-1
        if [[ ${array[0]} == "/SingleElectron"*"/" ]] ; then
            STREAM=0
        elif [[ ${array[0]} == "/SingleMuon"*"/"  ]] ; then
            STREAM=1
        fi

        echo ISDATA: $ISDATA
        echo STREAM: $STREAM

        MKTEMP=`mktemp`
        echo root://eoscms.cern.ch/${1} > $MKTEMP

        MKTEMPD=`mktemp -d`

        CMD="./NtupleProducer --file $MKTEMP --tree FlatTree/tree --outfile ${MKTEMPD}/output --noe ${array[2]} --xsec ${array[1]} --nmax -1 --isdata $ISDATA --stream $STREAM"
        #CMD="./NtupleProducer --file $MKTEMP --tree FlatTree/tree --outfile ${MKTEMPD}/output --noe ${array[2]} --xsec ${array[1]} --nmax 100 --isdata $ISDATA --stream $STREAM"
        #CMD="./NtupleProducer --file $MKTEMP --tree FlatTree/tree --outfile ${MKTEMPD}/output --noe ${array[2]} --xsec ${array[1]} --nmax 1000 --isdata $ISDATA --stream $STREAM"
        echo $CMD
        eval $CMD

        CMD="xrdcp -f ${MKTEMPD}/output.root root://eoscms.cern.ch/${2}"
        echo $CMD
        eval $CMD

        SFO="root://eoscms.cern.ch/${2}"
        SFO1=${SFO/output_/sf_}
        CMD="xrdcp -f ${MKTEMPD}/sf.root ${SFO1}"
        echo $CMD
        eval $CMD
        
        rm -rf $MKTEMP $MKTEMPD

    fi
done < ${3}

