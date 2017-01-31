rm EOS_input_*.txt EOS_output_*.txt

for i in FCNC6 FCNC7 FCNC8 FCNC9 FCNC10; do
    #eos mkdir /eos/cms/store/caf/user/mdjordje/Cirkovic/NPO/${i}
    for j in `eos ls /eos/cms/store/caf/user/mdjordje/Cirkovic/${i}`; do
        #eos mkdir /eos/cms/store/caf/user/mdjordje/Cirkovic/NPO/${i}/${j}_Ntuple
        DIR=/eos/cms/store/caf/user/mdjordje/Cirkovic/${i}/${j}_Ntuple
        if [[ "${DIR}"  == *"_Ntuple_Ntuple" ]]; then
            continue
        fi
        if [[ `eos ls ${DIR//_Ntuple/}`  == "" ]]; then
            continue
        fi
        #COMMAND="eos mkdir ${DIR}"
        #COMMAND="eos ls ${DIR//_Ntuple/} | wc -l"
        #echo $COMMAND
        #eval $COMMAND
        txt_in=EOS_input_${i}_`basename ${DIR//_Ntuple/}`.txt
        touch $txt_in
        for k in `eos ls ${DIR//_Ntuple/}`; do
            echo ${DIR//_Ntuple/}/${k} >> $txt_in
        done
        txt_out=${txt_in/_input_/_output_}
        cp $txt_in $txt_out
        sed -i "s/Cirkovic\//Cirkovic\/NPO\//g" $txt_out
        sed -i "s/\/output\_/\_Ntuple\/output\_/g" $txt_out
        #exit
    done
done

rm *_FCNC6_*_x.txt
