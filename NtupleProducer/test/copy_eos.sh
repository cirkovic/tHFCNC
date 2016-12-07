N=8
PSs=""
n=0
for i in /eos/cms/store/caf/user/mdjordje/Cirkovic/FCNC/ST_Hct_NEW /eos/cms/store/caf/user/mdjordje/Cirkovic/FCNC/ST_Hut_NEW /eos/cms/store/caf/user/mdjordje/Cirkovic/FCNC/ST_Hut_NEW /eos/cms/store/caf/user/mdjordje/Cirkovic/FCNC/TTbar_AntiTopLep_Hct_NEW /eos/cms/store/caf/user/mdjordje/Cirkovic/FCNC/TTbar_AntiTopLep_Hut_NEW /eos/cms/store/caf/user/mdjordje/Cirkovic/FCNC/TTbar_SM_NEW /eos/cms/store/caf/user/mdjordje/Cirkovic/FCNC/TTbar_TopLep_Hct_NEW /eos/cms/store/caf/user/mdjordje/Cirkovic/FCNC/TTbar_TopLep_Hut_NEW ; do
    O1=${i/FCNC/FCNC1}
    O2=${O1/_NEW/_FlatTree}
    for j in `eos ls ${i}`; do
        eos cp ${i}/${j} ${O2}/${j} &
        P="$!"
        if [ "$n" -lt "$N" ]; then
            PSs="$PSs $P"
            n=$(( n + 1 ))
        else
            wait $PSs
            PSs=""
            n=0
        fi
    done
done

