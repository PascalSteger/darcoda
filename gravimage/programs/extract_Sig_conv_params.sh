#!/bin/bash
# get all nu n(r) parameters, assuming pops == 1:

# call with output directory, phys_live.points is appended to that name

FILE=$1/phys_live.points
FILE2=$1/phys_live.points_cleaned
FILE3=$1/Sig_conv.params

nrho=21
nbeta=4
pops=1

sed 's/[ ]\+/ /g' $FILE | sed 's/^ //g' > $FILE2

echo -n > $FILE3
while read params; do
    endoffset=$nrho
    #paramsrhonr=$(echo $params|cut -d" " -f1-$endoffset)
    #echo "params rho nr: "
    #echo $paramsrhonr
    #echo $paramsrhonr|wc -w
    endoffset2=$(echo $endoffset+$nrho|bc)
    endoffset=$(echo $endoffset+1|bc)
    paramsnunr=$(echo $params|cut -d" " -f$endoffset-$endoffset2)
    #echo "params nu nr: "
    echo $paramsnunr >> $FILE3
    #echo $paramsnunr|wc -w
    endoffset3=$(echo $endoffset2+$nbeta|bc)
    endoffset2=$(echo $endoffset2+1|bc)
    #paramsbeta=$(echo $params|cut -d" " -f$endoffset2-$endoffset3)
    #echo "params nu beta: "
    #echo $paramsbeta
    #echo $paramsbeta|wc -w
    #exit 0
done < "$FILE2"
