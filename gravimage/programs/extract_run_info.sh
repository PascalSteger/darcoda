#!/bin/bash

# enable extended globbing
shopt -s extglob

cd /home/ast/read/dark/darcoda/gravimage

for inv in DT*
do
    for case in $inv/?([1,2])[0-9]/
    do
        #echo  $case
        var=$(ls $case/* |grep -E /201[0-9]|cut -d":" -f1)
        count=1
        for timestamp in $var
        do
            #echo $timestamp
            if [ ! -f $timestamp/ev.dat ]; then
                echo "File "$timestamp"/ev.dat not found!"
            else
                if [ ! -f $timestamp/programs/gl_params.py ]; then
                    echo "File "$timestamp"/programs/gl_params.py not found!"
                else
                    lines=$(wc -l $timestamp/ev.dat|cut -d" " -f1)
                    nbeta=$(grep "self.nbeta =" $timestamp/programs/gl_params.py | cut -d"=" -f2 | cut -d"#" -f1)
                    pops=$(grep "self.pops =" $timestamp/programs/gl_params.py | cut -d"=" -f2 | cut -d"#" -f1)
                    bins=$(grep "self.nipol =" $timestamp/programs/gl_params.py | cut -d"=" -f2 | cut -d"#" -f1)

                    echo -e $count"\t"$timestamp"\t"$lines"\t"$pops"\t"$bins"\t"$nbeta
                    count=$(echo $count"+1"|bc)
                fi
            fi
        done
    done
done