#!/bin/bash

if [ -d /home/ast ]; then
    cd /home/ast/read/dark/darcoda/gravimage
else
    cd /home/psteger/sci/darcoda/gravimage
fi

extract_active_runs.sh > active_runs

# enable extended globbing
shopt -s extglob

for inv in DT*
do
    for case in $inv/?([1,2])[0-9]/
    do
        #echo  $case
        var=$(ls $case* |grep -E /201[0-9]|cut -d":" -f1)
        count=1
        for timestamp in $var
        do
            #echo $timestamp
            if [ ! -f $timestamp/ev.dat ]; then
                echo "File "$timestamp"/ev.dat not found!"
            else
                if [ ! -f $timestamp/programs/gi_params.py ]; then
                    echo "File "$timestamp"/programs/gi_params.py not found!"
                else
                    lines=$(wc -l $timestamp/ev.dat|cut -d" " -f1)
                    nbeta=$(grep "self.nbeta =" $timestamp/programs/gi_params.py | cut -d"=" -f2 | cut -d"#" -f1)
                    pops=$(grep "self.pops =" $timestamp/programs/gi_params.py | cut -d"=" -f2 | head -n1 | cut -d"#" -f1)
                    bins=$(grep "self.nipol =" $timestamp/programs/gi_params.py | cut -d"=" -f2 | cut -d"#" -f1)
                    beta00=$(grep "self.beta00prior =" $timestamp/programs/gi_params.py | cut -d"=" -f2 | cut -d"#" -f1)
                    minbeta=$(grep "self.minbetastar =" $timestamp/programs/gi_params.py | cut -d"=" -f2 | cut -d"#" -f1)

                    # if found in active_runs, append "a"
                    found=$(grep $timestamp active_runs)
                    active=" "
                    if [ $timestamp = $found"" ]; then
                        active="a"
                    fi

                    # if already plotted, append "p"
                    plotted=" "
                    if [ -f $timestamp/output/prof_chi2_0.pdf ]; then
                        plotted="p"
                    fi

                    # add converged flag
                    conv=" "
                    if [ -f $timestamp/output/converged ]; then
                        conv="c"
                    fi
                    echo -e $count"\t"$timestamp"\t"$lines"\t"$pops"\t"$bins"\t"$nbeta"\t"$beta00"\t"$minbeta"\t"$active"\t"$plotted"\t"$conv

                    count=$(echo $count"+1"|bc)
                fi
            fi
        done
    done
done
