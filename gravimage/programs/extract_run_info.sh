#!/bin/bash


var=$(ls *|grep -E ^201[1-4]|cut -d":" -f1)
count=1
for i in $var
do
    if [ ! -f $i/ev.dat ]; then
        echo "File "$i"/ev.dat not found!"
    else
        if [ ! -f $i/programs/gl_params.py ]; then
            echo "File "$i"/programs/gl_params.py not found!"
        else
            lines=$(wc -l $i/ev.dat|cut -d" " -f1)
            nbeta=$(grep "self.nbeta =" $i/programs/gl_params.py | cut -d"=" -f2 | cut -d"#" -f1)
            pops=$(grep "self.pops =" $i/programs/gl_params.py | cut -d"=" -f2 | cut -d"#" -f1)
            bins=$(grep "self.nipol =" $i/programs/gl_params.py | cut -d"=" -f2 | cut -d"#" -f1)

            echo -e $count"\t"$i"\t"$lines"\t"$pops"\t"$bins"\t"$nbeta
            count=$(echo $count"+1"|bc)
        fi
    fi
done
