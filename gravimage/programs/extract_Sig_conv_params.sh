#!/bin/bash
# get all nu n(r) parameters, assuming pops == 1:

# call with output directory, phys_live.points is appended to that name

FILE=$1/phys_live.points
FILE2=$1/phys_live.points_cleaned
FILE3=$1/Sig_conv.params
FILE4=$1/Sig_conv.stats

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

# get important numbers via perl command

# perl -e 'use List::Util qw(max min sum); @a=();while(<>){$sqsum+=$_*$_; push(@a,$_)}; $n=@a;$s=sum(@a);$a=$s/@a;$m=max(@a);$mm=min(@a);$std=sqrt($sqsum/$n-($s/$n)*($s/$n));$mid=int @a/2;@srtd=sort @a;if(@a%2){$med=$srtd[$mid];}else{$med=($srtd[$mid-1]+$srtd[$mid])/2;};print "records:$n\nsum:$s\navg:$a\nstd:$std\nmed:$med\max:$m\min:$mm";'

# to be used like
# cat file_with_one_number_per_line|perl ........

echo -n >  $FILE4
for i in {1..21}
do
    #echo "say hello" $i
    stats=$(awk -v var=$i '{print $var}' $FILE3 | perl -e 'use List::Util qw(max min sum); @a=();while(<>){$sqsum+=$_*$_; push(@a,$_)}; $n=@a;$s=sum(@a);$a=$s/@a;$m=max(@a);$mm=min(@a);$std=sqrt($sqsum/$n-($s/$n)*($s/$n));$mid=int @a/2;@srtd=sort @a;if(@a%2){$med=$srtd[$mid];}else{$med=($srtd[$mid-1]+$srtd[$mid])/2;};print "$n,$s,$a,$std,$med,$m,$mm,";')

    myn=$(echo $myn $(echo $stats|cut -d',' -f1))
    mys=$(echo $mys $(echo $stats|cut -d',' -f2))
    mya=$(echo $mya $(echo $stats|cut -d',' -f3))
    mystd=$(echo $mystd $(echo $stats|cut -d',' -f4))
    mymed=$(echo $mymed $(echo $stats|cut -d',' -f5))
    mym=$(echo $mym $(echo $stats|cut -d',' -f6))
    mymm=$(echo $mymm $(echo $stats|cut -d',' -f7))
    #echo $stats >> $FILE4
    #vec=$(awk '{ print $$(echo $i) }' $FILE3)
    #echo $vec
    #echo "bogus"
    #exit 0
done


echo $mymm >> $FILE4
echo $mymed >> $FILE4
echo $mym >> $FILE4
