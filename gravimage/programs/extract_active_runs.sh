#!/bin/zsh

TODAY=$(date | cut -d" " -f2,3)
#echo $TODAY
ls -lt DT*/[0-12]/201*/phys_live.points | head -n 20| grep $TODAY | awk '{print $NF}' | cut -d"/" -f1-3|sort|uniq
