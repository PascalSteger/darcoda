#/usr/bin/env sh
#echo OFF
echo "Profiling $1 ..."
echo "with parameters " $*
python3 -m cProfile -o $1.pstats $*
echo "Generating image ..."
gprof2dot.py -n 0.01 -e 0.01 -f pstats $1.pstats | dot -Gcharset=latin1 -Tpng -o $1.png
echo "Displaying image ..."
gpicview $1.png

