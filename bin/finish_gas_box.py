#!/usr/bin/python
#take amr2map gas plot, add snapshot and scale info, draw colorbar

import os
import sys
import initialize as my

# check syntax
i=len(sys.argv)
if i!=6:
    print "Usage: finish_gas_box.py infile outfile snap aout radius[kpc/h]"
    exit(1)

os.nice(1)

from PIL import Image
import ImageDraw
import ImageFont

inm = Image.open(sys.argv[1])
x, y = inm.size

off = 5
im = Image.new("RGB", (x+2*off, y+2*off), "black")
draw = ImageDraw.Draw(im)
im.paste(inm, (off,off,x+off,y+off))

# assumption: box is across 1Mpc/h, with refinement in the inner 500kpc/h box
# this refinement box gets a scale bar
xstart = off + x/4
xend   = off + 3*x/4
ystart = yend = off + 6*y/8
draw.line((xstart, ystart, xend, yend), fill=0, width=3)


fout = sys.argv[2]
im.save(fout, "PNG")

sn   = sys.argv[3]
aout = float(sys.argv[4])
z    = 1.0/aout - 1.0
#text = "snap "+sn+",  a = "+str(aout) +",  z = "+str(z)
text = "  z = "+"%.2f"%z
my.txt2img(fout,text)

text = sys.argv[5]+" kpc/h"
my.txt2imgpos(fout,text,xstart,ystart)



del draw
