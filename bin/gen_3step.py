#!/usr/bin/python
#concatenate 6 pngs on one wall

import os
import sys
import initialize as my

# check syntax
i=len(sys.argv)
if i!=5:
    print "Usage: pngwall.py dm_contour gas_ov dm_prof outfile"
    exit(1)

os.nice(1)

import initialize as my

from PIL import Image
import ImageDraw
import ImageFont

rr = Image.open(sys.argv[3])

x, y = rr.size
x2=x/2; y2=y/2

LU = Image.open(sys.argv[1])
lu = LU.resize((x2, y2), Image.ANTIALIAS)
LB = Image.open(sys.argv[2])
lb = LB.resize((x2, y2), Image.ANTIALIAS)

off = 5
im = Image.new("RGB", (x+x2+3*off,y+2*off), "black")
draw = ImageDraw.Draw(im)

im.paste(lu, (off,off,x2+off,y2+off))
im.paste(lb, (off,y2+off,x2+off,y+off))
im.paste(rr, (x2+2*off,off,x+x2+2*off,y+off))

fil = sys.argv[4]
im.save(fil, "PNG")
