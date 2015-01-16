#!/usr/bin/python
#concatenate 6 pngs on one wall

import os
import sys
import initialize as my

# check syntax
i=len(sys.argv)
if i!=10:
    print "Usage: pngwall.py UL, UM, UR, BL, BM, BR, outfile, snap, aout"
    exit(1)

os.nice(1)

import initialize as my

from PIL import Image
import ImageDraw
import ImageFont

ul = Image.open(sys.argv[1])
um = Image.open(sys.argv[2])
ur = Image.open(sys.argv[3])
bl = Image.open(sys.argv[4])
bm = Image.open(sys.argv[5])
br = Image.open(sys.argv[6])
x0,y0= ul.size
x, y = um.size
x0 = x0*y/y0

#NEAREST, BILINEAR, BICUBIC, ANTIALIAS
ul2 = ul.resize((x0,y), Image.ANTIALIAS)
um2 = um.resize((x, y), Image.ANTIALIAS)
ur2 = ur.resize((x, y), Image.ANTIALIAS)
bl2 = bl.resize((x0,y), Image.ANTIALIAS)
bm2 = bm.resize((x, y), Image.ANTIALIAS)
br2 = br.resize((x, y), Image.ANTIALIAS)

off = 5
im = Image.new("RGB", (2*x+x0+4*off,2*y+3*off), "black")
draw = ImageDraw.Draw(im)
im.paste(ul2, (off,off,x0+off,y+off))

im.paste(um2, (x0+2*off,off,x+x0+2*off,y+off))
draw.line((x0+2*off, off+y/2, x+x0+2*off, off+y/2), fill=128)
draw.line((x0+x/2+2*off, off, x0+x/2+2*off, off+y), fill=128)

im.paste(ur2, (x0+x+3*off,off,x0+2*x+3*off,y+off))
draw.line((x0+x+3*off, off+y/2, x0+2*x+3*off, off+y/2), fill=128)
draw.line((x0+1.5*x+3*off, off, x0+1.5*x+3*off, off+y), fill=128)

im.paste(bl2, (off,y+2*off,x0+off,2*y+2*off))
im.paste(bm2, (x0+2*off,y+2*off,x0+x+2*off,2*y+2*off))

im.paste(br2, (x0+x+3*off,y+2*off,x0+2*x+3*off,2*y+2*off))
draw.line((x0+x+3*off, 2*off+1.5*y, x0+2*x+3*off, 2*off+1.5*y), fill=128)
draw.line((x0+1.5*x+3*off, 2*off+y, x0+1.5*x+3*off, 2*off+2*y), fill=128)


fil = sys.argv[7]
sn   = sys.argv[8]
aout = float(sys.argv[9])
z    = 1.0/aout - 1.0
text = "snap "+sn+",  a = "+str(aout) +",  z = "+str(z)
im.save(fil, "PNG")
my.txt2img(fil,text)
del draw
