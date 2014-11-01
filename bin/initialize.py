nprocessor=32

import os
import sys
import threading
import Image
import ImageDraw
import ImageFont
import MySQLdb
semaphore = threading.Semaphore(nprocessor)

import hashlib
import mys

def mv(a,b):
    os.system("mv "+a+" "+b+"&>/dev/null")

def read_a(snap):
    num5=str(snap).zfill(5)
    filename=mys.d(snap)+"/info_"+num5+".txt"
    if(not os.path.exists(filename)): exit(1)
    file = open(filename)
    i=0; z=100.0
    for line in file:
        i=i+1
        if(i==10):
            val = line.split()
            a = float(val[2])
            break
    file.close()
    return a

def read_z(snap):
    num5=str(snap).zfill(5)
    filename=mys.d(snap)+"/info_"+num5+".txt"
    if(not os.path.exists(filename)): exit(1)
    file = open(filename)
    i=0; z=100.0
    for line in file:
        i=i+1
        if(i==10):
            val = line.split()
            a = float(val[2])
            break
    file.close()
    return 1.0/a-1.0


def md5(snap,hid):
    m = hashlib.md5()
    m.update(str(snap*1e8+hid))
    return m.hexdigest() 

def print_usage():
    print "usage: master.py [action]"
    print "            0> prepare output folder structure"
    print "            1> ramses2gadget"
    print "            2> AHFstep, rename output, create halo file"
    print "            3> halo"
    print "            4> get spheres"
    print "            5> correct centering"
    print "            6> get particles in AOI"
    print "          7-9> get merger tree, evolution, interpolate missing snapshots"
    print "           10> plot DM density"
    print "           11> plot gas density"
    print "           12> plot star density"
    print "           13> plot DM/gas/star profile"
    print "           14> compile wall"
    print "           15> compile movies"
    print "if no action is given, all of them are performed"
    print "..."

def run_command(cmd):
    with semaphore:
        os.system(cmd)

def thread(cmd):
    print "threading: ",cmd
    threading.Thread(target=run_command, args=(cmd, )).start()

def runif(cmd1, cmd2, calc, show, runi):
    if(calc):
        cmd = cmd1 + " && " + cmd2
    else:
        cmd = cmd2
    if(show):
        print cmd
    if(runi):
        os.system(cmd)

        
def threadif(cmd1, cmd2, run1, run2, show, runi):
    cmd = ""
    if(run1):
        brun1 = "true"
    else:
        brun1 = "false"
    if(run2):
        brun2 = "true"
    else:
        brun2 = "false"
    cmd = "(("+brun1+" && "+cmd1+")||true)"
    cmd = cmd + " && "
    cmd = cmd + "("+brun2+" && "+cmd2+")"
    if(show):
        print cmd
    if(runi):
        thread(cmd)

def run(cmd):
    print cmd
    os.system(cmd)

def countlines(filename):
    lines = 0
    for line in open(filename):
        lines += 1
    return lines

def mkdir(path):
    os.system("mkdir -p "+path)

def txt2img(fi,text,bg="#ffffff",fg="#000000",font="FreeSans.ttf",FontSize=14):
    font_dir = "/usr/share/fonts/truetype/freefont/"
    img_name = fi#+".jpg"
    font_size = FontSize
    fnt = ImageFont.truetype(font_dir+font, font_size)
    lineWidth = 20
    img = Image.open(fi)
    # make an entirely black image
    imgbg = Image.new('RGBA', img.size, "#000000")
    # make a mask that masks out all
    mask = Image.new('L',img.size,"#000000")
    # setup to draw on the main image
    draw = ImageDraw.Draw(img)
    # setup to draw on the mask
    drawmask = ImageDraw.Draw(mask)
    # draw a line on the mask to allow some bg through
    drawmask.line((0, lineWidth/2, img.size[0],lineWidth/2),
                  fill="#999999", width=10)
    # put the (somewhat) transparent bg on the main
    img.paste(imgbg, mask=mask)
    # add some text to the main
    draw.text((10,0), text, font=fnt, fill=bg)      
    del draw
    img.save(img_name,"PNG",quality=100)#"JPEG",quality=100)

def txt2imgpos(fi,text,x,y,bg="#ffffff",fg="#000000",font="FreeSans.ttf",FontSize=14):
    font_dir = "/usr/share/fonts/truetype/freefont/"
    img_name = fi#+".jpg"
    font_size = FontSize
    fnt = ImageFont.truetype(font_dir+font, font_size)
    lineWidth = 20
    img = Image.open(fi)

    # setup to draw on the main image
    draw = ImageDraw.Draw(img)

    # add some text to the main
    draw.text((x,y), text, font=fnt, fill=bg)
    del draw
    img.save(img_name,"PNG",quality=100)#"JPEG",quality=100)

def open_file(file_name, mode):
        """Open a file."""
        try:
            the_file = open(file_name, mode)
            print "filesize: ",os.path.getsize(file_name)
        except(IOError), e:
            print "Unable to open the file", file_name, "Ending program.\n", e
            #raw_input("\n\nPress the enter key to exit.")
            sys.exit()
        else:
            return the_file

def file_exists(filename):
    return os.path.exists(filename)

def get_xyzr(snap):
    i=0
    f = open("mt/cen","r")
    for lines in f:
        i=i+1
        if(not i==snap):
            continue
        val = lines.split()
        x = float(val[0])
        y = float(val[1])
        z = float(val[2])
        r = float(val[3])

    return x,y,z,r

def sqlstart():
    connection = MySQLdb.connect('igloo.dhcp.phys','psteger','Pinux10','astro')
    cursor = connection.cursor ()
    return connection,cursor

def sqlstop(connection,cursor):
    cursor.close()
    connection.close()

def sql(cmd):
    co,cu=sqlstart()
#    print cmd
    try:
        cu.execute(cmd)
        out=cu.fetchall()
        sqlstop(co,cu)
        return out
    except MySQLdb.Error, e:
        print "Error %d: %s" % (e.args[0],e.args[1])
        sqlstop(co,cu)
        sys.exit(1)

def sqlout(cmd):
    row = sql(cmd)
    numrows=len(row)
    for i in range(numrows):
        print row[i][0]
