#!/usr/bin/env python2.7
import os, os.path, pdb
import sys, time, threading
import PIL
import PIL.Image, PIL.ImageDraw, PIL.ImageFont
import MySQLdb
import hashlib
import mysql as mys
nprocessor = 12
semaphore = threading.Semaphore(nprocessor)


def mv(a,b):
    os.system("mv "+a+" "+b+"&>/dev/null")
    return
## \fn mv(a, b)
# move a file a to b
# a string filename
# b string filename or folder

def read_a(snap):
    num5=str(snap).zfill(5)
    filename=mys.d(snap)+"/info_"+num5+".txt"
    if(not os.path.exists(filename)): exit(1)
    fi = open(filename)
    i=0;
    for line in fi:
        i=i+1
        if(i==10):
            val = line.split()
            a = float(val[2])
            break
    fi.close()
    return a
## \fn read_a(snap)
# read expansion factor of an output
# @param snap int snapshot

def read_z(snap):
    num5=str(snap).zfill(5)
    filename=mys.d(snap)+"/info_"+num5+".txt"
    if(not os.path.exists(filename)): exit(1)
    fi = open(filename)
    i=0
    for line in fi:
        i=i+1
        if(i==10):
            val = line.split()
            a = float(val[2])
            break
    fi.close()
    return 1.0/a-1.0
## \fn read_z(snap)
# read redshift of an output
# @param snap int snapshot number

def md5(snap,hid):
    m = hashlib.md5()
    m.update(str(snap*1e8+hid))
    return m.hexdigest()
## \fn md5(snap, hid)
# get (almost) unique name from snapshot and hid
# @param snap int snapshot number
# @param hid int TODO

def run_command(cmd):
    with semaphore:
        os.system(cmd)
## \fn run_command(cmd)
# run system command
# @param cmd string to be executed

def thread(cmd):
    #print("threading: ",cmd)
    threading.Thread(target=run_command, args=(cmd, )).start()
## \fn thread(cmd)
# start a command as a new thread
# @param cmd string of a command

def runif(cmd1, cmd2, calc, show, runi):
    if(calc):
        cmd = cmd1 + " && " + cmd2
    else:
        cmd = cmd2
    if(show):
        print(cmd)
    if(runi):
        os.system(cmd)
## \fn runif(cmd1, cmd2, calc, show, runi)
# run cmd1 if calc evaluates to true, else cmd2
# @param cmd1 string for os.system
# @param cmd2 string for os.system
# @param calc bool
# @param show bool
# @param runi bool

def threadif(cmd1, cmd2, run1, run2, show, runi):
    # cmd = ""
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
        print(cmd)
    if(runi):
        thread(cmd)
## \fn threadif(cmd1, cmd2, run1, run2, show, runi)
# run in background cmd1 or cmd2 based on run1, run2
# @param cmd1 string for os.system
# @param cmd2 string for os.system
# @param run1 bool
# @param run2 bool
# @param show bool
# @param runi bool

def run(cmd):
    print(cmd)
    os.system(cmd)
## \fn run(cmd)
# call os.system
# @param cmd string

def countlines(filename):
    lines = 0
    for line in open(filename):
        lines += 1
    return lines
## \fn countlines(filename)
# count lines in file, old way
# TODO: update with faster version of gravimage
# @param filename

def mkdir(path):
    os.system("mkdir -p "+path)
## \fn mkdir(path)
# create directory with parents
# @param path string

def txt2img(fi,text,bg="#ffffff",fg="#000000",font="FreeSans.ttf",FontSize=14):
    font_dir = "/usr/share/fonts/truetype/freefont/"
    img_name = fi#+".jpg"
    font_size = FontSize
    fnt = PIL.ImageFont.truetype(font_dir+font, font_size)
    lineWidth = 20
    img = PIL.Image.open(fi)
    # make an entirely black image
    imgbg = PIL.Image.new('RGBA', img.size, "#000000")
    # make a mask that masks out all
    mask = PIL.Image.new('L',img.size,"#000000")
    # setup to draw on the main image
    draw = PIL.ImageDraw.Draw(img)
    # setup to draw on the mask
    drawmask = PIL.ImageDraw.Draw(mask)
    # draw a line on the mask to allow some bg through
    drawmask.line((0, lineWidth/2, img.size[0],lineWidth/2),
                  fill="#999999", width=10)
    # put the (somewhat) transparent bg on the main
    img.paste(imgbg, mask=mask)
    # add some text to the main
    draw.text((10,0), text, font=fnt, fill=bg)
    del draw
    img.save(img_name,"PNG",quality=100)#"JPEG",quality=100)
## \fn txt2img(fi, text, bg, fg, font, FontSize)
# call pillow to draw text to canvas
# @param fi figure
# @param text string
# @param bg color
# @param fg color
# @param font string
# @param FontSize int

def txt2imgpos(fi, text, x, y, bg="#ffffff", fg="#000000", font="FreeSans.ttf", FontSize=14):
    font_dir = "/usr/share/fonts/truetype/freefont/"
    img_name = fi#+".jpg"
    font_size = FontSize
    fnt = PIL.ImageFont.truetype(font_dir+font, font_size)
    # lineWidth = 20
    img = PIL.Image.open(fi)
    # setup to draw on the main image
    draw = PIL.ImageDraw.Draw(img)
    # add some text to the main
    draw.text((x,y), text, font=fnt, fill=bg)
    del draw
    img.save(img_name,"PNG",quality=100)#"JPEG",quality=100)
## \fn txt2imgpos(fi, text, x, y, bg, fg, font, FontSize)
# draw text at specific position of figure
# @param fi figure
# @param text string
# @param x int
# @param y int
# @param bg color
# @param fg color
# @param font string
# @param FontSize int

def open_file(filename, mode):
        """Open a file."""
        if(not os.path.exists(filename)):
            print('file ',filename,' not found')
            sys.exit()
        try:
            the_file = open(filename, mode)
            print("filesize: ",os.path.getsize(filename))
        except(IOError) as e:
            print("Unable to open the file", the_file, "Ending program.\n", e)
            #raw_input("\n\nPress the enter key to exit.")
            sys.exit()
        else:
            return the_file
## \fn open_file(filename, mode)
# open file
# @param filename string
# @param mode string

def file_exists(filename):
    return os.path.exists(filename)
## \fn file_exists(filename)
# does file exist already?
# @param filename string

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
## \fn get_xyzr(snap)
# get center for a given snapshot
# @param snap int

def sqlstart():
    # ask for username and password
    #user = raw_input("Enter mySQL username: ")
    #somesecstring = raw_input(" .. and his little secret: ")
    user="root"
    somesecstring = "Sinux10"
    #user = "mangaloi_astro"
    #somesecstring = "Pinux10"
    #connection = MySQLdb.connect('108.167.143.108', user, somesecstring, user)
    connection = MySQLdb.connect('127.0.0.1', user, somesecstring, "mangaloi_astro")
    cursor = connection.cursor()
    return connection, cursor
## \fn sqlstart()
# start connection

def sqlstop(connection, cursor):
    cursor.close()
    connection.close()
    return
## \fn sqlstop(connection, cursor)
# @param connection
# @param cursor

def sqlcu(cmd, co, cu):
    try:
        cu.execute(cmd)
        co.commit() # THIS MAKES ALL CHANGES GO TO THE MYSQL SERVER! IF NOT USED: STALE DATA
        out = cu.fetchall()
        return out
    except MySQLdb.Error as e:
        print("Error %d: %s" % (e.args[0],e.args[1]))
        sys.exit(1)
## \fn sqlcu(cmd, co, cu)
# execute a SQL command
# @param cmd string to be executed, in MySQl format
# @param co connection
# @param cu cursor

def sql(cmd):
    co, cu = sqlstart()
    try:
        cu.execute(cmd)
        co.commit() # THIS MAKES ALL CHANGES GO TO THE MYSQL SERVER! IF NOT USED: STALE DATA
        out = cu.fetchall()
        sqlstop(co, cu)
        return out
    except MySQLdb.Error as e:
        print("Error %d: %s" % (e.args[0],e.args[1]))
        sqlstop(co,cu)
        sys.exit(1)
## \fn sql(cmd)
# execute cmd with sql, starting connection inside function
# @param cmd command, string

def sqlout(cmd):
    row = sql(cmd)
    numrows=len(row)
    for i in range(numrows):
        print(row[i][0])
    return
## \fn sqlout(cmd)
# show output
# @param cmd string for SQL to be executed within function

def done():
    while(threading.active_count()>1):
        time.sleep(1)
    print("all done, ",threading.active_count())
    print(" ")
    return
## \fn done()
# print done after end of threading
