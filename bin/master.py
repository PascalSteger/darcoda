#!/usr/bin/python
# generate movie of simulation#
# ASS: simulation performed
# ASS: sit in simulation output directory, sim snaps in output_*
# ASS: ramses2gadget = r2g, AHFstep, ../tools/bin in PATH

simdir="/scratch/psteger/sim/sim_aboley/"
simdir="/scratch/psteger/sim/nec_20111220/"
nstart = 2
nstop  = 68
num = 68  # total number of snapshots

import os
import sys
import initialize as my
import mys

brute = True; loop=True;
os.nice(1)

my.print_usage()

# get parameters
# no parameter: run all
# one param:    run specified
i = len(sys.argv)
if(i==1):
    for k in range(18):
        my.run("master.py "+str(k))
    sys.exit(0)
    
action = int(sys.argv[1])

if(action==-1):
    mys.clear()

if(action==0):
    print "prepare output folder structure, MySQL structure"
#    mys.setup()
#    mys.fill_sim(simdir,num,nstart,nstop)
    for ncounter in range(nstop-nstart+1):
        d=mys.d(nstart+ncounter)
#        d=simdir+"/output_"+str(ncounter+nstart).zfill(5)+"/";
        my.mkdir(d+'amr')
        my.mkdir(d+'grav')
        my.mkdir(d+'hydro')
        my.mkdir(d+'part')
        my.mv(d+"amr_*",d+"amr")
        my.mv(d+"grav_*",d+'/grav')
        my.mv(d+"hydro_*",d+'/hydro')
        my.mv(d+'part_*',d+'/part')

# convert to Gadget file format
if(action==1):
    for ncounter in range(nstop-nstart+1):
        nc = ncounter + nstart
        num5 = str(nc).zfill(5)
        d  = simdir + "/output_"+num5+"/"
        print "> ramses2gadget"
        my.mkdir(d+"/r2g")
        my.thread("r2g -i "+d)

# run AHFstep on all snapshots
if(action==2):
    print "> AHFstep, rename output, create halo file"
    for ncounter in range(nstop-nstart+1):
        nsnap = ncounter + nstart
        d  = mys.d(nsnap)

        f = open(d+"a.par","w")
        s = d+"r2g/r2g. "
        s += "61 1\nao\n16\n4\n4\n0\n0\n0\n0\n"
        f.write(s)
        f.close()

        cmd = "cd "+d+" && AHFstep a.par"
        cmd = cmd + " && mv "+d+"*_halos "+d+"halos"
        cmd = cmd + " && mv "+d+"*_centres "+"centres"
        cmd = cmd + " && mv "+d+"*_particles "+d+"particles"
        cmd = cmd + " && mv "+d+"*_profiles "+d+"profiles"
        cmd = cmd + " && mv "+d+"*_substructure "+d+"substructure"
        cmd = cmd + " && rm "+d+"halo"
        my.thread(cmd)

#fill AHF halo properties into mysql db, have column "hid" filled afterwards
if(action==3):
    for ncounter in range(nstop-nstart+1):
        nsnap = ncounter + nstart

        my.thread("fill_ahf_to_db.py "+str(nsnap))
        
# gen_spheres with centering on AHF centers
if(action==4):
    for ncounter in range(nstop-nstart+1):
        nsnap = ncounter + nstart
        d  = mys.d(nsnap)
        my.mkdir(d+"stars"); my.mkdir(d+"dm")
        cmd = "cd "+d+" && gen_spheres.py "+str(nsnap)+" 1"
        if(brute):
            my.thread(cmd)
        else:
            my.run(cmd)
        if(not loop):
            break

if(action==5):
    print "> find center with shrinking spheres / sph densities"
    # have xs,xs_star set in the end
    for ncounter in range(nstop-nstart+1):
        nsnap = ncounter + nstart

        # shrink_spheres
        print mys.get_nhalo(nsnap)
        for i in range(mys.get_nhalo(nsnap)):
            cmd = "shrink_sphere.py "+str(nsnap)+" "+str(i+1)
            my.thread(cmd)
            cmd = "shrink_sphere_stars.py "+str(nsnap)+" "+str(i+1)
            my.thread(cmd)

            #break;

if(action==6):
    print "> get particles in AOI"

    for ncounter in range(nstop-nstart+1):
        nsnap = ncounter + nstart
        d  = mys.d(nsnap)
        my.mkdir(d+"stars"); my.mkdir(d+"dm")
        cmd = "gen_spheres.py "+str(nsnap)+" 2"
        if(brute):
            my.thread(cmd)
        else:
            my.run(cmd)
        if(not loop):
            break

# generate all particles list (to track dm halos and star clusters)
if(action==7):
    for ncounter in range(nstop-nstart+1):
        nsnap = ncounter + nstart
        d  = mys.d(nsnap)
        # for all halos:
        cmd = "gen_particles_list.py "+str(nsnap)
        cmd += "&& gen_particles_dm.py "+str(nsnap)
        my.run(cmd)
#        if(not loop):
#            break

if(action==8):
    print "> get merger tree, evolution, interpolate missing snapshots"
    if(nstart==nstop):
        ns = str(nstart)
        cmd = "substructure.py "+ns
        my.thread(cmd)
    else:
        for nc in range(nstop-nstart+1):
            cmd = "merger_tree.py "+str(nstop-nc)+" "+str(nstop-1-nc)
            cmd+= "&& substructure.py "+str(nstop-nc)
            my.thread(cmd)
            if(not loop):
                break
            
if(action==9):
    # follow halos backwards, get position list in mt/cen for most massive halo
    #my.run("center_main.py "+str(nstart)+" "+str(nstop))
    # fixes missing snapshots by linear interpolation
    # interpolate missing snapshots by assuming last halo position
    #my.run("fix_shrink_spheres.py")
    my.run("fill_snap_xyz.py "+str(nstart)+" "+str(nstop))
    
if(action==10):
    print "> plot panels"
    f=1 # scaling of box wrt rvir
    calc = False; vis  = False;
    show = False; run  = False;
    loop = True;

    my.mkdir(simdir+"/ana")
    ddm  = simdir + "/ana/dm/";    my.mkdir(ddm)
    dgas = simdir + "/ana/gas/";   my.mkdir(dgas)
    dstar= simdir + "/ana/stars/"; my.mkdir(dstar)

    x,y,z,r,snap=mys.getxyzrsnap(nstart,nstop)
    xs,ys,zs,rs,snap=mys.getxyzrsnap_stars(nstart,nstop)
    print "done"
    for i in range(nstop-nstart+1):
        nc = nstop-i
        stri=str(nc).zfill(5)
        print "nc = ",nc
        d = mys.d(nc)
        # scale
        #[i]=0.002; rs[i]=0.002;

        sx = str(x[i]); sy = str(y[i]); sz = str(z[i]); sr = str(r[i])
        ssx= str(sx[i]);ssy= str(ys[i]);ssz= str(sz[i]);ssr= str(rs[i])
        bndry  = " "+sx+" "+sy+" "+sz+" "+sr+" "
        starbndrys =" -xs "+ssx+" -ys "+ssy+" -zs "+ssz+" -rs "+ssr+" "
        bndrys =" -xs "+sx+" -ys "+sy+" -zs "+sz+" -rs "+sr+" "
        bndryc =" -xc "+sx+" -yc "+sy+" -zc "+sz+" -rc "+sr+" "

        #bu: old bndry:
        #xmi= str(x[i]-f*r[i]); xma = str(x[i]+f*r[i])
        #ymi= str(y[i]-f*r[i]); yma = str(y[i]+f*r[i])
        #zmi= str(z[i]-f*r[i]); zma = str(z[i]+f*r[i])
        #for gas: zma = z;  for all else: -zma zma
        #bndry="-xmi "+xmi+" -xma "+xma+" -ymi "+ymi+" -yma "+yma+" -zmi "+zmi+" -zma "+zma

        cmd1 = "count_stars -inp "+d+bndryc+" >> star_counts.dat"
        #cmd2 = "plot_columns.py star_counts.dat 1 4"
        cmd2 = "calc_sfr.py star_counts.dat 0 4"
        my.threadif(cmd1,cmd2,calc,vis,show,run)

        # get DM particle positions for each halo
        cmd1 = "get_sphere_dm -inp "+d+bndryc+" > "+ddm+"dm_"+stri+".dat"
        cmd1+= "&& octreef "+ddm+"dm_"+stri+".dat > "+ddm+"rho_"+stri+".dat"
        cmd2 ="vis_part_proj_dm.py "+bndry+" "\
               +ddm+"dm_"+stri+".dat "+ddm+"dm_part_"+stri+".png"
        cmd2+="; vis_dm_contour.py "+bndry+" "\
               +ddm+"rho_"+stri+".dat "+ddm+"contour_"+stri+".png"
        my.threadif(cmd1,cmd2,calc,vis,show,run)


        # plot gas

        cmd1 ="amr2map -typ 1 -lma 17 -inp "+d+" -out "+dgas+"gas_boxall_"+stri
        cmd1 = cmd1 +".dat -dir z "+bndrys
        cmd2 = "map2img.py -l --colormap=hot "+dgas+"gas_boxall_"+stri+".dat "
        cmd2 = cmd2+"-o "+dgas+"gas_boxall_"+stri+".png"
        my.threadif(cmd1,cmd2,calc,vis,show,run)

        # plot pressure

        cmd1 = "amr2map -typ 5 -lma 17 -inp "+d
        cmd1 = cmd1+" -out "+dgas+"p_"+stri+".dat "+"-dir z "+bndrys
        cmd2 = "map2img.py -l --colormap=hot "+dgas+"p_"+stri+".dat "
        cmd2 = cmd2+" -o "+dgas+"p_"+stri+".png"
        my.threadif(cmd1,cmd2,calc,vis,show,run)

        # plot stars

        #cmd = "vis_parts.py stars.part"
        #cmd = "get_sphere_stars -inp "+d\
        #+" -xc "+sx+" -yc "+sy+" -zc "+sz+" -r "+sr+">ana/stars/stars_"+stri
        cmd1 = "metal2map -inp "+d+" -out "+dstar+"stars_"+stri+".dat "
        cmd1 = cmd1+"-nx 512 -ny 512 -dir z "+starbndrys+" -str true "
        cmd2 = "map2img.py -l --colormap=jet "+dstar+"stars_"+stri+".dat "
        cmd2 = cmd2+"-o "+dstar+"stars_"+stri+".png"
        my.threadif(cmd1,cmd2,calc,vis,show,run)

        if(not loop):
            break

if(action==11):
    for i in range(nstop-nstart+1):
        nc = nstop-i
        stri=str(nc).zfill(5)
        print "nc = ",nc
        d = mys.d(nc)
        nhalo=mys.get_nhalo(nc)

        for j in range(nhalo):
            # get DM particle positions for each halo
            cmd = "calc_rhalf_star.py "+stri+" "+str(j+1)
            my.thread(cmd)
        if(not loop):
            break
# gen_profs
if(action==13):
    print "            13> plot DM/gas/star profile"
    f=1
    ddm=simdir+"ana/dm/"
    x,y,z,r,snap=mys.getxyzrsnap(nstart,nstop)
    for i in range(nstop-nstart+1):
        nc = nstop-i
        print "nc = ",nc
        d = mys.d(nc)
        # scale
        r[i]=0.002

# TODO: use star positions
        sx = str(x[i]); sy = str(y[i]); sz = str(z[i]); sr = str(r[i])
        stri=str(nc).zfill(5)

        # get DM particle positions for each halo
#        cmd1 = "gen_prof_dm.py "+str(nc)+" "+sx+" "+sy+" "+sz+" "+sr;
        cmd1 = "gen_prof_stars.py "+str(nc)+" "+sx+" "+sy+" "+sz+" "+sr;
        cmd2 = "plot_prof_sph.py "+ddm+"prof_"+stri+".dat "\
               +ddm+"rho_"+stri+".dat "\
               +ddm+"prof_sph_"+stri+".png "+sx+" "+sy+" "+sz
        my.threadif(cmd1,cmd2,True,True,True,True)

if(action==14):
    print "generating wall"

    my.mkdir("ana/wall")
    show=False;run=False
    for ncounter in range(nstop-nstart+1):
        nc = ncounter + nstart
        num5 = str(nc).zfill(5)
        d  = simdir + "/output_"+num5+"/"

        my.run("head "+d+"/info_"+num5+".txt | tail -n1 | cut -d'=' -f2>"+d+"aexp")
        af = open(d+"aexp","r")
        a  = af.readline()

        # plot a wallpaper, with text
        cmd = "pngwall.py "
#        cmd = cmd + simdir+"ana/dm/dm_map_"+num5+".png "
        cmd = cmd + simdir+"ana/dm/contour_"+num5+".png "
        cmd = cmd + simdir+"ana/stars/stars_"+num5+".png "
        cmd = cmd + simdir+"ana/gas/gas_"+num5+".png "
        cmd = cmd + simdir+"ana/dm/dm_part_"+num5+".png "
        cmd = cmd + simdir+"ana/dm/prof_"+num5+".png "
        cmd = cmd + simdir+"ana/gas/p_"+num5+".png "
        cmd = cmd + simdir+"ana/wall/wall_"+num5+".png "
        cmd = cmd + num5 + " " + a
        if(show):print cmd
        if(run):my.thread(cmd)

        cmd = "gen_3step.py ana/dm/contour_"+num5+".png ana/gas/gas_"+num5+".png ana/dm/prof_sph_"+num5+".png ana/wall/3wall_"+num5+".png"
        if(True):print cmd
        if(True):my.thread(cmd)
        

if(action==15):
    print "> compile movies"
    men = "mencoder -ovc lavc -lavcopts vbitrate=5000 -mf type=png:fps=1 "
    my.thread( men +"-o ana/prof_new.avi mf://ana/dm/prof_*.png")
    my.thread( men +"-o ana/dm_new.avi mf://ana/dm/dm_part_*.png")
    my.thread( men +"-o ana/stars_new.avi mf://ana/stars/stars_*.png")
    my.thread( men +"-o ana/p_new.avi mf://ana/gas/p_*.png")
    my.thread( men +"-o ana/gas_new.avi mf://ana/gas/gas_*.png")
    my.thread( men +"-o ana/dm_map_new.avi mf://ana/dm/dm_map_*.png")
    my.thread( men +"-o ana/wall.avi mf://ana/wall/wall_*.png")

if(action==16):
    print "plot zoom in"
    calc = True;  vis=True;  show = True;  run  = False
    loop = False;
    
    my.mkdir(simdir+"/ana")
    in1 = open("mt/zoom")
    x=[]; y=[]; z=[]; r=[];
    for line in in1:
        val = line.split()
        x.append(float(val[0]))
        y.append(float(val[1]))
        z.append(float(val[2]))
        r.append(float(val[3]))
    in1.close()

    my.mkdir(simdir+"/ana/dm")
    my.mkdir(simdir+"/ana/gas")
    my.mkdir(simdir+"/ana/stars")

    print "len(x): ",len(x)
    for i in range(len(r)):
        nc = i
        d = simdir+"/output_"+str(nstart).zfill(5)+"/"
        sx = str(x[nc]); sy = str(y[nc]); sz = str(z[nc]); sr = str(r[nc])
        bndrys =" -xs "+sx+" -ys "+sy+" -zs "+sz+" -rs "+sr+" "
        stri =str(i)

        # plot gas

        cmd1 = "amr2map -typ 1 -lma 17 -inp "+d+" -out ana/gas/gas_"+stri
        cmd1 = cmd1 +".dat -dir z "+bndrys
        cmd2 = "map2img.py -l --colormap=hot ana/gas/gas_"+stri+".dat "
        cmd2 = cmd2+"-o ana/gas/gas_zoom_"+stri+".png"
        my.run(cmd1,cmd2,calc,vis,show,run)
        if(not loop):
            break

# plot phase diagram
if(action==18):
    for ncounter in range(nstop-nstart+1):
        nsnap = ncounter + nstart
        d  = mys.d(nsnap)
        pd = simdir+"phasediag/"
        my.mkdir(pd);

        fn = pd+"temp_"+str(nsnap).zfill(5)
        cmd = "rm "+fn+";get_temp -inp "+d+" -out "+fn
        cmd +=" -lma 21 -typ 18"
#        cmd ="plot_temp_rho.py "+fn+" "+fn+".png"
        print cmd

        if(brute):
            my.thread(cmd)
        else:
            my.run(cmd)
        if(not loop):
            break
