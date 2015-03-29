#!/usr/bin/env python2

## \file
# implement register with AHF related jobs

# (c) 2014 ETHZ, Pascal Steger, pascal@steger.aero

import wx
import pdb

import lib.mysql as mys
import lib.initialize as my
import get_halopos_rezoom

class TabAHF(wx.Panel):
    def __init__(self, parent):
        self.sim = mys.get_active_sim()
        """Tab with all things to setup for a new"""
        wx.Panel.__init__(self, parent=parent)
        self.SetBackgroundColour("darkred")
        sizer = wx.BoxSizer(wx.VERTICAL)
        # choose sim, nstart, nstop (same as in TabAnalysis)
        # run r2g
        btnR2G = wx.Button(self, label="ramses2gadget")
        sizer.Add(btnR2G, 0, wx.ALL, 3)
        self.Bind(wx.EVT_BUTTON, self.onRunR2G,btnR2G)
        # run AHF
        btnAHF = wx.Button(self, label="run AHF")
        sizer.Add(btnAHF, 0, wx.ALL, 3)
        self.Bind(wx.EVT_BUTTON, self.onRunAHF,btnAHF)
        # fill into DB
        btnFillAHF = wx.Button(self, label="fill AHF into DB")
        sizer.Add(btnFillAHF, 0, wx.ALL, 3)
        self.Bind(wx.EVT_BUTTON, self.onFillAHF, btnFillAHF)
        # show DM density and AHF halos
        btnHalos = wx.Button(self, label="show AHF halos")
        sizer.Add(btnHalos, 0, wx.ALL, 3)
        self.Bind(wx.EVT_BUTTON, self.onShowHalos,btnHalos)
        # get three most massive substructures
        #   for resimulation with zoom
        #   we want no structure with m>mvir/2 inside 5rvir
        btnGet3Halos = wx.Button(self, label="get 3 halos")
        sizer.Add(btnGet3Halos, 0, wx.ALL, 3)
        self.Bind(wx.EVT_BUTTON, self.onGet3Halos, btnGet3Halos)
        # generate spheres
        btnGenSpheres = wx.Button(self, label="gen spheres")
        sizer.Add(btnGenSpheres, 0, wx.ALL, 3)
        self.Bind(wx.EVT_BUTTON, self.onGenSpheres, btnGenSpheres)
        # center with shrinking spheres
        btnShrink = wx.Button(self, label="shrink spheres")
        sizer.Add(btnShrink, 0, wx.ALL, 3)
        self.Bind(wx.EVT_BUTTON, self.onShrink, btnShrink)
        # generate spheres with new center
        btnGenAOI = wx.Button(self, label="get sphere AOI")
        sizer.Add(btnGenAOI, 0, wx.ALL, 3)
        self.Bind(wx.EVT_BUTTON, self.onGenAOI, btnGenAOI)

        self.SetSizer(sizer)
    ## \fn __init__(self, parent)
    # fill tab for AHF control
    # @parent parent

    def onRunR2G(self, event):
        nstart,nstop=mys.get_range()
        print('nstart, nstop = ', nstart, nstop)
        for nc in range(nstop-nstart+1):
            nsnap = nc + nstart
            print('nsnap = ', nsnap)
            if(not mys.snap_exists(nsnap)):
                continue
            d = mys.d(nsnap)
            my.mkdir(d+"/r2g")
            if(mys.is_dmonly()):
                cmd = "mpirun -np 8 r2g -d "+d
            else:
                cmd = "mpirun -np 8 r2g -i "+d
            print(cmd)
            my.thread(cmd)
        self.SetBackgroundColour("darkgreen")
        my.done()
    ## \fn onRunR2G(self, event)
    # start ramses2gadget
    # @param event

    def onRunAHF(self,event):
        nstart,nstop=mys.get_range()
        print('nstart, nstop: ',nstart, nstop)
        for nc in range(nstop-nstart+1):
            nsnap = nc + nstart
            print('nsnap = ', nsnap)
            if(not mys.snap_exists(nsnap)):
                continue
            d  = mys.d(nsnap)
            f = open(d+"AHFinput","w")

            f.write("[AHF]\n\
\n\
# (stem of the) filename from which to read the data to be analysed\n\
ic_filename       = r2g/r2g.\n\
\n\
# what type of input file (cf. src/libio/io_file.h)\n\
ic_filetype       = 61\n\
\n\
# prefix for the output files\n\
outfile_prefix    = testahf\n\
\n\
# number of grid cells for the domain grid (1D)\n\
LgridDomain       = 4\n\
\n\
# number of grid cells for the domain grid (1D) (limits spatial resolution to BoxSize/LgridMax)\n\
LgridMax          = 16777216\n\
\n\
# refinement criterion on domain grid (#particles/cell)\n\
NperDomCell       = 5.0\n\
\n\
# refinement criterion on all higher resolution grids (#particles/cells)\n\
NperRefCell       = 5.0\n\
\n\
# particles with velocity v > VescTune x Vesc are considered unbound\n\
VescTune          = 1.5\n\
\n\
# minimum number of particles for a halo\n\
NminPerHalo       = 100\n\
\n\
# normalisation for densities (1: RhoBack(z), 0:RhoCrit(z))\n\
RhoVir            = 1\n\
\n\
# virial overdensity criterion (<0: let AHF calculate it); Rvir is defined via M(<Rvir)/Vol = Dvir * RhoVir\n\
Dvir              = 200\n\
#-1\n\
\n\
# maximum radius (in Mpc/h) used when gathering initial set of particles for each halo (should be larger than the largest halo expected)\n\
MaxGatherRad      = 0.2\n\
\n\
# the level on which to perform the domain decomposition (MPI only, 4=16^3, 5=32^3, 6=64^3, 7=128^3, 8=256^3, etc.)\n\
LevelDomainDecomp = 4\n\
\n\
# how many CPU's for reading (MPI only)\n\
NcpuReading       = 16\n\
\n\
############################### FILE SPECIFIC DEFINITIONS ###############################\n\
\n\
# NOTE: all these factors are supposed to transform your internal units to\n\
#           [x] = Mpc/h\n\
#           [v] = km/sec\n\
#           [m] = Msun/h\n\
#           [e] = (km/sec)^2\n\
\n\
[GADGET]\n\
GADGET_LUNIT      = 0.001\n\
GADGET_MUNIT      = 1.0E10\n")
            f.close()

            print("TODO: AHFinput written?")
            cmd = "cd "+d+" && AHF AHFinput"
            cmd += " && mv "+d+"*_halos "+d+"halos"
            cmd += " && mv "+d+"*_centres "+"centres"
            cmd += " && mv "+d+"*_particles "+d+"particles"
            cmd += " && mv "+d+"*_profiles "+d+"profiles"
            cmd += " && rm -f "+d+"halo"
            cmd += " && echo 'done'"
            print(cmd)
            my.thread(cmd)
        my.done()
    ## \fn onRunAHF(self, event)
    # start Amiga Halo Finder
    # @param event

    def onFillAHF(self, event):
        nstart,nstop=mys.get_range()
        for nc in range(nstop-nstart+1):
            nsnap = nc + nstart
            if(not mys.snap_exists(nsnap)):
                continue
            cmd = "fill_ahf_to_db.py " + str(nsnap)
            print(cmd)
            my.thread(cmd)
            exit
        my.done()
    ## \fn onFillAHF(self, event)
    # fill AHF output to MySQL database
    # @param event

    def onShowHalos(self, event):
        nstart,nstop = mys.get_range()
        nstart = nstop #! show only latest snapshot
        for nc in range(nstop-nstart+1):
            nsnap = nc + nstart
            if(not mys.snap_exists(nsnap)):
                continue
            cmd = "vis_AHF_centers.py "+str(nsnap)+" && xview "+mys.d(nsnap)+"ov.dat.png"
            print(cmd)
            my.thread(cmd)
            exit
    ## \fn onShowHalos(self, event)
    # plot AHF halo positions
    # @param event

    def onGet3Halos(self, event):
        get_halopos_rezoom.snap(mys.get_nstop())
        my.done()
    ## \fn onGet3Halos(self, event)
    # output positions of 3 most massive halos in simulation
    # @param event

    def onGenSpheres(self, event):
        nstart,nstop=mys.get_range()
        for nc in range(nstop-nstart+1):
            nsnap = nc + nstart
            if(not mys.snap_exists(nsnap)):
                continue
            d  = mys.d(nsnap)
            my.mkdir(d+"dm")
            if(not mys.is_dmonly()):
                my.mkdir(d+"stars");
            cmd = "cd "+d+" && gen_spheres.py "+str(nsnap)+" 1"
            my.run(cmd)
        my.done()
        #TODO: remove empty octree output files
    ## \fn onGenSpheres(self, event)
    # find particles inside virial radius
    # @param event

    def onShrink(self, event):
        nstart,nstop=mys.get_range()
        # have xs,xs_star set in the end
        for nc in range(nstop-nstart+1):
            nsnap = nc + nstart
            if(not mys.snap_exists(nsnap)):
                continue
            print(mys.get_nhalo(nsnap))
            for i in range(3):#range(mys.get_nhalo(nsnap)):
                if(not mys.exists_snap(nsnap)):
                    print("snapshot "+str(nsnap)+" missing")
                    continue
                cmd = "shrink_sphere.py "+str(nsnap)+" "+str(i+1)
                my.run(cmd)
                if(not mys.is_dmonly()):
                    cmd = "shrink_sphere_stars.py "+str(nsnap)+" "+str(i+1)
                    my.run(cmd)
                # TODO: if empty xs (shrinking sphere algorithm did not succeed: use xc)
        my.done()
    ## \fn onShrink(self, event)
    # shrinking sphere algorithm to find exact center
    # @param event

    def onGenAOI(self, event):
        nstart,nstop=mys.get_range()
        for nc in range(nstop-nstart+1):
            nsnap = nc + nstart
            if(not mys.snap_exists(nsnap)):
                continue
            d  = mys.d(nsnap)
            my.mkdir(d+"stars")
            my.mkdir(d+"dm")
            cmd = "gen_spheres.py "+str(nsnap)+" 2" # 1 at end: centers from AHF. 2: after shrinking sphere
            my.run(cmd)
        my.done()
    ## \fn onGenAOI(self, event)
    # generate area of interest
    # @param event
## \class TabAHF(wx.Panel)
# populate tab for running Amiga Halo Finder
# @param wx.Panel where to fill it
