#!/usr/bin/env python2

## \file
# implement MySQL related commands
# needs a MySQL server running on (another) computer, see act/lib/initialize for details

# (c) 2015 ETHZ, Pascal Steger, pascal@steger.aero

import os
import wx
import pdb

import lib.mysql as mys
import lib.initialize as my

class TabMySQL(wx.Panel):
    sim = "800nec"
    def __init__(self, parent):
        """Tab with all things to setup for a new simulation"""
        wx.Panel.__init__(self, parent=parent)
        self.SetBackgroundColour("darkred")

        sizer = wx.BoxSizer(wx.VERTICAL)
        grid = wx.GridBagSizer(hgap=5, vgap=5)
        # clear, so we can delete an old or precursor simulation
        btnClear = wx.Button(self, label="Clear all")
        sizer.Add(btnClear, 0, wx.ALL, 3)
        self.Bind(wx.EVT_BUTTON, self.onClear, btnClear)
        # setup db
        btnSetup = wx.Button(self, label="Setup DB")
        sizer.Add(btnSetup, 0, wx.ALL, 3)
        self.Bind(wx.EVT_BUTTON, self.onSetup, btnSetup)
        # setup simulation
        self.btnSimSetup = wx.Button(self, label="Setup Sim")
        grid.Add(self.btnSimSetup, pos=(0,0))
        self.Bind(wx.EVT_BUTTON, self.onSimSetup, self.btnSimSetup)
        print(mys.get_sims())
        self.sampleList = mys.get_sims()
        self.edithear = wx.ComboBox(self, size=(160, -1), choices=self.sampleList, style=wx.CB_DROPDOWN)
        grid.Add(self.edithear, pos=(0,1))
        self.Bind(wx.EVT_COMBOBOX, self.ChooseSim, self.edithear)
        self.Bind(wx.EVT_TEXT, self.ChooseSimNew,self.edithear)
        self.dmonly = wx.RadioButton ( self, -1, 'DM only', style = wx.RB_GROUP )
        self.hydro  = wx.RadioButton ( self, -1, 'hydro' )
        grid.Add ( self.dmonly, ( 0, 2 ) )
        grid.Add ( self.hydro, ( 0, 3 ) )
        self.Bind(wx.EVT_RADIOBUTTON, self.onDM,    self.dmonly )
        self.Bind(wx.EVT_RADIOBUTTON, self.onHydro, self.hydro )
        self.onDM(False)
        sizer.Add(grid, 0, wx.ALL, 3)
        # delete all from DB for a given sim
        self.btnSimDelete = wx.Button(self, label="Delete Sim")
        sizer.Add(self.btnSimDelete, 0, wx.ALL, 3)
        self.Bind(wx.EVT_BUTTON, self.onSimDelete, self.btnSimDelete)
        # folder structure
        btnFolder = wx.Button(self, label="Folder")
        sizer.Add(btnFolder, 0, wx.ALL, 3)
        self.Bind(wx.EVT_BUTTON, self.onFolderSetup, btnFolder)
        # folder structure
        btnMasses = wx.Button(self, label="Find Masses")
        sizer.Add(btnMasses, 0, wx.ALL, 3)
        self.Bind(wx.EVT_BUTTON, self.onMass, btnMasses)
        self.SetSizer(sizer)
    ## \fn __init__(self, parent)
    # constructor
    # @param parent window

    def onClear(self, event):
        mys.clear()
        print("cleared DB.")
        my.done()
    ## \fn onClear(self, event)
    # clear database
    # @param event

    def ChooseSimNew(self, event):
        self.sim = event.GetString()
        print("set sim to "+self.sim)
        my.done()
    ## \fn ChooseSimNew(self, event)
    # set simulation to chosen sim
    # @param event

    def ChooseSim(self, event):
        self.sim = event.GetString()
        print("set sim to "+self.sim)
    ## \fn ChooseSim(self, event)
    # set simulation string to newly entered string
    # @param event

    def onSetup(self, event):
        mys.setup()
        my.done()
    ## \onSetup(self, event)
    # setup mySQL connection, database
    # @param event

    def onDM(self, event):
        print(" ")
        print("simulation is dm only")
        self.dmonly=True
        my.done()
    ## \fn onDM(self, event)
    # set dmonly for DMonly simulations
    # @param event

    def onHydro(self, event):
        print(" ")
        print("simulation is hydro")
        self.dmonly=False
        my.done()
    ## \fn onHydro(self, event)
    # set hydro flag for simulation
    # @param event

    def onSimSetup(self, event):
        print(" ")
        print("setup DB with "+self.sim)
        # get simdir, max output
        self.simdir="/library/home/psteger/sci/sim/"+self.sim
        dirlist=os.listdir(self.simdir)
        mx=0
        for fname in dirlist:
            if(fname[:7]=='output_'):
                mx=max(int(fname[7:]),mx)
        self.nsnap=mx
        print(mx)
        mys.fill_sim(self.sim, self.nsnap, self.dmonly)
        mys.print_sim(self.sim)
        my.done()
    ## \fn onSimSetup(self, event)
    # store parameters in database for current simulation
    # @param event

    def onSimDelete(self, event):
        print(" ")
        print("delete DB with "+self.sim)
        my.sql('delete from halo where sim="'+self.sim+'"');
        my.sql('delete from sim where name="'+self.sim+'"');
        my.sql('delete from snapshot where sim="'+self.sim+'"');
        my.done()
    ## \fn onSimDelete(self, event)
    # purge simulation completely from DB
    # @param event

    def onMass(self, event):
        print(" ")
        print("find DM mass and star mass (if present)")
        nstop = mys.get_nstop()
        d = mys.d(nstop)
        print("m_dm = ")
        my.run_command("get_sphere_dm -inp "+d+" -xc 0.5 -yc 0.5 -zc 0.5 -rc 0.001|head -n1|cut -d'.' -f1,2|rev|cut -d' ' -f2-|rev")
        print("m_star = ")
        my.run_command("get_sphere_stars -inp "+d+" -xc 0.5 -yc 0.5 -zc 0.5 -rc 0.01|head -n1|cut -d'.' -f1,2|rev|cut -d' ' -f2-|rev")
        my.done()
    ## \fn onMass(self, event)
    # find DM mass and star mass
    # @param event

    def onFolderSetup(self, event):
        print(" ")
        print("folder setup")
        nstop = mys.getNmax()
        print('nstop = ',nstop)
        for nc in range(nstop):
            if(not mys.snap_exists(nc+1)):
                continue
            print("generating folders in output ",nc+1)
            d=mys.d(nc+1)
            print(d)
            my.mkdir(d+'amr');    my.mv(d+"amr_*",d+"amr")
            my.mkdir(d+'grav');   my.mv(d+"grav_*",d+'/grav')
            my.mkdir(d+'part');   my.mv(d+'part_*',d+'/part')
            if(not self.dmonly):
                print("create folder hydro as well")
                my.mkdir(d+'hydro');  my.mv(d+"hydro_*",d+'/hydro')
        self.SetBackgroundColour("darkgreen")
        my.done()
    ## \fn onFolderSetup(self, event)
    # prepare folder structure for further analysis
    # @param event

## \fn class TabMySQL
# setup MySQL database in this tab
