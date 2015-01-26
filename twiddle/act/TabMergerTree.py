#!/usr/bin/env python2

## \file
# implement all commands used for finding a merger tree

# (c) 2014 ETHZ, Pascal Steger, pascal@steger.aero


import wx
import lib.mysql as mys
import lib.initialize as my

class TabMergerTree(wx.Panel):
    def __init__(self, parent):
        """Tab with all things to setup for a new"""
        wx.Panel.__init__(self, parent=parent)
        self.SetBackgroundColour("darkred")

        sizer = wx.BoxSizer(wx.VERTICAL)
        grid = wx.GridBagSizer(hgap=5, vgap=5)

        # particle list
        btnPlist = wx.Button(self, label="particle list")
        sizer.Add(btnPlist, 0, wx.ALL, 3)
        self.Bind(wx.EVT_BUTTON, self.onPlist,btnPlist)

        # merger tree
        btnMtree = wx.Button(self, label="merger tree")
        sizer.Add(btnMtree, 0, wx.ALL, 3)
        self.Bind(wx.EVT_BUTTON, self.onMtree, btnMtree)

        # follow boss
        btnFboss = wx.Button(self, label="follow main halo")
        sizer.Add(btnFboss, 0, wx.ALL, 3)
        self.Bind(wx.EVT_BUTTON, self.onFboss, btnFboss)

        self.SetSizer(sizer)
    ## \fn TabMergerTree(wx.Panel)
    # populate tab for run of merger tree
    # @param wx.Panel where to draw to

    def onPlist(self, event):
        nstart, nstop = mys.get_range()
        for nc in range(nstop-nstart+1):
            nsnap = nc + nstart
            # for all halos:
            cmd = "gen_particles_list.py "+str(nsnap)
            cmd += "&& gen_particles_dm.py "+str(nsnap)
            my.run(cmd)
            # if(not loop):
            # break
        my.done()
    ## \fn onPlist(self, event)
    # list particles inside snapshot span
    # @param event click

    def onMtree(self, event):
        nstart,nstop=mys.get_range()
        if(nstart==nstop):
            ns = str(nstart)
            cmd = "substructure.py "+ns
            my.thread(cmd)
        else:
            for nc in range(nstop-nstart+1):
                cmd = "merger_tree.py "+str(nstop-nc)+" "+str(nstop-1-nc)
                cmd+= "&& substructure.py "+str(nstop-nc)
                my.thread(cmd)
        my.done()
    ## \fn onMtree(self, event)
    # find merger tree
    # @param event click

    def onFboss(self, event):
        nstart,nstop=mys.get_range()
        # follow halos backwards, get position list in mt/cen for most massive halo
        #my.run("center_main.py "+str(nstart)+" "+str(nstop))
        # fixes missing snapshots by linear interpolation
        # interpolate missing snapshots by assuming last halo position
        #my.run("fix_shrink_spheres.py")
        my.run("fill_snap_xyz.py "+str(nstart)+" "+str(nstop))
        my.done()
    ## \fn onFboss(self, event)
    # follow most massive halo backwards in time
    # @param event click
## \class TabMergerTree(wx.Panel)
# class with functionality to establish a merger tree
# @param wx.Panel
