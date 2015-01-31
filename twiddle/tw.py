#!/usr/bin/env python2.7

## \file
# GUI for all tools needed to work with ramses outputs
# TODO: check out pymses

# (c) 2015 ETHZ, Pascal Steger, pascal@steger.aero

import numpy as np
import wx

import act.TabMySQL as tms
import act.TabAHF as tah
import act.TabMergerTree as tmt
import act.TabAnalysis as tas
import act.lib.mysql as mys

class MainFrame(wx.Frame):
    def __init__(self):
        mys.setup_sim()
        mys.setup_halo()
        mys.setup_snapshot()
        wx.Frame.__init__(self, None, wx.ID_ANY, "twiddle with cosmological simulations", size=(600,400))
        panelnb = wx.Panel(self,-1)
        nb = wx.Notebook(panelnb)
#1
        tabAnalysis = tas.TabAnalysis(nb)
        nb.AddPage(tabAnalysis, "Analysis")
#2
        tabMtree = tmt.TabMergerTree(nb)
        nb.AddPage(tabMtree, "MergerTree")
#3
        tabAHF = tah.TabAHF(nb)
        nb.AddPage(tabAHF, "AHF")
#4
        tabMySQL = tms.TabMySQL(nb)
        nb.AddPage(tabMySQL, "MySQL")
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(nb, 1, wx.ALL|wx.EXPAND, 0)
        panelnb.SetSizer(sizer)
        #general things first
        panelsim = wx.Panel(self,-1)
        hbox = wx.BoxSizer(wx.VERTICAL)
        self.simList = np.array(mys.get_sims())
        sim = mys.get_active_sim()
        self.edithear = wx.ComboBox(self, choices=self.simList, value=sim, style=wx.CB_DROPDOWN)
        vbox1 = wx.BoxSizer(wx.HORIZONTAL)
        vbox1.Add(self.edithear,1,wx.ALL)
        hbox.Add(vbox1,1,wx.EXPAND)
        self.Bind(wx.EVT_COMBOBOX, self.ChooseSim, self.edithear)
        # choose snapshot span
        nstart = mys.get_nstart()
        nstop = mys.get_nstop()
        self.start = wx.TextCtrl(self, value=str(nstart), size=(-1,-1))
        self.stop = wx.TextCtrl(self, value=str(nstop), size=(-1,-1))
        vbox = wx.BoxSizer(wx.HORIZONTAL)
        vbox.Add(self.start, 1,wx.ALL)
        vbox.Add(self.stop, 1,wx.ALL)
        self.Bind(wx.EVT_TEXT, self.setStart, self.start)
        self.Bind(wx.EVT_TEXT, self.setStop, self.stop)
        hbox.Add(vbox,1,wx.EXPAND)
        vemp = wx.BoxSizer(wx.VERTICAL)
        hbox.Add(vemp,10,wx.EXPAND)
        panelsim.SetSizer(hbox)
        box = wx.BoxSizer(wx.HORIZONTAL)
        box.Add(panelsim,1,wx.EXPAND)
        box.Add(panelnb,3,wx.EXPAND)
        self.SetSizer(box)
        self.Layout()
        self.Show(True)

    def ChooseSim(self, event):
        sim = event.GetString()
        # print( 'set sim to:',sim)
        mys.set_active_sim(sim)
        #sm = mys.get_active_sim()
        # print('active sim:',sim)
        nstart = mys.get_nstart()
        nstop  = mys.get_nstop()
        self.start.SetValue(str(nstart))
        self.stop.SetValue(str(nstop))
    ## \fn ChooseSim(self, event)
    # pick a simulation from name
    # @param event

    def setStart(self, event):
        nstart = int(event.GetString())
        # print("chosen nstart: "+str(nstart))
        mys.set_nstart(nstart)
    ## \fn setStart(self, event)
    # set starting output number
    # @param event

    def setStop(self, event):
        nstop = int(event.GetString())
        # print("chosen nstart: "+str(nstop))
        mys.set_nstop(nstop)
    ## \fn setStop(self, event)
    # set stopping snapshot
    # @param event

## \class MainFrame()
# Frame that holds all other widgets

class MyApp(wx.App):
    def OnInit(self):
        frame = MainFrame()
        frame.Show(True)
        self.SetTopWindow(frame)
        return True
    ## \fn OnInit(self)
    # start window

## \class MyApp()
# basic application

def main():
    app = MyApp(False)
    app.MainLoop()
## \fn main()
# start application and wait for user input

if __name__ == "__main__":
    main()
