#!/usr/bin/env python2

## \file
# implement all analysis routines
# all done via external tools called from the command line

# (c) 2015 ETHZ, Pascal Steger, pascal@steger.aero

import wx
import lib.mysql as mys
import lib.initialize as my

amr2map="amr2map_8byte"

class TabAnalysis(wx.Panel):
    def __init__(self, parent):
        self.radio=0; self.movie=0; self.loop=1; self.calc=1; self.fix=0
        wx.Panel.__init__(self, parent=parent)
        #self.SetBackgroundColour("gray")
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        grid      = wx.GridBagSizer(hgap=5, vgap=5)
        hsiz1     = wx.BoxSizer(wx.HORIZONTAL)
        hsiz2     = wx.BoxSizer(wx.HORIZONTAL)
        radioList = [
            'density',    'pressure',  'temperature',
            'DM contour', 'metal',     'phase diag',
            'rho_dm',     'rho_gas',   'rho_star',
            'SFR']
        rb = wx.RadioBox(self, label="Output:",
                         choices=radioList, majorDimension=3,
                         style=wx.RA_SPECIFY_COLS)
        grid.Add(rb, pos=(4,0), span=(1,2))
        self.Bind(wx.EVT_RADIOBOX, self.EvtRadioBox, rb)
        hsiz1.Add(grid, 0, wx.ALL, 5)

        #calc things, not only visualize
        self.hagcalc = wx.CheckBox(self, label="calc")
        hsiz2.Add(self.hagcalc, 2, wx.ALL, 2)
        self.Bind(wx.EVT_CHECKBOX, self.checkCalc, self.hagcalc)
        self.hagcalc.SetValue(1)
        # compile plots into movie?
        self.hagmovie = wx.CheckBox(self, label="movie")
        hsiz2.Add(self.hagmovie, 2, wx.ALL, 2)
        self.Bind(wx.EVT_CHECKBOX, self.checkMovie, self.hagmovie)
        #set fix comoving scale, do not use virial radius
        self.hagfix = wx.CheckBox(self, label="fix r")
        hsiz2.Add(self.hagfix, 2, wx.ALL, 2)
        self.Bind(wx.EVT_CHECKBOX, self.checkFix, self.hagfix)
        # loop?
        self.hagloop = wx.CheckBox(self, label="loop")
        hsiz2.Add(self.hagloop, 2, wx.ALL, 2)
        self.Bind(wx.EVT_CHECKBOX, self.checkLoop, self.hagloop)
        self.hagloop.SetValue(1)
        self.buttonPlot =wx.Button(self, label="Plot!")
        self.Bind(wx.EVT_BUTTON, self.OnPlot,self.buttonPlot)
        hsiz2.Add(self.buttonPlot, 2, wx.ALL, 2)
        mainSizer.Add(hsiz1, 0, wx.ALL, 5); mainSizer.Add(hsiz2, 0, wx.ALL, 5)
        # maximal gas level
        lmaval = mys.get_lma()
        self.lma = wx.TextCtrl(self, value=str(lmaval), size=(-1,-1))
        vbox = wx.BoxSizer(wx.HORIZONTAL)
        vbox.Add(self.lma, 1,wx.ALL)
        self.Bind(wx.EVT_TEXT, self.setLMA, self.lma)
        mainSizer.Add(vbox, 0, wx.ALL, 5)
        self.SetSizerAndFit(mainSizer)
    ## \fn __init__(self, parent)
    # populate window
    # @param event

    def EvtRadioBox(self, event):
        print('plot: %d\n' % event.GetInt())
        self.radio=event.GetInt()
    ## \fn EvtRadioBox(self, event)
    # radio button for event
    # @param event

    def OnPlot(self, event):
        print("plot")
        f=1 # scaling of sphere wrt rvir
        vis  = True; show = True; run  = True;
        my.mkdir(mys.simdir()+"/ana")
        ddm  = mys.simdir() + "/ana/dm/";    my.mkdir(ddm)
        dgas = mys.simdir() + "/ana/gas/";   my.mkdir(dgas)
        dstar= mys.simdir() + "/ana/stars/"; my.mkdir(dstar)
        dpd  = mys.simdir() + "/phasediag/"; my.mkdir(dpd)
        nstart,nstop=mys.get_range()
        x,y,z,r,snap=mys.getxyzrsnap(nstop)#TODO: nstart,nstop
        xs,ys,zs,rs,snap=mys.getxyzrsnap_stars(nstart,nstop)
        #x,y,z,r,snap=mys.mt_xyzrsnap(nstart,nstop)
        #xs,ys,zs,rs,snap=mys.mt_xyzrsnap_stars(nstart,nstop)
        for i in range(nstop-nstart+1):
            nc = nstop-i
            if(not mys.snap_exists(nc)):
                continue
            stri=str(nc).zfill(5)
            print("nc = ",nc)
            print(x[i], xs[i])
            print(y[i], ys[i])
            print(z[i], zs[i])
            print(r[i], rs[i])
            d = mys.d(nc)
            # scale
            if self.fix:
                r[i]=0.002; rs[i]=0.002;

            sx = str(x[i])
            sy = str(y[i])
            sz = str(z[i])
            sr = str(r[i])
            ssx= str(sx[i])
            ssy= str(ys[i])
            ssz= str(sz[i])
            ssr= str(rs[i])
            lma = str(mys.get_lma())
            bndry  = " "+sx+" "+sy+" "+sz+" "+sr+" "
            bndryc =" -xc "+sx+" -yc "+sy+" -zc "+sz+" -rc "+sr+" "

            #calc all pix, dep. on which option was clicked
            if self.radio==0:
                print("gas density")
                ofname = dgas+"gas_boxall_"+stri+".png"
                cmd1 = amr2map+" -typ 1 -lma "+lma+" -inp "+d+" -out "+dgas+"gas_boxall_"+stri
                cmd1 = cmd1 +".dat -dir z "+bndryc
                cmd2 = "map2img.py -l --colormap=hot "+dgas+"gas_boxall_"+stri+".dat "
                cmd2 = cmd2+"-o "+ofname+" && xview "+ofname
                print(dgas+"gas_boxall_"+stri+".png")
                my.threadif(cmd1,cmd2,self.calc,vis,show,run)

            if self.radio==1:
                print("gas pressure")
                if(mys.is_dmonly()):
                    print("not available, dm only simulation!")
                    exit
                cmd1 = amr2map + " -typ 5 -lma "+lma+" -inp "+d
                cmd1 = cmd1+" -out "+dgas+"p_"+stri+".dat "+"-dir z "+bndryc
                cmd2 = "map2img.py -l --colormap=hot "+dgas+"p_"+stri+".dat "
                ofname = dgas+"p_"+stri+".png"
                cmd2 = cmd2+" -o "+ofname+" && xview "+ofname
                my.threadif(cmd1,cmd2,self.calc,vis,show,run)

            if self.radio==2:
                print("gas temperature")
                cmd1 = amr2map+" -typ 18 -lma "+lma+" -inp "+d
                cmd1 = cmd1+" -out "+dgas+"temp_"+stri+".dat "+"-dir z "+bndryc
                cmd2 = "map2img.py -l --colormap=hot "+dgas+"temp_"+stri+".dat "
                ofname = dgas+"temp_"+stri+".png"
                cmd2 = cmd2+" -o "+ofname+" && xview "+ofname
                my.threadif(cmd1,cmd2,self.calc,vis,show,run)

            if self.radio==3:
                print("DM contours")
                # get DM particle positions for each halo
                cmd1 = "get_sphere_dm -inp "+d+bndryc+" > "+ddm+"dm_"+stri+".dat"
                cmd1+= "&& octreef "+ddm+"dm_"+stri+".dat > "+ddm+"rho_"+stri+".dat"
                cmd2 ="vis_part_proj_dm.py "+bndry+" "\
                       +ddm+"dm_"+stri+".dat "+ddm+"dm_part_"+stri+".png"
                ofname=ddm+"rho_"+stri+".dat "+ddm+"contour_"+stri+".png"
                cmd2+="; vis_dm_contour.py "+bndry+" "\
                       +ofname+" && xview "+ofname
                my.threadif(cmd1,cmd2,self.calc,vis,show,run)

            if self.radio==4:
                print("metallicity")
                #cmd = "vis_parts.py stars.part"
                #cmd = "get_sphere_stars -inp "+d\
                #+" -xc "+sx+" -yc "+sy+" -zc "+sz+" -r "+sr+">ana/stars/stars_"+stri
                cmd1 = "metal2map -inp "+d+" -out "+dstar+"stars_"+stri+".dat "
                cmd1 = cmd1+"-nx 512 -ny 512 -dir z "+bndryc
                cmd2 = "map2img.py -l --colormap=jet "+dstar+"stars_"+stri+".dat "
                ofname = dstar+"stars_"+stri+".png"
                cmd2 = cmd2+"-o "+ofname+" && xview "+ofname
                my.threadif(cmd1,cmd2,self.calc,vis,show,run)

            if self.radio==5:
                print("phase diagram")
                fn = dpd+"temp_"+str(nc).zfill(5)
                cmd1 = "get_temp -inp "+d+" -out "+fn
                cmd1+=  " -lma "+lma+" -typ 18"
                ofname=fn+".png"
                cmd2 = "plot_temp_rho.py "+fn+" "+ofname+" && xview "+ofname
                my.threadif(cmd1,cmd2,self.calc,vis,show,run)

            if self.radio==6:
                print("rho_dm(r)")
                cmd1 = "gen_prof_dm.py "+str(nc)+" "+sx+" "+sy+" "+sz+" "+sr;
                ofname = ddm+"prof_"+stri+".png "
                cmd2 = "plot_prof_sph.py "+ddm+"prof_"+stri+".dat "\
                       +ddm+"rho_"+stri+".dat "\
                       +ofname+sx+" "+sy+" "+sz+" "+str(mys.get_z(nc))+"&& xview "+ofname
                my.threadif(cmd1,cmd2,self.calc,vis,show,run)


            if self.radio==7:
                print("rho_gas(r)")
                cmd1 = "gen_prof_gas.py "+str(nc)+" "+sx+" "+sy+" "+sz+" "+sr+" "+lma;
                ofname = dgas+"prof_"+stri+".png "
                cmd2 = "plot_prof_sph.py "+dgas+"prof_"+stri+".dat "\
                       +dgas+"rho_"+stri+".dat "\
                       +ofname+sx+" "+sy+" "+sz+" && xview "+ofname
                my.threadif(cmd1,cmd2,self.calc,vis,show,run)


            if self.radio==8:
                print("rho_star(r)")
                cmd1 = "gen_prof_stars.py "+str(nc)+" "+sx+" "+sy+" "+sz+" "+sr;
                ofname = dstar+"prof_"+stri+".png "
                cmd2 = "plot_prof_sph.py "+dstar+"prof_"+stri+".dat "\
                       +dstar+"rho_"+stri+".dat "\
                       +ofname+sx+" "+sy+" "+sz+" && xview "+ofname
                my.threadif(cmd1,cmd2,self.calc,vis,show,run)

            if self.radio==9:
                print("SFR")
                cmd1 = "count_stars -inp "+d+bndryc+" >> "+dstar+"star_counts.dat"
                #cmd2 = "plot_columns.py star_counts.dat 1 4"
                cmd2 = "calc_sfr.py "+dstar+"star_counts.dat 0 4"
                my.threadif(cmd1,cmd2,self.calc,vis,show,run)

            if(not self.loop):
                break
        my.done()
        #TODO: if movie, concat and save as avi
    ## \fn OnPlot(self, event)
    # call plotting routine based on radio button selection
    # @param event

    def checkMovie(self, event):
        self.movie=self.hagmovie.GetValue()
    ## \fn checkMovie(self, event)
    # get value of movie
    # @param event

    def checkLoop(self, event):
        self.loop=self.hagloop.GetValue()
    ## \fn checkLoop(self, event)
    # get value of loop
    # @param event

    def checkFix(self, event):
        self.fix=self.hagfix.GetValue()
    ## \fn checkFix(self, event)
    # get value of fix
    # @param event

    def checkCalc(self, event):
        self.calc=self.hagcalc.GetValue()
    ## \fn checkCalc(self, event)
    # get value of calc, whether to calculate anew
    # @param event

    def setLMA(self, event):
        lma = int(event.GetString())
        mys.set_lma(lma)
        print('maximum level lma = ', lma)
    ## \fn setLMA(self, event)
    # get value of maximum level to refine plot to
    # @param eventt
## \class TabAnalysis(wx.Panel)
# tab for starting analysis
# @param wx.Panel where to draw to
