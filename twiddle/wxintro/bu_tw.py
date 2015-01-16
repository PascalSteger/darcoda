#!/usr/bin/python
import wx

class panel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)

        mainSizer = wx.BoxSizer(wx.VERTICAL)
        grid = wx.GridBagSizer(hgap=5, vgap=5)
        hsiz1 = wx.BoxSizer(wx.HORIZONTAL)
        hsiz2  = wx.BoxSizer(wx.HORIZONTAL)

        self.quote = wx.StaticText(self, label="Twiddle with a snapshot.")
        grid.Add(self.quote, pos=(0,0))

        # status message panel
        self.logger = wx.TextCtrl(self, size=(200,300), style=wx.TE_MULTILINE | wx.TE_READONLY)
        self.buttonPlot =wx.Button(self, label="Plot!")
        self.Bind(wx.EVT_BUTTON, self.OnPlot,self.buttonPlot)
        self.buttonExit =wx.Button(self, label="Exit")
        self.Bind(wx.EVT_BUTTON, self.OnExit,self.buttonExit)

        # choose snapshot span
        self.lblname = wx.StaticText(self, label="snapshots:")
        grid.Add(self.lblname, pos=(1,0))

        self.start = wx.TextCtrl(self, value="1", size=(40,-1))
        grid.Add(self.start, pos=(1,1))
        self.Bind(wx.EVT_TEXT, self.EvtText, self.start)
        self.Bind(wx.EVT_CHAR, self.EvtChar, self.start)
        
        self.stop = wx.TextCtrl(self, value="270", size=(40,-1))
        grid.Add(self.stop, pos=(1,2))
        self.Bind(wx.EVT_TEXT, self.EvtText, self.stop)
        self.Bind(wx.EVT_CHAR, self.EvtChar, self.stop)

        # choose simulation
        self.sampleList = ['sim_aboley', 'nec_20111220']
        self.lblhear = wx.StaticText(self, label="simulation:")
        grid.Add(self.lblhear, pos=(3,0))
        self.edithear = wx.ComboBox(self, size=(95, -1), choices=self.sampleList, style=wx.CB_DROPDOWN)
        grid.Add(self.edithear, pos=(3,1))
        self.Bind(wx.EVT_COMBOBOX, self.ChooseSim, self.edithear)
        self.Bind(wx.EVT_TEXT, self.EvtText,self.edithear)

        grid.Add((10, 40), pos=(2,0))         # add a spacer to the sizer
         # compile plots into movie?
        self.movie = wx.CheckBox(self, label="movie")
        grid.Add(self.movie, pos=(4,0), span=(1,2), flag=wx.BOTTOM, border=5)
        self.Bind(wx.EVT_CHECKBOX, self.EvtCheckBox, self.movie)

        # what to plot
        radioList = ['star formation', 'DM contour', 'rho_gas',
                     'p_gas', 'metal', 'rho_rad prof', 'phase diag']
        rb = wx.RadioBox(self, label="Output:", pos=(20, 210),
                         choices=radioList, majorDimension=3, style=wx.RA_SPECIFY_COLS)
        grid.Add(rb, pos=(5,0), span=(1,2))
        self.Bind(wx.EVT_RADIOBOX, self.EvtRadioBox, rb)

        hsiz1.Add(grid, 0, wx.ALL, 5)
        hsiz1.Add(self.logger)
        mainSizer.Add(hsiz1, 0, wx.ALL, 5)
        #mainSizer.Add(self.buttonPlot, 0, wx.CENTER)
        hsiz2.Add(self.buttonPlot, 2, wx.ALL, 2)
        hsiz2.Add(self.buttonExit, 2, wx.ALL, 2)
        
        mainSizer.Add(hsiz2, 0, wx.ALL, 5)
        
        self.SetSizerAndFit(mainSizer)

    def EvtRadioBox(self, event):
        self.logger.AppendText('plot: %d\n' % event.GetInt())
    def ChooseSim(self, event):
        self.logger.AppendText('simulation: %s\n' % event.GetString())
    def OnPlot(self,event):
        self.logger.AppendText("Plot!\n")# %event.GetId())
    def OnExit(self,event):
        self.Close() #TODO: check
    def EvtText(self, event):
        self.logger.AppendText('Text: %s\n' % event.GetString())
    def EvtChar(self, event):
        self.logger.AppendText('Char: %d\n' % event.GetKeyCode())
        event.Skip()
    def EvtCheckBox(self, event):
        self.logger.AppendText('movie: %d\n' % event.Checked())



app = wx.App(False)
frame = wx.Frame(None, title="Twiddle")
nb = wx.Notebook(frame)

nb.AddPage(panel(nb), "Analysis")
#nb.AddPage(panel(nb), "2nd sim")
#nb.AddPage(panel(nb), "3rd sim")
frame.Show()
app.MainLoop()

# app = wx.App(False)
# frame = wx.Frame(None)
# panel = ExamplePanel(frame)
# frame.Show()
# app.MainLoop()
