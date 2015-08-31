#!/usr/bin/env python
import matplotlib
matplotlib.use('TkAgg')



from numpy import arange, sin, pi,log10,max,min,cos,isnan, meshgrid,sqrt,abs
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import pyPLUTO as pp
import string
import time

from Tkinter import *
import sys
import os



class App:
    def __init__(self,master):

        # create toplevel window
        frame = Frame(master)
        
        frame.grid(ipadx=10,ipady=10)

        try:
            sys.argv[1]
        except:
            self.datatype = None
        else:
            self.datatype = sys.argv[1].split('--')[1]
            
        if self.datatype == 'hdf5':
            print "GUI currently doesnot support pyPLUTO AMR Reader!!"
            sys.exit()
        
        self.I = pp.Image()
        self.Tool = pp.Tools()
        


        self.lb1=Label(frame, text="Nstep").grid(row=0,column=0)
        
        
        
        self.enstep = Entry(frame,width=8)
        self.enstep.grid(row=0,column=1)
        self.enstep.insert(0, "0")

        self.LoadedNstep = StringVar()
        self.PresentTime = StringVar()

        self.myData = self.loaddata()
        self.varkeys = self.myData.vars
        self.wdir = self.myData.wdir
    

        if self.myData.n3 != 1:
            self.Geom = '3D'
        elif self.myData.n3 == 1 and self.myData.n2 != 1:
            self.Geom = '2D'
        else:
            self.Geom = '1D'
            
        

        self.ldatabutton=Button(frame,text="Load data",command=self.loaddata)
        self.ldatabutton.grid(row=0,column=2)

############### MARK THE CUTS #################################

        self.ex1 = Entry(frame,width=5)
        self.ex1.grid(row=2,column=0)
        self.ex1.insert(0, "x1")

        self.ex2 = Entry(frame,width=5)
        self.ex2.grid(row=2,column=1)
        self.ex2.insert(0, "x2")
        
        self.ex3 = Entry(frame,width=5)
        self.ex3.grid(row=2,column=2)
        self.ex3.insert(0, "x3")

        if self.Geom == '2D':
            self.ex3.config(state='disabled')

        if self.Geom == '1D':
            self.ex3.config(state='disabled')
            self.ex2.config(state='disabled')
            self.ex1.config(state='disabled')
        

        # place a graph somewhere here
        self.f = Figure(figsize=(7,7), dpi=100)
        self.a = self.f.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.f, master=root)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(row=0,column=3,columnspan=10,rowspan=10,sticky=E)

        #self.toolbar = NavigationToolbar2TkAgg(self.canvas,tl)
        #self.toolbar.update()
        #self.canvas._tkcanvas.grid(row=60,column=15,sticky=E)

        self.v = StringVar()
        self.v.set("None")

################ VARIABLES TO PLOT #################################
        for i in ['bx1s', 'bx2s', 'bx3s']:
            try:
                self.varkeys.remove(i)
            except ValueError:
                pass
        
        for j in range(len(self.varkeys)):
            self.ldata = Radiobutton(frame,text=self.varkeys[j],variable=self.v,value=self.varkeys[j],command=self.getmyvar)
            self.ldata.grid(row=3+j,column=0,sticky=W)
        

################ SLICES CHOICE #################################
            
        self.slvar = StringVar()
        self.slvar.set("Choose Slice")
        if self.Geom == '3D' :
            SliceList = ("Along x1","Along x2","Along x3","Along x1-x2","Along x2-x3","Along x3-x1")
        elif self.Geom == '2D' :
            SliceList = ("Along x1", "Along x2", "Along x1-x2")
        else:
            SliceList = ()

        for j in range(len(SliceList)):
            self.sldata = Radiobutton(frame,text=SliceList[j],variable=self.slvar,value=SliceList[j],command=self.setslice)
            self.sldata.grid(row=3+j,column=1,sticky=W)

############### PLOT PROPERTIES #################################
            
        self.logvar = IntVar()
        self.chkb = Checkbutton(frame,text="Log  ",variable=self.logvar,onvalue=1,offvalue=0,command=self.logchkcall)
        self.chkb.grid(row=3,column=2,sticky=W)#(row=15,column=0,sticky=W)

        self.polarvar = IntVar()
        self.polchkb = Checkbutton(frame,text="Polar",variable=self.polarvar,onvalue=1,offvalue=0,command=self.polchkcall)
        self.polchkb.grid(row=4,column=2,sticky=W)#(row=15,column=1)
        if self.Geom == '1D':
            self.polchkb.config(state='disabled')
            self.polarvar.set(0)

        self.preaspect = IntVar()
        self.aspectb = Checkbutton(frame,text="Aspect",variable=self.preaspect,onvalue=1,offvalue=0,command=self.aspchkcall)
        self.aspectb.grid(row=5,column=2,sticky=W)#(row=15,column=2)
        if self.Geom == '1D':
            self.aspectb.config(state='disabled')

        

        

        
################ X and Y LABELS #################################
        
        self.lb2=Label(frame,text="Labels").grid(row=22,column=0)
    
        self.xlb = Entry(frame,width=15)
        self.xlb.grid(row=22,column=1)
        self.xlb.insert(0, "xlabel")

        self.ylb = Entry(frame,width=15)
        self.ylb.grid(row=22,column=2)
        self.ylb.insert(0, "ylabel")

############### X and Y RANGE#######################

        self.lb2a=Label(frame,text="XRange").grid(row=24,column=0)
        self.lb2b=Label(frame,text="YRange").grid(row=26,column=0)
        self.lb2c=Label(frame,text="VarRange").grid(row=28,column=0)
        
        self.xrmin = Entry(frame,width=15)
        self.xrmin.grid(row=24,column=1)
        self.xrmin.insert(0,'')
        self.xrmax = Entry(frame,width=15)
        self.xrmax.grid(row=24,column=2)
        self.xrmax.insert(0,'')
        
        self.yrmin = Entry(frame,width=15)
        self.yrmin.grid(row=26,column=1)
        self.yrmin.insert(0,'')
        self.yrmax = Entry(frame,width=15)
        self.yrmax.grid(row=26,column=2)
        self.yrmax.insert(0,'')

        self.varmin = Entry(frame,width=15)
        self.varmin.grid(row=28,column=1)
        self.varmin.insert(0,'')
        self.varmax = Entry(frame,width=15)
        self.varmax.grid(row=28,column=2)
        self.varmax.insert(0,'')
        if self.Geom == '1D':
            self.yrmin.config(state='disabled')
            self.yrmax.config(state='disabled')


################ CONTOURS #################################

        self.lb3=Label(frame,text="Contours").grid(row=16,column=0)
        
        self.contvar = IntVar()
        self.chkb = Checkbutton(frame,text="Contour",variable=self.contvar,onvalue=1,offvalue=0,command=self.contchkcall)
        self.chkb.grid(row=6,column=2,sticky=W)#(row=16,column=0,sticky=W)

        self.plcont = StringVar()
        self.contkeys = ["None"]
        if "bx3" in self.varkeys:
            for item in self.varkeys:
                self.contkeys.append(item)
            self.contkeys.append("x1*bx3")
            if "Ax3" in self.varkeys:
                self.contkeys.append("x1*Ax3")
        else:
            for item in self.varkeys:
                self.contkeys.append(item)        
        self.plcont.set("None")
        self.contmenu = OptionMenu(frame, self.plcont,*self.contkeys)
        self.contmenu.grid(row=16,column=1)

        self.xlevb = Entry(frame,width=15)
        self.xlevb.grid(row=16,column=2,sticky=W)
        self.xlevb.insert(0, "Levels")
        self.xlevb.config(state='disabled')
        self.contmenu.config(state='disabled')

        if self.Geom == '1D':
            self.chkb.config(state = 'disabled')
        
################ ARROWS #################################

        self.lb4=Label(frame,text="Arrows").grid(row=19,column=0)

        self.arrowvar = IntVar()
        self.arrowchkb = Checkbutton(frame,text="Arrows",variable=self.arrowvar,onvalue=1,offvalue=0,command=self.arrchkcall)
        self.arrowchkb.grid(row=7,column=2,sticky=W)#(row=16,column=0,sticky=W)

        self.arrspb = Entry(frame,width=15)
        self.arrspb.grid(row=19,column=2,sticky=W)
        self.arrspb.insert(0, "20")

    
        self.plarr = StringVar()
        self.arrkeys = ["None"]
        self.arrkeys.append("Vp")
        self.arrkeys.append("Vp_norm")
        if "bx1" in self.varkeys:
            self.arrkeys.append("Bp")
            self.arrkeys.append("Bp_norm")
        self.plarr.set("None")
        self.arrmenu = OptionMenu(frame,self.plarr,*self.arrkeys)
        self.arrmenu.grid(row=19,column=1)
        self.arrmenu.config(state='disabled')
        self.arrspb.config(state='disabled')

        if self.Geom == '1D':
            self.arrowchkb.config(state = 'disabled')
        

                
################ VARIOUS PLOTTING BUTTONS #################################

        self.pltbutton=Button(frame,text="Plot",command=self.plotfinal)
        self.pltbutton.grid(row=36,column=0)
        if self.Geom == '1D':
            self.pltbutton.config(state='active')
        else:
            self.pltbutton.config(state='disabled')

        self.surfbutton=Button(frame,text="Surface",command=self.plotsurface)
        self.surfbutton.grid(row=36,column=1)
        self.surfbutton.config(state='disabled')
        #if self.Geom == '1D':
        #    self.surfbutton.config(state='disabled')
        
        
        self.clrbutton=Button(frame,text="Clear",command=self.plotclear)
        self.clrbutton.grid(row=36,column=2)

################ INFORMATION #################################

        self.lbinf0 =  Label(frame,text="Information",font=("Times",12,"bold"))
        self.lbinf0.grid(row=47,column=0,sticky=W,columnspan=3)

        self.lbinf1a = Label(frame,text="Dir :",font=("Times",10,"bold")).grid(row=49,column=0,sticky=W,columnspan=3)
        self.lbinf1 =  Label(frame,text=self.wdir).grid(row=50,column=0,sticky=W,columnspan=3)
        self.lbinf2a = Label(frame,text="Domain :",font=("Times",10,"bold")).grid(row=51,column=0,sticky=W,columnspan=3)
        self.lbinf2 = Label(frame,text="n1 x n2 x n3 =  %d x %d x %d " % (self.myData.n1,self.myData.n2,self.myData.n3)).grid(row=52,column=0,sticky=W,columnspan=3)
        self.lbinf3a = Label(frame,text="Time Status",font=("Times",10,"bold")).grid(row=53,column=0,sticky=W,columnspan=3)
        self.lbinf4 = Label(frame,text="Nlast = %d"% pp.nlast_info(w_dir=self.wdir,datatype=self.datatype)['nlast']).grid(row=54,column=0,sticky=W,columnspan=3)
        self.lbinf5 = Label(frame,textvariable = self.LoadedNstep).grid(row=55,column=0,sticky=W,columnspan=3)
        self.lbinf6 = Label(frame,textvariable = self.PresentTime).grid(row=56,column=0,sticky=W,columnspan=3)
    

################ VARIOUS FUNCTIONS #################################
        
    
                                
    def loaddata(self):
        try:
            int(self.enstep.get().strip().split()[0])
        except (ValueError, IndexError):
            print "Specify the proper value of Nstep"
        else:
            mynstep=int(self.enstep.get())
            self.D = pp.pload(mynstep,datatype=self.datatype)
            self.LoadedNstep.set("Loaded Nstep = "+self.enstep.get())
            self.PresentTime.set("Present Time = "+str(self.D.SimTime) + " [cu]")
            return self.D


    def getmyvar(self):
        try:
            self.v.get() != "None"
        except KeyError:
            print "Specify the variable to plot"
        else:
            self.myvar=self.v.get()
        
    def logchkcall(self):
        self.logchk = self.logvar.get()

    def contchkcall(self):
        self.contchk = self.contvar.get()
        if self.contchk == 1:
            self.contmenu.config(state='normal')
            self.xlevb.config(state='normal')
        else:
            self.contmenu.config(state='disabled')
            self.xlevb.config(state='disabled')

    def arrchkcall(self):
        self.arrchk = self.arrowvar.get()
        if self.arrchk == 1:
            self.arrmenu.config(state='normal')
            self.arrspb.config(state='normal')
        else:
            self.arrmenu.config(state='disabled')
            self.arrspb.config(state='disabled')

    def aspchkcall(self):
        self.aspchk=self.preaspect.get()

    def polchkcall(self):
        self.polchk = self.polarvar.get()
        
    def setslice(self):
        self.slicename=self.slvar.get()
        if self.slicename == "Along x1" or self.slicename == "Along x2" or self.slicename == "Along x3":
            self.surfbutton.config(state='disabled')
            self.arrowchkb.config(state = 'disabled')
            self.arrowvar.set(0)
            self.chkb.config(state = 'disabled')
            self.contvar.set(0)
            self.pltbutton.config(state='active')
            self.polchkb.config(state='disabled')
            self.polarvar.set(0)
        else:
            self.pltbutton.config(state='disabled')
            self.arrowchkb.config(state = 'normal')
            self.chkb.config(state = 'normal')
            self.surfbutton.config(state='active')
            self.polchkb.config(state='normal')
        
        if self.slicename == "Along x2-x3":
            self.polchkb.config(state='disabled')
            self.polarvar.set(0)
            
        
    
    def plotclear(self):

        self.f.clf()
        self.a = self.f.add_subplot(111)        
        self.canvas.show()

    def plotfinal(self):
        if self.getplotvar() == True:
            self.a.axis([self.getxaxisrange()[0],self.getxaxisrange()[1],self.getvarrange()[0],self.getvarrange()[1]])
            self.a.plot(self.x,self.var)
            self.a.set_aspect('auto')
        self.a.set_xlabel(self.xlb.get())
        self.a.set_ylabel(self.ylb.get())
        self.canvas.show()

    def plotsurface(self):
        tdum = time.time()
        self.plotclear()
        
        if self.preaspect.get() == 1:
            self.a.set_aspect('equal')
        else:
            self.a.set_aspect('auto')

        
        if self.polarvar.get() == 1:
            if self.drawpolar() == True:
                self.a.axis([self.getxaxisrange()[0],self.getxaxisrange()[1],self.getyaxisrange()[0],self.getyaxisrange()[1]])
                self.image = self.a.imshow(self.SphData[self.myvar], origin='lower',extent=self.extent, interpolation='nearest',cmap="jet", vmin=self.getvarrange()[0],vmax=self.getvarrange()[1])
                self.f.colorbar(self.image)
        else:
            if self.getsurfvar() == True:
                self.a.axis([self.getxaxisrange()[0],self.getxaxisrange()[1],self.getyaxisrange()[0],self.getyaxisrange()[1]])
                self.image=self.a.pcolormesh(self.x,self.y,self.var,cmap='jet',vmin=self.getvarrange()[0],vmax=self.getvarrange()[1])
                self.f.colorbar(self.image)
        

        if self.contvar.get() == 1:
            try:
                self.plcont.get() != "None"
            except KeyError:
                print "Specify the variable for Contour"
            else:
                self.drawcontour()
                self.contlevlist=[]
                self.contlevstr = string.split(self.xlevb.get(),',')
                try:
                    if self.contlevstr[0] == 'log':
                        self.flevel = self.contlevstr[1]
                        self.varcont = log10(self.varcont)
                    else:
                        self.flevel = self.contlevstr[0]
                    
                    float(self.flevel)
                    self.contlevlist = [float(self.flevel)]
                    
                except:
                    self.contlevlist = 5
                else:
                    for j in range(1,len(self.contlevstr)):
                        self.contlevlist.append(float(self.contlevstr[j]))


                self.cs1 = self.a.contour(self.xcont,self.ycont,self.varcont,self.contlevlist,colors="w")
                self.a.clabel(self.cs1,inline=True)
        

        if self.arrowvar.get() == 1:
            try:
                self.plarr.get() != "None"
            except KeyError:
                print "Specify the variable for plotting the arrow"
            else:
                self.drawarrow()
                self.a.quiver(self.xcong, self.ycong, self.xveccong, self.yveccong,color='w')
        
        self.a.set_xlabel(self.xlb.get())
        self.a.set_ylabel(self.ylb.get())
        self.canvas.show()
        


    def getvarrange(self):
        try:
            float(self.varmin.get())
        except:
            if self.polarvar.get() != 1:
                self.varminval = min(self.var)
            else:
                self.varminval = min(self.SphData[self.myvar][self.isnotnan].flat)#self.minPl
        else:
            self.varminval = float(self.varmin.get())

        try:
            float(self.varmax.get())
        except:
            if self.polarvar.get() != 1:
                self.varmaxval = max(self.var)
            else:
                self.varmaxval = max(self.SphData[self.myvar][self.isnotnan].flat)#self.maxPl
        else:
            self.varmaxval = float(self.varmax.get())

        return [self.varminval,self.varmaxval]
        

    def getxaxisrange(self):
        try:
            float(self.xrmin.get())
        except:
            if self.polarvar.get() != 1:
                self.xminval = min(self.x)
            else:
                self.xminval = min(self.R.flat)
        else:
            self.xminval = float(self.xrmin.get())

        try:
            float(self.xrmax.get())
        except:
            if self.polarvar.get() != 1:
                self.xmaxval = max(self.x)
            else:
                self.xmaxval = max(self.R.flat)
        else:
            self.xmaxval = float(self.xrmax.get())

        return [self.xminval,self.xmaxval]

    def getyaxisrange(self):
        try:
            float(self.yrmin.get())
        except:
            if self.polarvar.get() != 1:
                self.yminval = min(self.y)
            else:
                self.yminval = min(self.Z.flat)
        else:
            self.yminval = float(self.yrmin.get())

        try:
            float(self.yrmax.get())
        except:
            if self.polarvar.get() != 1:
                self.ymaxval = max(self.y)
            else:
                self.ymaxval = max(self.Z.flat)
        else:
            self.ymaxval = float(self.yrmax.get())

        return [self.yminval,self.ymaxval]


    def getplotvar(self):
        self.sucess = False
        if self.logvar.get() == 1:
            self.var = log10(self.D.__getattribute__(self.myvar))
        else:
            self.var = self.D.__getattribute__(self.myvar)

        if self.Geom == '1D':
            self.x = self.D.x1
            self.sucess = True
        else:
            if self.slicename == "Along x1":
                self.x = self.D.x1
                if self.D.n3 == 1:
                    try:
                        int(self.ex2.get().strip().split()[0])
                    except (ValueError, IndexError):
                        print "Specify the value of x2 cut"
                    else:
                        self.var = self.var[:,int(self.ex2.get())]
                        self.sucess = True
                else:
                    try:
                        int(self.ex2.get().strip().split()[0])
                        int(self.ex3.get().strip().split()[0])
                    except (ValueError, IndexError):
                        print "Specify the value of x2 or x3 cut"
                    else:
                        self.var = self.var[:,int(self.ex2.get()),int(self.ex3.get())]
                        self.sucess = True
            
            elif self.slicename == "Along x2":
                self.x = self.D.x2
                if self.D.n3 == 1:
                    try:
                        int(self.ex1.get().strip().split()[0])
                    except (ValueError, IndexError):
                        print "Specify the value of x1 cut"
                    else:
                        self.var = self.var[int(self.ex1.get()),:]
                        self.sucess = True
                else:
                    try:
                        int(self.ex1.get().strip().split()[0])
                        int(self.ex3.get().strip().split()[0])
                    except (ValueError, IndexError):
                        print "Specify the value of x1 or x3 cut"
                    else:
                        self.var = self.var[int(self.ex1.get()),:,int(self.ex3.get())]
                        self.sucess = True

            else:
                self.x = self.D.x3
                try:
                    int(self.ex1.get().strip().split()[0])
                    int(self.ex2.get().strip().split()[0])
                except (ValueError, IndexError):
                    print "Specify the value of x1 or x2 cut"
                else:
                    self.var = self.var[int(self.ex1.get()),int(self.ex2.get()),:]
                    self.sucess = True

        return self.sucess

    def getsurfvar(self):
        self.sucess = False
        if self.logvar.get() == 1:
            self.var = log10(self.D.__getattribute__(self.myvar))
        else:
            self.var = self.D.__getattribute__(self.myvar)
        
        if self.slicename == "Along x1-x2":
            self.x = self.D.x1
            self.y = self.D.x2
            xmineed = (abs(self.x-self.getxaxisrange()[0])).argmin()
            xmaneed = (abs(self.x-self.getxaxisrange()[1])).argmin()
            ymineed = (abs(self.y-self.getyaxisrange()[0])).argmin()
            ymaneed = (abs(self.y-self.getyaxisrange()[1])).argmin()
            self.x = self.x[xmineed:xmaneed]
            self.y = self.y[ymineed:ymaneed]
            if self.D.n3 == 1:
                self.var = self.var[xmineed:xmaneed,ymineed:ymaneed].T
                self.sucess = True
            else:
                try:
                    int(self.ex3.get().strip().split()[0])
                except (ValueError, IndexError):
                    print "Specify the value of x3 cut"
                else:
                    self.var = self.var[xmineed:xmaneed,ymineed:ymaneed,int(self.ex3.get())].T
                    self.sucess = True
                
        elif self.slicename == "Along x2-x3":
            self.x = self.D.x2
            self.y = self.D.x3
            xmineed = (abs(self.x-self.getxaxisrange()[0])).argmin()
            xmaneed = (abs(self.x-self.getxaxisrange()[1])).argmin()
            ymineed = (abs(self.y-self.getyaxisrange()[0])).argmin()
            ymaneed = (abs(self.y-self.getyaxisrange()[1])).argmin()
            self.x = self.x[xmineed:xmaneed]
            self.y = self.y[ymineed:ymaneed]
            try:
                int(self.ex1.get().strip().split()[0])
            except (ValueError, IndexError):
                print "Specify the value of x1 cut"
            else:
                self.var = self.var[int(self.ex1.get()),xmineed:xmaneed,ymineed:ymaneed].T
                self.sucess = True
        else:
            self.x = self.D.x1
            self.y = self.D.x3
            xmineed = (abs(self.x-self.getxaxisrange()[0])).argmin()
            xmaneed = (abs(self.x-self.getxaxisrange()[1])).argmin()
            ymineed = (abs(self.y-self.getyaxisrange()[0])).argmin()
            ymaneed = (abs(self.y-self.getyaxisrange()[1])).argmin()
            self.x = self.x[xmineed:xmaneed]
            self.y = self.y[ymineed:ymaneed]
            try:
                int(self.ex2.get().strip().split()[0])
            except (ValueError, IndexError):
                print "Specify the value of x2 cut"
            else:
                self.var = self.var[xmineed:xmaneed,int(self.ex2.get()),ymineed:ymaneed].T
                self.sucess = True
        
        return self.sucess

    def drawpolar(self):
        self.sucess = False
        if self.slicename == "Along x1-x2":
            if self.D.n3 == 1:
                self.R,self.Z,self.SphData = self.I.getSphData(self.D,w_dir=self.wdir,datatype=self.datatype, rphi=False)
                self.sucess = True
            else:
                try:
                    int(self.ex3.get().strip().split()[0])
                except (ValueError, IndexError):
                    print "Specify the value of x3 cut"
                else:
                    self.R,self.Z,self.SphData = self.I.getSphData(self.D,w_dir=self.wdir,datatype=self.datatype, rphi=False,x3cut=int(self.ex3.get()))
                    self.sucess = True

        if self.slicename == "Along x3-x1":
            try:
                int(self.ex2.get().strip().split()[0])
            except (ValueError, IndexError):
                print "Specify the value of x2 cut"
            else:
                self.R,self.Z,self.SphData = self.I.getSphData(self.D,w_dir=self.wdir,datatype=self.datatype, rphi=True, x2cut=int(self.ex2.get()))
                self.sucess = True

        if self.sucess == True:
            self.extent=(min(self.R.flat),max(self.R.flat),min(self.Z.flat),max(self.Z.flat))
            self.dRR=max(self.R.flat)-min(self.R.flat)
            self.dZZ=max(self.Z.flat)-min(self.Z.flat)

            self.isnotnan=-isnan(self.SphData[self.myvar])
            self.maxPl=max(self.SphData[self.myvar][self.isnotnan].flat)
            self.minPl=min(self.SphData[self.myvar][self.isnotnan].flat)
            self.normrange=False
            if self.minPl<0:
                self.normrange=True
            if self.maxPl>-self.minPl:
                self.minPl=-self.maxPl
            else:
                self.maxPl=-self.minPl	  
            if (self.normrange and self.myvar !='rho' and self.myvar !='prs'):
                self.SphData[self.myvar][-1][-1]=self.maxPl
                self.SphData[self.myvar][-1][-2]=self.minPl
            if self.logvar.get() == 1:
                self.SphData[self.myvar] = log10(self.SphData[self.myvar])

        return self.sucess


    def drawcontour(self):
        if self.polarvar.get() != 1:
            if self.slicename == "Along x1-x2":
                self.xcont = self.D.x1
                self.ycont = self.D.x2
                self.Xmesh, self.Ymesh = meshgrid(self.D.x1.T,self.D.x2.T)
                if self.D.n3 == 1:
                    if self.plcont.get() == 'x1*Ax3':
			self.varcont = self.Xmesh*(self.D.Ax3.T)
                    elif self.plcont.get() == 'x1*bx3':
			self.varcont = self.Xmesh*(self.D.bx3.T)
                    else:
                        self.varcont = self.D.__getattribute__(self.plcont.get())[:,:].T              
                else:
                    if self.plcont.get() == 'x1*Ax3':
			self.varcont = self.Xmesh*(self.D.Ax3[:,:,int(self.ex3.get())].T)
                    elif self.plcont.get() == 'x1*bx3':
			self.varcont = self.Xmesh*(self.D.bx3[:,:,int(self.ex3.get())].T)
                    else:
                        self.varcont = self.D.__getattribute__(self.plcont.get())[:,:,int(self.ex3.get())].T
                        
            elif self.slicename == "Along x2-x3":
                self.xcont = self.D.x2
                self.ycont = self.D.x3
                self.varcont = self.D.__getattribute__(self.plcont.get())[int(self.ex1.get()),:,:].T
            else:
                self.xcont = self.D.x1
                self.ycont = self.D.x3
                self.varcont = self.D.__getattribute__(self.plcont.get())[:,int(self.ex2.get()),:].T
        else:
            self.xcont = self.R
            self.ycont = self.Z
            if self.plcont.get() == 'x1*Ax3':
                self.varcont = self.R*(self.SphData['Ax3'])
            elif self.plcont.get() == 'x1*bx3':
                self.varcont = self.R*(self.SphData['bx3'])
            else:
                if self.logvar.get() == 1 and self.plcont.get() == self.myvar:
                    self.varcont = 10**(self.SphData[self.plcont.get()])
                else:
                    self.varcont = self.SphData[self.plcont.get()]

    def drawarrow(self):
        if self.polarvar.get() != 1:
            if self.slicename == "Along x1-x2":
                self.Xmesh, self.Ymesh = meshgrid(self.D.x1.T,self.D.x2.T)
                self.xcong = self.Tool.congrid(self.Xmesh,2*(int(self.arrspb.get()),),method='linear')
                self.ycong = self.Tool.congrid(self.Ymesh,2*(int(self.arrspb.get()),),method='linear')
                if self.plarr.get() == 'Vp' or self.plarr.get() =='Vp_norm':
                    if self.D.n3 == 1:
                        self.vel1 = self.D.vx1[:,:].T
                        self.vel2 = self.D.vx2[:,:].T
                    else:
                        self.vel1 = self.D.vx1[:,:,int(self.ex3.get())].T
                        self.vel2 = self.D.vx2[:,:,int(self.ex3.get())].T
                        
                    self.xveccong = self.Tool.congrid(self.vel1,2*(int(self.arrspb.get()),),method='linear')
                    self.yveccong = self.Tool.congrid(self.vel2,2*(int(self.arrspb.get()),),method='linear')
                    self.normVp = sqrt(self.xveccong**2 + self.yveccong**2)
                    if self.plarr.get() == 'Vp_norm':
                        self.xveccong = self.xveccong/self.normVp
                        self.yveccong = self.yveccong/self.normVp
                if self.plarr.get() == 'Bp' or self.plarr.get() =='Bp_norm':
                    if self.D.n3 == 1:
                        self.mag1 = self.D.bx1[:,:].T
                        self.mag2 = self.D.bx2[:,:].T
                    else:
                        self.mag1 = self.D.bx1[:,:,int(self.ex3.get())].T
                        self.mag2 = self.D.bx2[:,:,int(self.ex3.get())].T
                        
                    self.xveccong = self.Tool.congrid(self.mag1,2*(int(self.arrspb.get()),),method='linear')
                    self.yveccong = self.Tool.congrid(self.mag2,2*(int(self.arrspb.get()),),method='linear')
                    self.normVp = sqrt(self.xveccong**2 + self.yveccong**2)
                    if self.plarr.get() == 'Bp_norm':
                        self.xveccong = self.xveccong/self.normVp
                        self.yveccong = self.yveccong/self.normVp
                
            elif self.slicename == "Along x2-x3":
                self.Xmesh, self.Ymesh = meshgrid(self.D.x2.T,self.D.x3.T)
                self.xcong = self.Tool.congrid(self.Xmesh,2*(int(self.arrspb.get()),),method='linear')
                self.ycong = self.Tool.congrid(self.Ymesh,2*(int(self.arrspb.get()),),method='linear')
                if self.plarr.get() == 'Vp' or self.plarr.get() =='Vp_norm':
                    self.vel1 = self.D.vx2[int(self.ex1.get()),:,:].T
                    self.vel2 = self.D.vx3[int(self.ex1.get()),:,:].T
                    self.xveccong = self.Tool.congrid(self.vel1,2*(int(self.arrspb.get()),),method='linear')
                    self.yveccong = self.Tool.congrid(self.vel2,2*(int(self.arrspb.get()),),method='linear')
                    self.normVp = sqrt(self.xveccong**2 + self.yveccong**2)
                    if self.plarr.get() == 'Vp_norm':
                        self.xveccong = self.xveccong/self.normVp
                        self.yveccong = self.yveccong/self.normVp
                if self.plarr.get() == 'Bp' or self.plarr.get() =='Bp_norm':
                    self.mag1 = self.D.bx2[int(self.ex1.get()),:,:].T
                    self.mag2 = self.D.bx3[int(self.ex1.get()),:,:].T
                    self.xveccong = self.Tool.congrid(self.mag1,2*(int(self.arrspb.get()),),method='linear')
                    self.yveccong = self.Tool.congrid(self.mag2,2*(int(self.arrspb.get()),),method='linear')
                    self.normVp = sqrt(self.xveccong**2 + self.yveccong**2)
                    if self.plarr.get() == 'Bp_norm':
                        self.xveccong = self.xveccong/self.normVp
                        self.yveccong = self.yveccong/self.normVp
            else:
                self.Xmesh, self.Ymesh = meshgrid(self.D.x1.T,self.D.x3.T)
                self.xcong = self.Tool.congrid(self.Xmesh,2*(int(self.arrspb.get()),),method='linear')
                self.ycong = self.Tool.congrid(self.Ymesh,2*(int(self.arrspb.get()),),method='linear')
                if self.plarr.get() == 'Vp' or self.plarr.get() =='Vp_norm':
                    self.vel1 = self.D.vx1[:,int(self.ex2.get()),:].T
                    self.vel2 = self.D.vx3[:,int(self.ex2.get()),:].T
                    self.xveccong = self.Tool.congrid(self.vel1,2*(int(self.arrspb.get()),),method='linear')
                    self.yveccong = self.Tool.congrid(self.vel2,2*(int(self.arrspb.get()),),method='linear')
                    self.normVp = sqrt(self.xveccong**2 + self.yveccong**2)
                    if self.plarr.get() == 'Vp_norm':
                        self.xveccong = self.xveccong/self.normVp
                        self.yveccong = self.yveccong/self.normVp
                if self.plarr.get() == 'Bp' or self.plarr.get() =='Bp_norm':
                    self.mag1 = self.D.bx1[:,int(self.ex2.get()),:].T
                    self.mag2 = self.D.bx3[:,int(self.ex2.get()),:].T
                    self.xveccong = self.Tool.congrid(self.mag1,2*(int(self.arrspb.get()),),method='linear')
                    self.yveccong = self.Tool.congrid(self.mag2,2*(int(self.arrspb.get()),),method='linear')
                    self.normVp = sqrt(self.xveccong**2 + self.yveccong**2)
                    if self.plarr.get() == 'Bp_norm':
                        self.xveccong = self.xveccong/self.normVp
                        self.yveccong = self.yveccong/self.normVp
        else:
            self.xcong = self.Tool.congrid(self.R,2*(int(self.arrspb.get()),),method='linear')
            self.ycong = self.Tool.congrid(self.Z,2*(int(self.arrspb.get()),),method='linear')
            if self.plarr.get() == 'Vp' or self.plarr.get() =='Vp_norm':
                if self.slicename == "Along x1-x2":
                    self.vel1 = self.SphData['v1c']
                    self.vel2 = self.SphData['v2c']
                else:
                    self.vel1 = self.SphData['v1c']
                    self.vel2 = self.SphData['v3c']
                
                self.xveccong = self.Tool.congrid(self.vel1,2*(int(self.arrspb.get()),),method='linear')
                self.yveccong = self.Tool.congrid(self.vel2,2*(int(self.arrspb.get()),),method='linear')
                self.normVp = sqrt(self.xveccong**2 + self.yveccong**2)
                if self.plarr.get() == 'Vp_norm':
                    self.xveccong = self.xveccong/self.normVp
                    self.yveccong = self.yveccong/self.normVp
            if self.plarr.get() == 'Bp' or self.plarr.get() =='Bp_norm':
                
                if self.slicename == "Along x1-x2":
                    self.mag1 = self.SphData['b1c']
                    self.mag2 = self.SphData['b2c']
                else:
                    self.mag1 = self.SphData['b1c']
                    self.mag2 = self.SphData['b3c']
                
                self.xveccong = self.Tool.congrid(self.mag1,2*(int(self.arrspb.get()),),method='linear')
                self.yveccong = self.Tool.congrid(self.mag2,2*(int(self.arrspb.get()),),method='linear')
                self.normVp = sqrt(self.xveccong**2 + self.yveccong**2)
                if self.plarr.get() == 'Bp_norm':
                    self.xveccong = self.xveccong/self.normVp
                    self.yveccong = self.yveccong/self.normVp
        
    def epssave(self):
        self.f.savefig(self.myvar+'_'+self.enstep.get()+'.eps')
    def pngsave(self):
        self.f.savefig(self.myvar+'_'+self.enstep.get()+'.png')
    def pdfsave(self):
        self.f.savefig(self.myvar+'_'+self.enstep.get()+'.pdf')
    def jpgsave(self):
        self.f.savefig(self.myvar+'_'+self.enstep.get()+'.jpg')
    
    



    
    
            
root=Tk()
app=App(root)
root.title("pyPLUTO")

menubar = Menu(root)
savemenu = Menu(menubar,tearoff=0)
savemenu.add_command(label='EPS',command=app.epssave)
savemenu.add_command(label='PDF',command=app.pdfsave)
savemenu.add_command(label='PNG',command=app.pngsave)
savemenu.add_command(label='JPG',command=app.jpgsave)
menubar.add_cascade(label="Save As", menu=savemenu)



#menubar.add_command(label='Plot',command = app.plotfinal)
#menubar.add_command(label='Surface',command=app.plotsurface)
#menubar.add_command(label='Clear',command=app.plotclear)
menubar.add_command(label='Quit',command=root.quit)

root.config(menu=menubar)

root.mainloop()   

