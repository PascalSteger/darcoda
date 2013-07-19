#!/usr/bin/python2.7
# class to generate all filenames, independent of which investigation is being performed

import gl_params as gp
import gl_helper as gh
import pdb

if gp.geom == 'sphere':
    import physics_sphere as phys
elif gp.geom == 'disc':
    import physics_disc as phys



    
def get_case(cas):
    # Set number of tracers to look at
    # want to set ntracer = 3e3              # case 1
    #             ntracer = 3e4              # case 2
    
    ntracer = 0
    if cas == 1:
        ntracer = 3000
    elif cas == 2:
        ntracer = 30000
    return ntracer




    
class Files:
    'Common base class for all filename sets'
    # massfile, nufile[], sigfile[], outname = Files(investigate)


    
    def __init__ (self):
        self.machine = ''
        self.dir = ''
        self.progdir = ''
        self.set_dir(gp.machine) # 'darkside' or 'local'
        self.massfile = ''; self.analytic = '';               self.surfdenfiles = []
        self.nufiles  = [];        self.sigfiles = [];        self.posvelfiles = []
        self.ntracer, self.nstr1 = self.set_ntracer(gp.cas)
        
        if gp.investigate == 'hernquist':
            self.set_hernquist()
        elif gp.investigate == 'walker':
            self.set_walker()
        elif gp.investigate == 'triaxial':
            self.set_triaxial()
        elif gp.investigate == 'fornax':
            self.set_fornax()
        elif gp.investigate == 'sim':
            self.set_disc_sim()
        elif gp.investigate == 'simple':
            self.set_disc_simple()

        self.outdir, self.outname = self.get_outname()
        return




    
    def set_dir(self, machine):
        if machine == 'darkside':
            self.machine = '/home/ast/read/dark/dwarf_data/'
        elif machine == 'local':
            self.machine = '/home/psteger/sci/dwarf_data/'
        self.progdir = self.machine + 'programs'
        return




        
    def set_ntracer(self,cas):
        ntracer = get_case(cas)
        self.ntracer = ntracer
        self.nstr = str(ntracer)
        return self.ntracer, self.nstr



    



        
    def get_sim_name(self):
        if gp.pops == 1:
            sim = 'unit_hern_1_'
        elif gp.pops == 2:
            sim = 'dual_unit_hern_1_'
        return sim




        
    def set_hernquist(self):
        '''set all filenames for Hernquist case'''
        self.dir = self.machine + 'data_hernquist/'
        sim = self.get_sim_name()
        self.massfile = self.dir+'enclosedmass/'+sim+'enclosedmass.txt'
        if gp.pops == 1:
            if gp.checksigma and gp.cas == 0:
                self.nufiles.append(self.dir+'densityfalloff/'+sim+'falloffnotnorm.txt') # all comp.
                self.nufiles.append(self.dir+'densityfalloff/'+sim+'falloffnotnorm.txt') # first comp.
                self.sigfiles.append(self.dir+'velocitydispersionlos/'+sim+'veldisplos.txt') # all comp.
                self.sigfiles.append(self.dir+'velocitydispersionlos/'+sim+'veldisplos.txt') # first comp.
            else:
                self.nufiles.append(self.dir+'densityfalloff/' +sim+'falloffnotnorm_'+\
                                    self.nstr1+'_'+self.nstr2+'.txt') # all comp.
                self.nufiles.append(self.dir+'densityfalloff/' +sim+'falloffnotnorm_'+\
                                    self.nstr1+'_'+self.nstr2+'.txt') # first comp.
                self.sigfiles.append(self.dir+'/velocitydispersionlos/' +sim+'veldisplos_'\
                                     +self.nstr1+'_'+self.nstr2+'.txt') # all comp.
                self.sigfiles.append(self.dir+'/velocitydispersionlos/' +sim+'veldisplos_'\
                                     +self.nstr1+'_'+self.nstr2+'.txt') # first comp.
                
        elif gp.pops == 2: # before: _*_[1,2].txt
            self.nufiles.append(self.dir+'densityfalloff/'+sim+'falloffnotnorm.txt') # all comp.
            self.nufiles.append(self.dir+'densityfalloff/'\
                                +sim+'falloffnotnorm_'+self.nstr1+'_0.txt')
            self.sigfiles.append(self.dir+'/velocitydispersionlos/' +sim+'veldisplos_'+self.nstr1+'_'+nstr2+'.txt') # all comp.
            self.sigfiles.append(self.dir+'velocitydispersionlos/'\
                                 +sim+'veldisplos_'+self.nstr1+'_0.txt')
            self.nufiles.append(self.dir+'densityfalloff/'\
                                +sim+'falloffnotnorm_0_'+self.nstr2+'.txt')
            self.sigfiles.append(self.dir+'velocitydispersionlos/'\
                                 +sim+'veldisplos_0_'+self.nstr2+'.txt')
        return





        

    def set_walker(self):
        '''derive filenames from Walker&Penarrubia parameters'''
        self.dir = self.machine + 'data_walker/'
        if gp.walkercase == 0:
            gamma_star1 =   0.1;    gamma_star2 =   1.0 # 1. or 0.1
            beta_star1  =   5.0;    beta_star2  =   5.0 # fixed to 5
            r_star1     = 1000.;    r_star2     = 1000. # 500 or 1000
            r_a1        =   1.0;    r_a2        =   1.0
            gamma_DM    = 0 # 0 or 1

        elif gp.walkercase == 1:
            gamma_star1 =   1.0;    gamma_star2 =   1.0 # 1. or 0.1
            beta_star1  =   5.0;    beta_star2  =   5.0 # fixed to 5
            r_star1     =  500.;    r_star2     = 1000. # 500 or 1000
            r_a1        =   1.0;    r_a2        =   1.0
            gamma_DM    = 0 # core

        elif gp.walkercase == 2:
            gamma_star1 =   1.0;    gamma_star2 =   1.0 # 1. or 0.1
            beta_star1  =   5.0;    beta_star2  =   5.0 # fixed to 5
            r_star1     =  500.;    r_star2     = 1000. # 500 or 1000
            r_a1        =   1.0;    r_a2        =   1.0
            gamma_DM    = 1 # cusp

            
        alpha_DM    = 1;    beta_DM     = 3;
        rho0        = 1 # to be read from the corresp. gs*mpc3.dat file
        r_DM        = 1000                    # fixed to 1000pc
        
        AAA = gh.myfill(100*gamma_star1)     #100
        BBB = gh.myfill(10*beta_star1)       #050
        CCC = gh.myfill(r_star1/10)          #100
        DDD = gh.myfill(100*r_a1)            #100
        EEE = {0: "core",
             1: "cusp"
             }[gamma_DM]                 #core
        FFF = gh.myfill(100*gamma_star2)     #010
        GGG = gh.myfill(10*beta_star2)       #050
        HHH = gh.myfill(r_star2/10)          #100
        III = gh.myfill(100*r_a2)            #100
        JJJ = EEE                         #core
        NNN = gh.myfill(3)                   #003    # realization (1..10)
        
        self.dir = self.dir+ "c1_"+AAA+"_"+BBB+"_"+CCC+"_"+DDD+"_"+EEE+"_c2_"+FFF+"_"+GGG+\
                  "_"+HHH+"_"+III+"_"+JJJ+"_"+NNN+"_6d/"

        self.massfile = self.dir+'enclosedmass/enclosedmass_0.txt'
        self.analytic = self.dir+'samplepars'
        print 'analytic set to ',self.analytic
        
        self.nufiles.append(self.dir+'nu/nunotnorm_0.txt') # all comp.
        self.sigfiles.append(self.dir+'siglos/siglos_0.txt')
        if gp.pops == 1:
            self.nufiles.append(self.dir+'nu/nunotnorm_0.txt') # first and only comp.
            self.sigfiles.append(self.dir+'siglos/siglos_0.txt')
        elif gp.pops == 2:
            self.nufiles.append(self.dir+'nu/nunotnorm_1.txt') # first comp.
            self.sigfiles.append(self.dir+'siglos/siglos_1.txt')
            self.nufiles.append(self.dir+'nu/nunotnorm_2.txt') # second comp.
            self.sigfiles.append(self.dir+'siglos/siglos_2.txt')
        self.outdir, self.outname = self.get_outname()
        return
    



    def set_triaxial(self):
        self.dir = self.machine + '/data_triaxial/'
        if gp.triaxcase == 0:           # core
            casename = 'StarsInCore'
        else:
            casename = 'StarsInCusp'
        if gp.projcase == 1:            # along X
            proj = 'X'
        elif gp.projcase == 2:          # along Y
            proj = 'Y'
        elif gp.projcase == 3:          # along Z
            proj = 'Z'
        elif gp.projcase == 4:          # along intermediate axis
            proj = 'I'
        basename = casename + proj
        self.dir = self.dir + basename + '/'

        pre = self.dir
        self.massfile = pre + 'enclosedmass/enclosedmass_0.txt'
        self.nufiles.append(pre  + 'nu/nunotnorm_0.txt')
        self.sigfiles.append(pre + 'siglos/siglos_0.txt')
        self.nufiles.append(pre  + 'nu/nunotnorm_0.txt') # first and only comp.
        self.sigfiles.append(pre + 'siglos/siglos_0.txt')
        if gp.pops == 2:
            print 'TODO: 2 tracer populations for triaxial dataset'
            pdb.set_trace()
        self.outdir, self.outname = self.get_outname()
        return basename






    def set_fornax(self):
        self.dir = self.machine + '/data_obs/for/'

        self.massfile = self.dir+'enclosedmass.txt'
        self.nufiles.append(self.dir+'densityfalloff.txt')
        self.sigfiles.append(self.dir+'velocitydispersionlos.txt')
        if gp.pops == 1:
            self.nufiles.append(self.dir+'densityfalloff.txt') # first and only comp.
            self.sigfiles.append(self.dir+'velocitydispersionlos.txt')
        if gp.pops == 2:
            self.nufiles.append(self.dir+'densityfalloff_1.txt') # first comp.
            self.sigfiles.append(self.dir+'velocitydispersionlos_1.txt')
            self.nufiles.append(self.dir+'densityfalloff_2.txt') # second comp.
            self.sigfiles.append(self.dir+'velocitydispersionlos_2.txt')
        self.outdir, self.outname = self.get_outname()
        return



    def get_outname(self):
        import datetime
        bname = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        if gp.investigate == 'walker':
            bname = bname + '_case_' + str(gp.walkercase)
            bname = bname + '_' + str(get_case(gp.cas))
        if (gp.gprior>0) : bname = bname + '_gprior'
        if (gp.cprior>=0) : bname = bname + '_cprior'
        
        if gp.mirror  : bname = bname + '_mirr'
        if gp.nulog: bname = bname + '_nulog'
        if gp.denslog: bname = bname + '_denslog'
        if gp.lbprior : bname = bname + '_lb'
        
        if gp.deltaprior: bname = bname + '_delta' 
        bname = bname + '_mslope' if (gp.mprior<0) else bname + '_mconst'
        
        if gp.bprior: bname = bname + '_bprior'
        if gp.investigate == 'hernquist':
            bname = bname+'_'+self.nstr1+'_'+self.nstr2+'_'+str(gp.nipol)
        
        if gp.sprior:    bname = bname + '_sprior'
        if gp.uselike:   bname = bname + '_uselike'
        if gp.constdens: bname = bname + '_constdens'
        if gp.adderrors: bname = bname + '_errors'
        if gp.rprior:    bname = bname + '_rprior'
        if gp.quadratic: bname = bname + '_quad'
        if gp.monotonic: bname = bname + '_mono'
        return self.dir+bname+'/', bname




    def get_scale_file(self,i):
        # get rcore, dens0, dens0pc, totmass, maxvlos from file par_0.txt
        return self.dir+'scale_'+str(i)+'.txt'


    def get_ntracer_file(self,i):
        return self.dir+'ntracer_'+str(i)+'.txt'


    def set_disc_sim(self):
        self.dir = self.machine + 'data_disc_sim/mwhr/'
        # self.posvelfiles.append(self.dir + 'sim/mwhr_r8500_ang'+gp.patch+'_stars.txt') # [PS] perhaps we need all components as the first entry in this list, with convention: 0. all 1. pop, 2. pop, 3. pop = background
        # self.nufiles.append(self.dir + 'nu/mwhr_r8500_ang'+gp.patch+'_falloff_stars.txt') # again all components
        # self.sigfiles.append(self.dir +  'siglos/mwhr_r8500_ang'+gp.patch+'_dispvel_stars.txt') # all comp.
        # self.surfdenfiles.append(self.dir + 'surfden/mwhr_r8500_ang'+gp.patch+'_surfaceden.txt') # overall surface density?

        self.posvelfiles.append(self.dir + 'sim/mwhr_r8500_ang'+gp.patch+'_stars.txt') # first comp.
        self.nufiles.append(self.dir + 'nu/mwhr_r8500_ang'+gp.patch+'_falloff_stars.txt') # first comp
        self.sigfiles.append(self.dir +  'siglos/mwhr_r8500_ang'+gp.patch+'_dispvel_stars.txt') # first comp.
        self.surfdenfiles.append(self.dir + 'surfden/mwhr_r8500_ang'+gp.patch+'_surfaceden.txt') # baryonic surface density 
        
        if gp. pops ==2:
            # TODO: stars?
            self.posvelfiles.append(self.dir + 'sim/mwhr_r8500_ang'+gp.patch+'_dm.txt') # second comp.
            self.nufiles.append(self.dir + 'nu/mwhr_r8500_ang'+gp.patch+'_falloff_dm.txt') # second comp.
            self.sigfiles.append(self.dir +  'siglos/mwhr_r8500_ang'+gp.patch+'_dispvel_dm.txt') # second comp.
            self.surfdenfiles.append(self.dir + 'surfden/mwhr_r8500_ang'+gp.patch+'_surfacedenDM.txt') # DM surface density

        return

    
    def set_disc_simple(self):
        self.dir = self.machine + 'data_disc_simple/'
        return


    def get_outfiles(self):
        pre = self.outdir + self.outname
        outplot = pre + '.png'
        outdat  = pre + '.dat'
        outtxt  = pre + '.txt'
        return outplot, outdat, outtxt

    def get_outpng(self):
        pre = self.outdir + self.outname
        return pre + '.png'

    def get_outdat(self):
        pre = self.outdir + self.outname 
        return pre + '.dat'

    def get_outtxt(self):
        pre = self.outdir + self.outname 
        return pre + '.txt'



    def get_outprofs(self):
        pre = self.outdir + self.outname
        profnus = []; profdeltas = []; profsigs = []
        profM = pre +'.profM'
        profdens = pre + '.profdens'
        profnus.append( pre + '.profnu1')
        profdeltas.append( pre + '.profdelta1')
        profsigs.append( pre + '.profsig1')
        if gp.pops == 2:
            profnus.append( pre + '.profnu2')
            profdeltas.append( pre + '.profdelta2')
            profsigs.append( pre + '.profsig2')
        return profM, profdens, profnus, profdeltas, profsigs
