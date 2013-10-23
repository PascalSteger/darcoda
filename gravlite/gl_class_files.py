#!/usr/bin/env python3

##
# @file
# class to generate all filenames, independent of which investigation

# (c) 2013 ETHZ, psteger@phys.ethz.ch

import pdb
import gl_params as gp
import gl_helper as gh

def get_case(cas):
    ## Set number of tracers to look at
    # want to set ntracer = 3e3              # case 1
    #             ntracer = 3e4              # case 2
    
    ntracer = 0
    if cas == 1:
        ntracer = 3000
    elif cas == 2:
        ntracer = 10000
    return ntracer

## Common base class for all filename sets
class Files:
    # massfile, nufile[], sigfile[], outname = Files(investigate)
    
    ## constructor
    def __init__ (self):
        ## set which computer we are working on
        self.machine = ''
        ## set base directory
        self.dir = ''
        ## relative path to the 'programs' directory
        self.progdir = ''
        self.set_dir(gp.machine) # 'darkside' or 'local'
        ## file with 2D summed masses
        self.massfile = ''
        ## file with analytic values for Walker models
        self.analytic = ''
        ## file with 2D surface density
        self.surfdenfiles = []
        ## files with tracer densities
        self.nufiles  = []
        ## files with velocity dispersions
        self.sigfiles = []
        ## files with centered positions and velocities for the tracer particles
        self.posvelfiles = []
        ## files for the fourth order moment of the LOS velocity
        self.kappafiles = [];
        ## number of tracers
        self.ntracer, self.nstr1 = self.set_ntracer(gp.cas)
        ## TODO: unknown parameters
        self.params = []
        
        if gp.investigate == 'hernquist':
            self.set_hernquist()
        elif gp.investigate == 'walker':
            self.set_walker()
        elif gp.investigate == 'gaia':
            self.set_gaia()
        elif gp.investigate == 'triaxial':
            self.set_triaxial()
        elif gp.investigate == 'fornax':
            self.set_fornax()
        elif gp.investigate == 'sim':
            self.set_disc_sim()
        elif gp.investigate == 'simple':
            self.set_disc_simple()

        ## directory and basename of all output files
        self.outdir, self.outname = self.get_outname()
        return

    ## set local directory
    # depending on which computer is used
    def set_dir(self, machine):
        if machine == 'darkside':
            self.machine = '/home/ast/read/dark/dwarf_data/'
        elif machine == 'local':
            self.machine = '/home/psteger/sci/dwarf_data/'
        self.progdir = self.machine + 'programs'
        return

    ## set number of tracers
    # based on the case (3k, 30k tracers) we are working on
    def set_ntracer(self,cas):
        ntracer = get_case(cas)
        ## number of tracers
        self.ntracer = ntracer
        ## string of the same quantity
        self.nstr = str(ntracer)
        return self.ntracer, self.nstr
        
    ## get simulation names in the Hernquist case
    def get_sim_name(self):
        if gp.pops == 1:
            sim = 'unit_hern_1_'
        elif gp.pops == 2:
            sim = 'dual_unit_hern_1_'
        return sim

    ## set all filenames for Hernquist case
    def set_hernquist(self):
        self.dir = self.machine + 'data_hernquist/'
        sim = self.get_sim_name()
        self.massfile = self.dir+'enclosedmass/'+sim+'enclosedmass_0.txt'
        self.nufiles.append(self.dir+'densityfalloff/'+sim+'falloffnotnorm_0.txt') # all comp.
        self.nufiles.append(self.dir+'densityfalloff/'+sim+'falloffnotnorm_1.txt') # first comp.
        self.sigfiles.append(self.dir+'velocitydispersionlos/'+sim+'veldisplos_0.txt') # all comp.
        self.sigfiles.append(self.dir+'velocitydispersionlos/'+sim+'veldisplos_1.txt') # first comp.
        self.kappafiles.append(self.dir+'velocitykurtosislos/'+sim+'kappalos_0.txt') # all comp.
        self.kappafiles.append(self.dir+'velocitykurtosislos/'+sim+'kappalos_1.txt') # first comp.
        if gp.pops == 2:
            self.nufiles.append(self.dir+'densityfalloff/' +sim+'falloffnotnorm_2.txt')
            self.sigfiles.append(self.dir+'/velocitydispersionlos/' +sim+'veldisplos_2.txt')
            self.kappafiles.append(self.dir+'/kappalos/' +sim+'kappalos_2.txt')

                
        elif gp.pops == 2: # before: _*_[1,2].txt
            self.nufiles.append(self.dir+'densityfalloff/'+sim+'falloffnotnorm.txt') # all comp.
            self.nufiles.append(self.dir+'densityfalloff/'\
                                +sim+'falloffnotnorm_'+self.nstr1+'_0.txt')
            self.sigfiles.append(self.dir+'/velocitydispersionlos/' +sim+'veldisplos_'+self.nstr1+'_'+nstr2+'.txt') # all comp.
            self.sigfiles.append(self.dir+'velocitydispersionlos/'\
                                 +sim+'veldisplos_'+self.nstr1+'_0.txt')
            self.kappafiles.append(self.dir+'kappalos/'\
                                 +sim+'kappalos_'+self.nstr1+'_0.txt')
            self.nufiles.append(self.dir+'densityfalloff/'\
                                +sim+'falloffnotnorm_0_'+self.nstr2+'.txt')
            self.sigfiles.append(self.dir+'velocitydispersionlos/'\
                                 +sim+'veldisplos_0_'+self.nstr2+'.txt')
            self.kappafiles.append(self.dir+'kappalos/'\
                                 +sim+'kappalos_0_'+self.nstr2+'.txt')
        return


    ## derive filenames from gaia case
    def set_gaia(self):
        self.dir = self.machine + 'data_gaia/'
        beta_star1 = 5; r_DM = 1000
        if gp.case == 1:
            gamma_star1=0.1;r_star1=100;r_a1=100;gamma_DM=1;rho0=0.064
        elif gp.case == 2:
            gamma_star1=0.1;r_star1=250;r_a1=250;gamma_DM=0;rho0=0.400
        elif gp.case == 3:
            gamma_star1=0.1;r_star1=250;r_a1=np.inf;gamma_DM=1;rho0=0.064
        elif gp.case == 4:
            gamma_star1=0.1;r_star1=1000;r_a1=np.inf;gamma_DM=0;rho0=0.400
        elif gp.case == 5:
            gamma_star1=1.0;r_star1=100;r_a1=100;gamma_DM=1;rho0=0.064
        elif gp.case == 6:
            gamma_star1=1.0;r_star1=250;r_a1=250;gamma_DM=0;rho0=0.400
        elif gp.case == 7:
            gamma_star1=1.0;r_star1=250;r_a1=np.inf;gamma_DM=1;rho0=0.064
        elif gp.case == 8:
            gamma_star1=1.0;r_star1=1000;r_a1=np.inf;gamma_DM=0;rho0=0.400

        self.params = [beta_star1, r_DM, gamma_star1, r_star1,r_a1,gamma_DM,rho0]
            
        AAA = gh.myfill(100*gamma_star1)     #100
        BBB = gh.myfill(10*beta_star1)       #050
        CCC = gh.myfill(100*r_star1/r_DM)    #100
        DDD = gh.myfill(100*r_a1/r_star1)     #100
        EEEE = {0: "core",                     
             1: "cusp"
             }[gamma_DM]                 #core
        FFFF = gh.myfill(1000*rho0,4)
        self.dir = self.dir+'gs'+AAA+'_bs'+BBB+'_rcrs'+CCC+'_rarc'\
            +DDD+'_'+EEEE+'_'+FFFF+'mpc3_df/'
        ## new variable to hold the .dat input file
        self.datafile = self.dir + 'dat'
        
        self.massfile = self.dir+'enclosedmass/enclosedmass_0.txt'
        self.nufiles.append(self.dir+'nu/nunotnorm_0.txt') # all comp.
        self.sigfiles.append(self.dir+'siglos/siglos_0.txt')
        self.kappafiles.append(self.dir+'kappalos/kappalos_0.txt')
        self.nufiles.append(self.dir+'nu/nunotnorm_0.txt') # first and only
        self.sigfiles.append(self.dir+'siglos/siglos_0.txt')
        self.kappafiles.append(self.dir+'kappalos/kappalos_0.txt')

        
    ## derive filenames from Walker&Penarrubia parameters
    def set_walker(self):
        # TODO: include parameter storage in self.params
        self.dir = self.machine + 'data_walker/'
        if gp.case == 0:
            gamma_star1 =   0.1;    gamma_star2 =   1.0 # 1. or 0.1
            beta_star1  =   5.0;    beta_star2  =   5.0 # fixed to 5
            r_star1     = 1000.;    r_star2     = 1000. # 500 or 1000
            r_a1        =   1.0;    r_a2        =   1.0
            gamma_DM    = 0 # 0 or 1

        elif gp.case == 1:
            gamma_star1 =   1.0;    gamma_star2 =   1.0 # 1. or 0.1
            beta_star1  =   5.0;    beta_star2  =   5.0 # fixed to 5
            r_star1     =  500.;    r_star2     = 1000. # 500 or 1000
            r_a1        =   1.0;    r_a2        =   1.0
            gamma_DM    = 0 # core

        elif gp.case == 2:
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
        print('analytic set to ',self.analytic)
        
        self.nufiles.append(self.dir+'nu/nunotnorm_0.txt') # all comp.
        self.sigfiles.append(self.dir+'siglos/siglos_0.txt')
        self.kappafiles.append(self.dir+'kappalos/kappalos_0.txt')
        if gp.pops == 1:
            self.nufiles.append(self.dir+'nu/nunotnorm_0.txt') # first and only comp.
            self.sigfiles.append(self.dir+'siglos/siglos_0.txt')
            self.kappafiles.append(self.dir+'kappalos/kappalos_0.txt')
        elif gp.pops == 2:
            self.nufiles.append(self.dir+'nu/nunotnorm_1.txt') # first comp.
            self.sigfiles.append(self.dir+'siglos/siglos_1.txt')
            self.kappafiles.append(self.dir+'kappalos/kappalos_1.txt')
            self.nufiles.append(self.dir+'nu/nunotnorm_2.txt') # second comp.
            self.sigfiles.append(self.dir+'siglos/siglos_2.txt')
            self.kappafiles.append(self.dir+'kappalos/kappalos_2.txt')
        self.outdir, self.outname = self.get_outname()
        return
    
    ## set all parameters for working on the triaxial data
    def set_triaxial(self):
        self.dir = self.machine + '/data_triaxial/'
        if gp.case == 0:           # core
            casename = 'StarsInCore'
        elif gp.case == 1:
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
        self.kappafiles.append(pre + 'kappalos/kappalos_0.txt')
        self.nufiles.append(pre  + 'nu/nunotnorm_0.txt') # first and only comp.
        self.sigfiles.append(pre + 'siglos/siglos_0.txt')
        self.kappafiles.append(pre + 'kappalos/kappalos_0.txt')
        if gp.pops == 2:
            print('TODO: 2 tracer populations for triaxial dataset')
            pdb.set_trace()
        self.outdir, self.outname = self.get_outname()
        return basename

    ## set all variables in the case we work with Fornax observational data
    def set_fornax(self):
        self.dir = self.machine + '/data_obs/for/'

        self.massfile = self.dir+'enclosedmass.txt'
        self.nufiles.append(self.dir+'densityfalloff.txt')
        self.sigfiles.append(self.dir+'velocitydispersionlos.txt')
        self.kappafiles.append(self.dir+'kappalos.txt')
        if gp.pops == 1:
            self.nufiles.append(self.dir+'densityfalloff.txt') # first and only comp.
            self.sigfiles.append(self.dir+'velocitydispersionlos.txt')
            self.kappafiles.append(self.dir+'kappalos.txt')
        if gp.pops == 2:
            self.nufiles.append(self.dir+'densityfalloff_1.txt') # first comp.
            self.sigfiles.append(self.dir+'velocitydispersionlos_1.txt')
            self.kappafiles.append(self.dir+'kappalos_1.txt')
            self.nufiles.append(self.dir+'densityfalloff_2.txt') # second comp.
            self.sigfiles.append(self.dir+'velocitydispersionlos_2.txt')
            self.kappafiles.append(self.dir+'kappalos_2.txt')
        self.outdir, self.outname = self.get_outname()
        return

    ## determine output directory and filenames
    def get_outname(self):    
        import datetime
        bname = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        if gp.investigate == 'walker':
            bname = bname + '_case_' + str(gp.case)
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
            bname = bname+'_'+self.nstr1+'_'+str(gp.nipol)
        
        if gp.sprior:    bname = bname + '_sprior'
        if gp.uselike:   bname = bname + '_uselike'
        if gp.constdens: bname = bname + '_constdens'
        if gp.adderrors: bname = bname + '_errors'
        if gp.rprior:    bname = bname + '_rprior'
        if gp.quadratic: bname = bname + '_quad'
        if gp.monotonic: bname = bname + '_mono'
        return self.dir+bname+'/', bname

    ## return filename for storing scaling properties (half-mass radius, ...)
    def get_scale_file(self,i):
        return self.dir+'scale_'+str(i)+'.txt'

    ## get filename with attached tracer information
    def get_ntracer_file(self,i):
        return self.dir+'ntracer_'+str(i)+'.txt'

    ## set all properties for disc case        
    def set_disc_sim(self):
        self.dir = self.machine + 'data_disc_sim/mwhr/'
        # self.posvelfiles.append(self.dir + 'sim/mwhr_r8500_ang'+gp.patch+'_stars.txt') # TODO: we need all components as the first entry in this list, with convention: 0. all 1. pop, 2. pop, 3. pop = background
        # self.nufiles.append(self.dir + 'nu/mwhr_r8500_ang'+gp.patch+'_falloff_stars.txt') # again all components
        # self.sigfiles.append(self.dir +  'siglos/mwhr_r8500_ang'+gp.patch+'_dispvel_stars.txt') # all comp.
        # self.surfdenfiles.append(self.dir + 'surfden/mwhr_r8500_ang'+gp.patch+'_surfaceden.txt') # overall surface density?

        self.posvelfiles.append(self.dir + 'sim/mwhr_r8500_ang'+gp.patch+'_stars.txt') # first comp.
        self.nufiles.append(self.dir + 'nu/mwhr_r8500_ang'+gp.patch+'_falloff_stars.txt') # first comp
        self.sigfiles.append(self.dir +  'siglos/mwhr_r8500_ang'+gp.patch+'_dispvel_stars.txt') # first comp.
        self.kappafiles.append(self.dir +  'kappalos/mwhr_r8500_ang'+gp.patch+'_kappa_stars.txt') # first comp.
        self.surfdenfiles.append(self.dir + 'surfden/mwhr_r8500_ang'+gp.patch+'_surfaceden.txt') # baryonic surface density 
        
        if gp. pops ==2:
            # TODO: stars?
            self.posvelfiles.append(self.dir + 'sim/mwhr_r8500_ang'+gp.patch+'_dm.txt') # second comp.
            self.nufiles.append(self.dir + 'nu/mwhr_r8500_ang'+gp.patch+'_falloff_dm.txt') # second comp.
            self.sigfiles.append(self.dir +  'siglos/mwhr_r8500_ang'+gp.patch+'_dispvel_dm.txt') # second comp.
            self.kappafiles.append(self.dir +  'kappalos/mwhr_r8500_ang'+gp.patch+'_kappa_dm.txt') # second comp.
            self.surfdenfiles.append(self.dir + 'surfden/mwhr_r8500_ang'+gp.patch+'_surfacedenDM.txt') # DM surface density

        return

    ## set all properties if looking at simple disc (generated on the fly, no data input files needed)
    def set_disc_simple(self):
        self.dir = self.machine + 'data_disc_simple/'
        return

    ## get all output filenames
    def get_outfiles(self):
        pre = self.outdir # + self.outname + '.'
        outplot = pre + 'png'
        outdat  = pre + 'dat'
        outtxt  = pre + 'txt'
        return outplot, outdat, outtxt

    ## get output filename for png picture
    def get_outpng(self):
        pre = self.outdir # + self.outname + '.'
        return pre + 'png'

    ## get output filename for data
    def get_outdat(self):
        pre = self.outdir # + self.outname + '.'
        return pre + 'dat'

    ## get output filename for ASCII output
    def get_outtxt(self):
        pre = self.outdir # + self.outname + '.'
        return pre + 'txt'

    ## get output filenames for profiles
    def get_outprofs(self):
        pre = self.outdir # + self.outname + '.'
        profnus = []; profdeltas = []; profsigs = []; profkaps = []
        profM = pre +'profM'
        profdens = pre + 'profdens'
        profnus.append( pre + 'profnu1')
        profdeltas.append( pre + 'profdelta1')
        profsigs.append( pre + 'profsig1')
        profkaps.append( pre + 'profkap1')
        if gp.pops == 2:
            profnus.append( pre + 'profnu2')
            profdeltas.append( pre + 'profdelta2')
            profsigs.append( pre + 'profsig2')
            profkaps.append( pre + 'profkap2')
        return profM, profdens, profnus, profdeltas, profsigs, profkaps

# Local Variables:
# py-master-file: "gravlite.py"
# End:

