;============================= setcommon ==============================
;This routine sets up common block definitions for all my IDL sofware.

PRO setcommon

;======================================================================
;Maths,Physics and IDL constants:
COMMON constants,G,Pi,seed,pc,pl,eflag
G = DOUBLE(1.)                  ;Gravitational constant
Pi = DOUBLE(3.141592654)        ;Pi

;Random number generator:
seed = 1001L                    ;Random number seed

;Text output:
pc = 0                          ;Print column
pl = 14                         ;Print line

;Error flagging (0 is no error):
eflag = 0


;======================================================================
;Numerical integration:
COMMON integrate,intpnts
intpnts = 1e2                   ;Number of integration points


;======================================================================
;Interpolation and lookup-tables:
COMMON interpolate,rmeshbinmax,rmaxtablesize,phitol,rmeshpnts,$
  rminlt,rmaxlt,rlookup_table,rminpnt,rmaxpnt,rbuffer,$
  Qmeshbinmax,Qmaxtablesize,Ftol,Qmeshpnts,Qfactor,$
  Qminlt,Qmaxlt,Qlookup_table,Qmaxpnt,Qbuffer
  
;Look up table for f(Q):
Qmeshbinmax = 1e3               ;Max no. of mesh bins per log-interval
Qmaxtablesize = 1e6             ;Max length of lookup table array
Ftol = 1.e-3                    ;Tolerance for adaptive mesh points.
                                ;If gradient in dist. func. > Ftol
                                ;then Qmeshbin -> Qmeshbin*Qmeshbin
Qminlt = DOUBLE(1.e-11)         ;Minimum value of Q in lookup table
Qmaxlt = DOUBLE(0)              ;If zero, this is calculated by 
                                ;the program.
Qlookup_table = dblarr(1)       ;lookup_table array (size set later...)
Qbuffer = 25                    ;Buffer size for Q interpolate array. 
                                ;Large is more accurate but also
                                ;slower. Do NOT make this too small

;Look up table for r(phi):
rmeshbinmax = 1e3               ;Max no. of mesh bins per log-interval
rmaxtablesize = 1e6             ;Max length of lookup table array
phitol = 1.e-3                  ;Tolerance for adaptive mesh points.
rmaxlt = DOUBLE(1.e28)          ;Max and min values of r in lookup table
rminlt = DOUBLE(1.e-12)         ;important that these are large/small enough!
rlookup_table = dblarr(1)       ;lookup_table array (size set later...)
rbuffer = 25                    ;Buffer size for r interpolate array. 
                                ;Large is more accurate but also
                                ;slower. Do NOT make this too small


COMMON vcompfunc,Vran_table,Pran_table,vpnts,vminlt,vfactor
;Lookup up table for the V comparison function:
Vran_table = dblarr(1)          ;Look up table arrays defined later.
Pran_table = dblarr(1)          ;
vpnts = 1e1                     ;Number of points in lookup table
vminlt = DOUBLE(1e-12)          ;Smallest value of v in table
vfactor = DOUBLE(2)             ;Constant difference between comparison
                                ;function and true distribution. Large
                                ;is slow but more accurate.

;======================================================================
;Unit conversion:
COMMON units,Msun,parsec,kiloparsec,kilometer,year,G_si,c_si
Msun = 1.989e30                 ;Solar mass (kg)
parsec = 3.08567802e16          ;Parsec (m)
kiloparsec = 1000*parsec        ;Kiloparsec (m)
kilometer = 1000.               ;Kilometer (m)
year = 60.*60.*24.*365.         ;Year (s)
G_si = 6.672e-11                ;Gravitational constant (SI)
c_si = 3.e8                     ;Speed of light (SI)


;======================================================================
;Some useful double precision numbers:
COMMON numbers,zero,one,two,three,four,five,six,seven,eight,nine,$
  ten,eleven,twelve,thirteen,fourteen,fifteen,sixteen,seventeen
zero = DOUBLE(0)
one = DOUBLE(1)
two = DOUBLE(2)
three = DOUBLE(3)
four = DOUBLE(4)
five = DOUBLE(5)
six = DOUBLE(6)
seven = DOUBLE(7)
eight = DOUBLE(8)
nine = DOUBLE(9)
ten = DOUBLE(10)
eleven = DOUBLE(11)
twelve = DOUBLE(12)
thirteen = DOUBLE(13)
fourteen = DOUBLE(14)
fifteen = DOUBLE(15)
sixteen = DOUBLE(16)
seventeen = DOUBLE(17)

END
