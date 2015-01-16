from scipy import *
from numpy import *
from matplotlib import *
from matplotlib.pyplot import *
from matplotlib.axis import *
from scipy.optimize import leastsq

def residuals(p, y, x): 
	err = y-peval(x,p) 
	return err

def peval(x, p): 
	return p[0]*(1-exp(-(p[2]*x)**p[4])) + p[1]*(1-exp(-(p[3]*(x))**p[5] ))

filename=('tgdata.dat')
data = loadtxt(filename)

y = data[:,1]
x = data[:,0]

A1_0=4
A2_0=3
k1_0=0.5
k2_0=0.04
n1_0=2
n2_0=1
pname = (['A1','A2','k1','k2','n1','n2'])
p0 = array([A1_0 , A2_0, k1_0, k2_0,n1_0,n2_0])
plsq = leastsq(residuals, p0, args=(y, x), maxfev=2000)
plot(x,y,'x')#,'title "Meas" with points')
plot(x,peval(x,plsq[0]))#,'title "Fit" with lines lt -1')
#gca().yaxis((0, 7))
legend('right bottom Left')
#xtitle('Time [h]')
#ytitle('Hydrogen release [wt. %]')
grid("off")
savefig('kinfit.png')

print "Final parameters"
for i in range(len(pname)):
	print "%s = %.4f " % (pname[i], p0[i])
