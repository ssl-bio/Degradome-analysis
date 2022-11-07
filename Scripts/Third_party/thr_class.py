from scipy.special import gammaincc,gamma,gammaln
from scipy.optimize import brentq
from scipy.integrate import quad
from scipy.stats import norm
from scipy import interpolate
from math import *
import struct
import numpy as np
import sys

# Parameter setting the transition in the scaling function from linear to log
lam=10
flam=float(lam)

# Scaling function used for binning data used in the regression to obtain the
# prior distribution
def l_scale(x):
    if x <= lam:
        return x
    else:
        return int(flam*(1.+log(x/flam))+0.5)

# Inverse scaling function
def il_scale(x):
    if x <= lam:
        return x
    else:
        return flam*exp((float(x)-flam)/flam)

# Bin size in the scaling function
def l_size(x):
    if x < lam:
        return 1.
    else:
        return float(int(il_scale(x+0.5))-int(il_scale(x-0.5)))

# Inverse cumulative normal distribution
def phi(x):
    return (1.0+erf(x/sqrt(2.0)))/2.0

# Class for computing tables of read count thresholds, used in the Bayesian
# analysis to determine significant differences between the control and test
# read counts
class thresh_tab:
    def __init__(self):
        pass

    # Loads in a previously computed threshold table
    def load_table(self,filename):
        fi=open(filename,'r')

        # Read in header
        fi.readline()
        self.mdl_a=float(fi.readline())
        self.mdl_b=float(fi.readline())
        self.mdl_c=float(fi.readline())
        self.mdl_d=float(fi.readline())
        fi.readline()
        fi.readline()
        self.pthresh=float(fi.readline())
        self.psi=norm.ppf(self.pthresh)
        fi.readline()
        fi.readline()

        # Read in table dimensions
        te=fi.readline().split(' ')
        self.fs=int(te[0])
        self.gs=int(te[1])
        fs=self.fs
        gs=self.gs
        fi.readline()
        fi.readline()

        # Read in the zero cut information
        te=fi.readline().split(' ')
        self.nFz=int(te[0])
        self.Fcut=float(te[1])
        if self.nFz>0:
            self.Fxcut=np.zeros((self.nFz+1))
        fi.readline()
        fi.readline()

        # Read in the ratio data points
        self.Fl=np.zeros((fs))
        i=0
        while i<self.nFz:
            te=fi.readline().split(' ')
            self.Fl[i]=float(te[0])
            self.Fxcut[i]=float(te[1])
            i+=1
        while i<fs:
            self.Fl[i]=float(fi.readline())
            i+=1
        fi.readline()
        fi.readline()

        # Read in the control count data points
        self.xl=np.zeros((gs))
        for i in range(gs):
            self.xl[i]=float(fi.readline())
        fi.readline()
        fi.readline()

        # Read in the matrix of thresholds
        self.tha=np.zeros((fs,gs))
        for i in range(fs):
            self.tha[i,:]=np.array([fi.readline().split(' ')])
        fi.close()

        self.make_interp()

    # Constructs a bicubic interpolation of threshold values from the
    # sampled grid of thresholds
    def make_interp(self):

        if self.nFz>0:
            self.xz_inter=interpolate.interp1d(self.Fl[:self.nFz+1],self.Fxcut,kind='cubic')
            self.fxz_inter=interpolate.interp2d(self.xl,self.Fl[:self.nFz+1],self.tha[:self.nFz+1,:],kind='cubic')
        self.f_inter=interpolate.interp2d(self.xl,self.Fl[self.nFz:],self.tha[self.nFz:,:],kind='cubic')

    # Saves the table of thresholds
    def save_table(self,filename):

        # Open table and write header
        fi=open(filename,'w')
        fi.write("# Model parameters\n")
        fi.write(str(self.mdl_a)+"\n"+str(self.mdl_b)+"\n"+str(self.mdl_c)+"\n"+str(self.mdl_d))
        fi.write("\n\n# Confidence level\n"+str(self.pthresh))
        fi.write("\n\n# Table dimensions\n"+str(self.fs)+" "+str(self.gs))
        fi.write("\n\n# Points and cutoff for zero region interpolation\n"+str(self.nFz)+" "+str(self.Fcut))

        # Output data
        fi.write("\n\n# Log ratio data points (plus optional zero information)\n")
        i=0
        while i<self.nFz:
            fi.write(str(self.Fl[i])+" "+str(self.Fxcut[i])+"\n")
            i+=1
        while i<self.fs:
            fi.write(str(self.Fl[i])+"\n")
            i+=1
        fi.write("\n# Control count data points\n")
        for i in range(self.gs):
            fi.write(str(self.xl[i])+"\n")
        fi.write("\n# Thresholds\n")
        for i in range(self.fs):
            fi.write(" ".join([str(self.tha[i,j]) for j in range(self.gs)])+"\n")
        fi.close()

    # Computes the threshold for a given test count x, and multiplicative
    # factor F
    def thr(self,x,F):

        # Check that F inputs are within the range of the table
        l=exp(self.Fl[0]);u=exp(self.Fl[-1])
        if F<l or F>u:
            print "Ratio of",F,"out of table range from",l,"to",u
            sys.exit()

        # Search for case where the (x,F) pair is in the zero region
        lF=log(F)
        if F<self.Fcut:
            xz=self.xz_inter(lF)
            if x<xz:
                return 0
            xi=x-xz
            if xi>self.xl[-1]:
                return self.ycr_asymp(x,F)
            return self.fxz_inter(xi,lF)[0]+x*F
        else:
            xi=x
            if xi>self.xl[-1]:
                return self.ycr_asymp(x,F)
            return self.f_inter(xi,lF)[0]+x*F

    # Upper incomplete Gamma function
    def mygi(self,s,x):
        if s<0:
            return (self.mygi(s+1,x)-x**s*exp(-x))/s
        else:
            return gammaincc(s,x)*gamma(s)

    # Normalized upper incomplete Gamma function
    def myngi(self,s,x):
        if s<0:
            return self.myngi(s+1,x)-x**s*exp(-x)/gamma(s+1)
        else:
            return gammaincc(s,x)

    # Control function to rootfind on
    def rootf(self,x):
        a=self.expo
        return (a-1)*x**(a-1)*self.mygi(1-a,x)-self.prob0

    # Fit prior probability based on the read count distribution. chd_l[0] is
    # the control count distribution, and chd_l[1] is the test count
    # distribution.
    def fit_prior(self,chd_l,tot_sites=32348300.):
        na=["control","test"]

        # Fit two priors: 0=control and 1=test
        for i in [0,1]:
            cc=sum(chd_l[i])

            # Give error in case when more non-zero read counts have been
            # observed than there are possible sites
            if cc>tot_sites:
                sys.stderr.write("Error: total non-zero read counts in "+na[i]+" data set\nexceed estimated total number of possible sites\n")
                sys.exit()

            # Perform linear regression on the log-log data
            s=0.
            sx=0.;sy=0.
            sxx=0.;sxy=0.;syy=0.
            nbin=0
            ccou=0
            csz=0.
            clx=0.
            for j in range(1,len(chd_l[i])):
                nbin+=1
                ccou+=chd_l[i][j]
                csz+=l_size(j)
                clx+=log(il_scale(j))
                if ccou>4:
                    x=clx/nbin
                    y=log(float(ccou)/csz)
                    nbin=0
                    ccou=0
                    csz=0.
                    clx=0.

                    s+=1;sx+=x;sy+=y
                    sxx+=x*x;sxy+=x*y;syy+=y*y
            ni=1.0/s
            sx*=ni;sy*=ni
            sxx=sxx*ni-sx*sx
            sxy=sxy*ni-sx*sy
            syy=syy*ni-sy*sy

            # Print status message about the goodness of fit, which should be
            # close to 1 if the data fits a power law well
            sys.stderr.write("%s data set:\nr = %g\n" % (na[i],sxy/sqrt(sxx*syy)))

            # Determine the power-law slope, and determine the lower cutoff
            # based on the probability of seeing a zero count
            self.expo=-sxy/sxx
            self.prob0=1-float(cc)/tot_sites
            cut=brentq(self.rootf,1e-6,1)

            # Print status messages
            if i==0:
                self.mdl_a=self.expo
                self.mdl_b=cut
                sys.stderr.write("a = %g\nb = %g\n\n" % (self.mdl_a,self.mdl_b))
            else:
                self.mdl_c=self.expo
                self.mdl_d=cut
                sys.stderr.write("c = %g\nd = %g\n\n" % (self.mdl_c,self.mdl_d))

    # Sets the parameters in the prior distribution
    def set_prior(self,a,tot_a,c,tot_c,tot_sites=32348300.):
        self.expo=a
        self.prob0=1.-tot_a/tot_sites
        self.mdl_a=a
        self.mdl_b=brentq(self.rootf,1e-6,1)
        self.expo=c
        self.prob0=1.-tot_c/tot_sites
        self.mdl_c=c
        self.mdl_d=brentq(self.rootf,1e-6,1)
        sys.stderr.write("a = %g\nb = %g\n\n" % (self.mdl_a,self.mdl_b))
        sys.stderr.write("c = %g\nd = %g\n\n" % (self.mdl_c,self.mdl_d))

    # Make the table
    def make_table(self,pthresh):
        self.pthresh=pthresh
        self.psi=norm.ppf(pthresh)
        self.xl=[0,0.25,0.5]+[int(399940.*i*i*i*i*i/60.**5+i) for i in range(1,61)]
        self.gs=len(self.xl)
        self.fs=41
        Flo=0.001
        Fhi=1000
        lFlo=log(Flo)
        lFhi=log(Fhi)
        self.tha=np.zeros((self.fs,self.gs))

        # Calculate the zero region (if any)
        if self.prob(0,0,Flo)>self.pthresh:

            if self.prob(0,0,Fhi)>self.pthresh:
                sys.stderr.write("Characteristics of zero region invalid\n")
                sys.exit()

            # Use a bisection search to find the value of F where ycr=0 for x=0
            lo=Flo
            hi=Fhi
            while hi-lo>0.5e-10:
                mid=0.5*(hi+lo)
                if self.prob(0,0,mid)>self.pthresh:
                    lo=mid
                else:
                    hi=mid
            self.Fcut=0.5*(lo+hi)
            lFcut=log(self.Fcut)
            self.nFz=int((self.fs-2)*(lFcut-lFlo)/(lFhi-lFlo))+1
            self.Fl=np.concatenate((np.linspace(lFlo,lFcut,self.nFz,False),np.linspace(lFcut,lFhi,self.fs-self.nFz)))

            # For values of F lower than the critical value, find the values of
            # x such that ycr=0, marking the boundary of the zero region
            self.Fxcut=np.zeros((self.nFz+1))
            for i in range(0,self.nFz):
                lo=0
                hi=1
                F=exp(self.Fl[i])
                while(self.prob(hi,0,F)>pthresh): hi*=2
                while hi-lo>0.5e-10:
                    mid=0.5*(hi+lo)
                    if self.prob(mid,0,F)>self.pthresh:
                        lo=mid
                    else:
                        hi=mid
                self.Fxcut[i]=0.5*(lo+hi)

        else:

            # There is no zero region with the table range. Set the zero region
            # cutoff to be below the lower limit of F values, so it is never
            # triggered
            self.Fl=np.linspace(lFlo,lFhi,self.fs)
            self.Fcut=0.9*self.Fl[0]
            self.nFz=0

        # Calculate the thresholds
        for i in range(self.fs):
            F=exp(self.Fl[i])
            if i<self.nFz:
                xz=self.Fxcut[i]
            else:
                xz=0
            fac=1
            for j in range(self.gs):
                x=xz+self.xl[j]*fac
                self.tha[i,j]=self.ycr(x,F)-x*F
        self.make_interp()

    # Function to calculate the critical value of y, using a bisection search
    # on the probability function
    def ycr(self,x,F):
        if self.prob(x,0,F)>self.pthresh: return 0.0
        else:

            # Use a bisection search to find the threshold
            lo=0.
            hi=1.+1.1*F*x
            while self.prob(x,hi,F)<self.pthresh:
                lo=hi
                hi*=2.0
            er=1e-10*hi
            while hi-lo>er:
                mid=0.5*(hi+lo)
                if self.prob(x,mid,F)<self.pthresh:
                    lo=mid
                else:
                    hi=mid
            return 0.5*(lo+hi)

    # Asymptotic formula for the critical value of y, based on approximating
    # the posterior distributions as gamma distributions
    def ycr_asymp(self,x,F):
        al=x-self.mdl_a
        psisq=self.psi*self.psi
        return self.mdl_c+(2*al*F+psisq+self.psi*sqrt(4*al*F*(1+F)+psisq))/2.

    # Probability function
    def prob(self,x,y,F):

        # Get model constants
        a=self.mdl_a;b=self.mdl_b
        c=self.mdl_c;d=self.mdl_d

        # Set constants that appear in the integrand
        self.opxma=1.+x-a
        self.ymc=y-c
        self.gln=gammaln(1+y-c)
        self.igf=1./self.myngi(1+x-a,b)
        self.Fi=1./F

        # Estimate an integration interval to use, where the bulk of the non-zero
        # probability density is
        dp=max(d,b*F)
        ra=25.*sqrt(y)
        lo=y-ra
        hi=y+ra
        if lo<dp:
            lo=dp
        if hi<100.:
            hi=100.
        mi=(lo+hi)*0.5

        # Compute the integral as the sum of two sub-intevals surrounding the mode
        tol=1e-9*x*F
        if tol<1.49e-8: tol=1.49e-8
        (ti,ep)=quad(self.fu,lo,mi,limit=1000,epsabs=tol)
        (ti2,ep)=quad(self.fu,mi,hi,limit=1000,epsabs=tol)
        return (ti+ti2)/abs(self.myngi(1.+y-c,d))

    # Function to integrate
    def fu(self,th):
        return exp(-th+log(th)*self.ymc-self.gln)*(1.-self.myngi(self.opxma,th*self.Fi)*self.igf)

        # Test code for approximating integrand as a multivariate gaussian
        # al=self.opxma-1
        # be=self.ymc
        # if al<=0.01: al=0.01
        # if be<=0.01: be=0.01
        # return 1/sqrt(2*pi*be)*exp(-(th-be)*(th-be)/(2.*be))*phi((th*self.Fi-al)/sqrt(al))
