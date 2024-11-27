import numpy as np
import pandas as pd
from scipy.constants import c as vlight
from scipy.integrate import quad
from scipy.interpolate import CubicSpline

from scipy.misc import derivative
import camb
from scipy.interpolate import interp1d

# now in km/s
vlight = vlight / 1000

class Redshift_Remapping:
    def __init__(self, omega_m, omega_lambda, w0, wa, H0, a, b, c, d, der0, der1):
        self.omega_m = omega_m
        self.omega_lambda = omega_lambda
        self.H0 = H0
        
        self.omegak = 1 - (self.omega_m + self.omega_lambda)
        
        self.z0 = 0
        self.z1 = 0.5
        self.z2 = 1
        self.z3 = 2.5
        
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        
        self.der0 = der0
        self.der1 = der1
        
        self.w0 = w0
        self.wa = wa

    def alfa(self, zobs):
        x = np.array([self.z0, self.z1, self.z2, self.z3])
        y = np.array([self.a, self.b, self.c, self.d])
        
        slope_start = self.der0
        slope_end = self.der1
        
        cs = CubicSpline(x, y, bc_type=((1, slope_start), (1, slope_end)))
        return cs(zobs)
    
    def alfap(self, zobs):
        alfa_z = lambda z: self.alfa(z)
        derivative_alfa = derivative(alfa_z, zobs, dx=1e-11)
        return derivative_alfa
    
    def theta(self, zobs):
        den = 1 + zobs
        num = 1 + zobs + self.alfa(zobs) * zobs
        return np.log(num / den)

    def zflrw(self, zobs):
        return zobs * (1 + self.alfa(zobs))

    def E(self, x):
        z = 1 + x
        pt1 = self.omega_m * z**3
        pt2 = self.omegak * z**2
        
        a = 1 / (1 + x)
        pt3 = self.omega_lambda * pow(a, -3 * (1 + self.w0 + self.wa)) * np.exp(-3 * self.wa * (1 - a))
        return np.sqrt(pt1 + pt2 + pt3)
    
    def H(self, x):
        return self.H0 * self.E(x)


    def DC(self, zobs):
        integrand = lambda x: 1 / self.E(x)
        upplim = zobs * (1 + self.alfa(zobs))
        res = []
        for i in range(len(upplim)):
            res.append(quad(integrand, 0, upplim[i])[0])
        res = np.array(res)
        return self.DH() * res

    def DH(self):
        return vlight / (self.H0)

    def DM(self, zobs):
        DC = self.DC(zobs)
        DH = self.DH()
        if self.omegak > 0:
            a = DH / np.sqrt(self.omegak)
            b = np.sqrt(self.omegak) * DC / DH
            return a * np.sinh(b)
        elif self.omegak == 0:
            return DC
        else:
            a = DH / np.sqrt(abs(self.omegak))
            b = np.sqrt(abs(self.omegak)) * DC / DH
            return a * np.sin(b)

    def DL(self, zobs):
        return (1 + zobs) * self.DM(zobs)


    def magnitude(self, zobs):
        x=self.DL(zobs)
        return 5*np.log10(x)+25

