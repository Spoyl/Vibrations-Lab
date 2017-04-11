import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt
import cmath as cm

flr_brdth = 0.199
flr_wdth = 0.009
density = 2720
E = 68.9*10**9
hv = [0.086, 0.063]
b = 0.025
d = 0.0039
m1 = 1.32
m2 = 0.97
clmn_wdth = 0.025
clmn_hght = 0.2

Kc = 12772.12
k_array = np.array([[Kc, Kc],[Kc, Kc]])
k_factors = np.array([[8, -4],[-4, 4]])
m_array = np.array([[m1, 0],[0, m2]])
f_array = np.arange(0, 100, 2)
lossf = 0.01

m_balance = 0.0012
eccen = 0.008

X1 = []
X2 = []
phase1 = []
phase2 = []
def flr_mass(b, w, h, rho):
    return round(rho*b*w*h, 3)

def I(b, d):
    return (b*d**3.)/12.

def K(E, I, l):
    return 12.*(E*I)/(l**3.)

def IL(resp1, resp2):
    return 20.*np.log10(resp1/resp2) 

def K_Finder(b, d, h):
    inertia = I(b, d)
    dis = K(E, inertia, h)
    return dis

def absrbr_f(b, d, hv, m):
    for i in hv:
        kh = K_Finder(b, d, i)
        wn = np.sqrt(kh/m)/(2*np.pi)
        print(wn)

def omega_n(k_array, k_factors, m_array):
    k = k_array*k_factors
    wn2, B = la.eig(k, m_array)
    return (wn2**0.5)/(2*np.pi), B

def MAC():
    mode1 = np.array([[0.747], [-0.664]])
    response1 = np.array([39.4, -15.8])
    mac = (abs(np.dot(mode1.transpose(), response1))**2)/(
                (np.dot(mode1.transpose(), mode1))*(np.dot(
                    response1.transpose(), response1)))
    return mac


def plotGraph():
    Kcomp = Kc*complex(1., lossf)
    for f in f_array:
        omega = 2.*np.pi*f
        complex_matrix = [[8.*Kcomp-(omega**2)*m1, -4.*Kcomp],
                          [-4.*Kcomp, 4.*Kcomp-(omega**2)*m2]]
        d = la.det(complex_matrix)
        x1 = (4.*Kcomp*m_balance*eccen*omega**2)/d
        x2 = ((8.*Kcomp-(omega**2)*m1)*m_balance*eccen*omega**2)/d
        X1.append(abs(x1))
        X2.append(abs(x2))
        phase1.append(cm.phase(x1))
        phase2.append(cm.phase(x2))

    plt.plot(f_array, phase1, label="Lower")    
    plt.plot(f_array, phase2, label="Upper")
    plt.title("Phase Response of Floor Vibrations")
    plt.xlabel("Frequency, Hz")
    plt.ylabel("Phase, rad")
    #plt.minorticks_on()
    plt.legend()
    plt.grid(which = "major")
    #plt.grid(which = "minor")
    plt.savefig("phase response.png")
    plt.show()
    
    plt.plot(f_array, X1, label = "Lower")
    plt.plot(f_array, X2, label = "Upper")
    plt.title("Magnitude Response of Floor Vibrations")
    plt.xlabel("Frequency, Hz")
    plt.ylabel("X, m")
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend()
    plt.minorticks_on()    
    plt.grid(which = "major")
    plt.grid(which = "minor")
    plt.savefig("magnitude response.png")
    plt.show()