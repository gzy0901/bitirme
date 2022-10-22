#created by Gulsah Zeynep Yigit 
#kodu açmak için code
#çalıştırmak için python

from numpy import exp
from mpmath import *

mp.dps = 300

G = 6.6743e-8 # cm^3/g.s
M_sun = 1.989e+33 # g
M = 1.4*M_sun
M_dot = 1e+16
pi = 3.1415
h = 6.6261e-27 # cm^2.g/s
cosi = 1
sigma = 5.6704e-5 #sigma g/s^3.K^4
R_in = 1e+8 #cm
R_out = 1000*R_in
ly = 9.461e+18 #lightyear in cm
D = 9*ly #nearest magnetar
pc = 3.086e+18 #parsec in cm
c = 3e+10 #cm/s
k = 1.3807e-16 #cm^2.g/s^2.K
x = R_out/R_in
T_s = (3*G*M*M_dot)/(8*pi*R_in*sigma)
#T_x = T_s*x**(-3/4)*(1-x(-1/2))**(1/4) this gives an error because of the complexity of the expression
T_x = T_s*pow(x, (-3/4))*pow((1-pow(x, (-1/2))), (1/4))
lower_limit = R_in
upper_limit = R_out
sub_interval = 50
v = 1e16 #frequency

def main():
    global v
    
    #empty list
    frequency = []
    flux = []
    
    while (v < 1.5*1e16):
        
        v*=1.5
   
        frequency.append(v)
        
        flux.append(f(x))

    # Call trapezoidal() method and get result
    result = simpson13(lower_limit, upper_limit, sub_interval)

    print(f"Integration result by Simpson's 1/3 method is: {result}" )
    
    print(frequency)
    print(flux)
    
    # for i in frequency:
    #     print (i)
    # for i in flux:
    #     print (i)

        
# Define function to integrate
def f(x):
    return(mpf((4*pi*h*cosi*v**3*R_in**2)/(c**2*D**2)*(x/exp(h*v/k*T_x)-1)))

# Implementing Simpson's 1/3 
def simpson13(a,b,n):
    if n<1 or a>b:
        print ("Hatali veri!")
    elif n%2 == 1:
        print ("n cift degil!")
    else:
        h = (b-a)/n
        s = f(a) + f(b)
        for i in range (1,n):
            katsayi = 2*(i%2+1)
            x = a+i*h
            s = s+katsayi*f(x)
    return h*s/3.0

    
main()
