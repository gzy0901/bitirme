#created by Gulsah Zeynep Yigit 

from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
from numpy import exp
from math import log
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
x_out = R_out/R_in

#empty list
frequency = []
Flux = []

v = np.logspace(16, 20, num=pow(10,5), endpoint=True, base=10.0, dtype=None, axis=0)
v_list = list(v)
print(type(v_list))

for item in v_list:

    T_s = pow(((3*G*M*M_dot)/(8*pi*R_in*sigma)), (1/4))
    
    T_x = T_s*pow(x, (-3/4))*pow((1-pow(x, (-1/2))), (1/4))

    # Define function to integrate
    def f(x):
        return(mpf(4*pi*h*cosi*v**3*R_in**2)/(c**2*D**2)*(x/exp(h*v/k*T_x)-1))

    Flux_int = integrate.quad(f, 1, 1000)
    
    frequency.append(v)
    Flux.append(Flux_int)




#writing to text file
file = open("acc_disc_list.txt", "w")
for index in range(len(frequency)):
    file.write(str(frequency[index]) + " " + str(Flux[index]) + "\n")
file.close()

plt.plot(frequency, Flux) 
  
# naming the x axis 
plt.xlabel('frequency (hv/kT_out) (Hz)') 
# naming the y axis 
plt.ylabel('Flux_v') 

plt.show() 

#multiplying freq values with h/kT

frequency_hkt= [item * (h/k*T_x) for item in frequency]

#changing the freq values
#float_frequency = []
#for item in frequency_hkt:
    #float_frequency.append(float(item))

#taking log
frequency_log = [log(x) for x in frequency]

#changing the flux value

Flux_abs = [item * (-10**16) for item in Flux]

#taking log
Flux_log = [log(x) for x in Flux_abs]

print (frequency_log)
print (Flux_log)
#writing the log values to txt file

file = open("acc_disc_log_list.txt", "w")
for index in range(len(frequency_log)):
    file.write(str(frequency_log[index]) + " " + str(Flux_log[index]) + "\n")
file.close()

# plotting the points  
plt.plot(frequency_log, Flux_log) 

# naming the x axis 
plt.xlabel('frequency (hv/kT_out)') 
# naming the y axis 
plt.ylabel('Flux_v') 

plt.show() 
