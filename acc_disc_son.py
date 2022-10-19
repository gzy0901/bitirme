from numpy import exp

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
v = 1e16 #frequency

T_s = (3*G*M*M_dot)/(8*pi*R_in*sigma)
#T_x = T_s*x**(-3/4)*(1-x(-1/2))**(1/4) this gives an error because of the complexity of the expression
T_x = T_s*pow(x, (-3/4))*pow((1-pow(x, (-1/2))), (1/4))

# Define function to integrate
def f(x):
    return (4*pi*h*cosi*v**3*R_in**2)/(c**2*D**2)*(x/exp(h*v/k*T_x)-1)

# Implementing Simpson's 1/3 
def simpson13(x0,xn,n):
    # calculating step size
    h = (xn - x0) / n
    
    # Finding sum 
    integration = f(x0) + f(xn)
    
    for i in range(1,n):
        k = x0 + i*h
        
        if i%2 == 0:
            integration = integration + 2 * f(k)
        else:
            integration = integration + 4 * f(k)
    
    # Finding final integration value
    integration = integration * h/3
    
    return integration
    
# Input section

lower_limit = R_in
upper_limit = R_out
sub_interval = int(input("Enter number of sub intervals: "))

# Call trapezoidal() method and get result
result = simpson13(lower_limit, upper_limit, sub_interval)
print("Integration result by Simpson's 1/3 method is: %0.6f" % (result) )
