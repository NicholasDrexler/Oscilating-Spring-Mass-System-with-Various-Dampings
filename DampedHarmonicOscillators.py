"""
Oscilating Spring-Mass System with Various Dampings
Written by: Nicholas Drexler 

For the non-physicist / non-engineers:

    We consider the sum of forces acting on a spring-mass system in the
    direction of the bob's axial motion.
    These are the friction force, spring force and newtons second law.
    The sum of these forces can be written:
        b * (velocty) - K * (position) = (mass)*(acceleration)
    We can then express this as a Linear 2nd Order Differential Equation,
    where acceleration is the 2nd derivative of position and velocity is the
    1st derivative of position. 
 
    The system is thus directly affect by b, which is the coeffient of 
    friction (damping coefficient) for the spring, K which is the 
    spring constant which is how stiff the spring is, and the mass of the bob. 
    
    All differential equations require initial conditions for the position
    and for the velocity. The amount dt, is an iterator that refines detail. 

INSTRUCTIONS
Only adjust the quantities in the SETTINGS section of the code. 

Depending on those settings, the code will abutomatically determine if the
system is underdamped, overdamped, critically damped, or undamped, and 
produce the proper graph with dampening type in the title and the 
differential equation in the label. 

"""

# VERSION 1

import numpy as np
import matplotlib.pyplot as plt

### SETTINGS #############################################
### ADJUST ONLY THESE TO PRODUCE A GRAPH #################
mass = 1.0          # Mass of the bob in Kg
b = 1.0             # Damping constant
K = 3.0             # Spring Constant
x_i_cond = 1        # need initial conditions (t=0)
v_i_cond = 0        # need initial conditions (t=0)
t = 0               # Initial time
dt = 0.1            # Time differential
duration = 10       # Duration of observation
##########################################################

# internal maths
w = np.sqrt(K/mass)                 
discriminant = 4 * mass * K
graph_type = str

# useful indicators
T = 2 * np.pi * np.sqrt(mass/K)     
f = 1/T
bsqr = b**2

# case-swap logic
nodamp = False
overdamp = False
underdamp = False
critdamp = False

def case_check():
    """
    Checks the b^2 relation to discriminant
    and evaluates which case we are in.
    """ 
    
    global nodamp
    global overdamp
    global underdamp
    global critdamp
    
    # check the case
    if bsqr > discriminant:
        nodamp = False
        overdamp = True
        underdamp = False
        critdamp = False
    elif bsqr < discriminant:
        nodamp = False
        overdamp = False
        underdamp = True
        critdamp = False
    elif bsqr == discriminant:
        nodamp = False
        overdamp = False
        underdamp = False
        critdamp = True

def overdamp_math():
    """
    Solves the Differential Equation, for the overdamped case
    considering only 2 possible coeffient solutions
    overdamp: b^2 > 4 * mass * K
    EQN: x = c1*e^(r1*t) + c2*e^(r2*t)
    """
    
    global r1
    global r2
    global c1
    global c2
    
    r1 = (-b + np.sqrt(b**2 - discriminant))/ (2*mass)
    r2 = (-b - np.sqrt(b**2 - discriminant)) / (2*mass)
    
    eq1 = np.array([1,1,x_i_cond])          # this must always be the case
    eq2 = np.array([r1,r2,v_i_cond])        # this will always be the case 
    
    # solving the system of equations
    if eq1[0] != 0 and eq2[0] != 0:
        eq1_s = eq1 * (-1*eq2[0])
        eq2_s = eq2 * eq1[0]
        eq_sum = eq1_s + eq2_s
    elif eq1[0] == 0 or eq2[0] == 0:
        eq_sum = eq1_s + eq2_s
            
    c2 = eq_sum[2] / (eq_sum[0] + eq_sum[1])
    c1 = x_i_cond - c2

def underdamp_math():
    """
    Solves the Differential Equation, for the underdamped case
    considering only 2 possible coeffient solutions
    overdamp: b^2 < 4 * mass * K
    EQN: x(t) = x = e^(-bt/2m) * ( c1*cos(wd*t) +  c2*sin(wd*t) )
    Occurs when the roots are real and different
    """
    
    global wd # damped angualr frequency
    wd = np.sqrt(np.abs(b**2 - discriminant)) / (2*mass)

    eq1 = np.array([np.cos(0),np.sin(0),x_i_cond])
    eq1_list = eq1.tolist()
    eq1_list.pop(-1)
    for i in eq1:
        if i == 0:
             continue
        elif i != 0:
             global c1
             c1 = i
    
    global c2
    c2 = (c1*b)/(2*mass*wd)
    
def critdamp_math():        
    """
    Solves the Differential Equation, for the critically damped case
    considering only 2 possible coeffient solutions
    overdamp: b^2 = 4 * mass * K
    EQN: x(t) = e^(-bt/2m) * (c1 + c2*t)
    """
    global c1
    global c2
    c1 = x_i_cond               # always true since t=0 for initial conditions
    c2 = v_i_cond + (c1 * b)/(2*mass)

    
# check to see if this is the nodamp version first
if b == 0:
    nodamp = True
else:
    nodamp = False
    case_check()

# we call the correct function based on values in settings
xvals = []
tvals = []
if nodamp == True:
    graph_type = "No Damping"
    while t < duration: 
        x = x_i_cond * np.cos(w*t)
        t += dt
        xvals.append(x)
        tvals.append(t) 
elif overdamp == True:
    graph_type = "Overdamping"
    overdamp_math()
    while t < duration: 
        x = c1*np.exp(r1*t) + c2*np.exp(r2*t)
        t += dt
        xvals.append(x)
        tvals.append(t)      
elif underdamp == True:
    graph_type = "Underdamping"
    underdamp_math()
    while t < duration: 
        x = np.exp(-(b*t)/(2*mass)) * ( (c1)*np.cos(wd*t) +  (c2)*np.sin(wd*t) )
        t += dt
        xvals.append(x)
        tvals.append(t)       
elif critdamp == True:
    graph_type = "Critical damping"
    critdamp_math()
    while t < duration: 
        x = np.exp(-(b*t)/(2*mass)) * (c1 + c2*t)
        t += dt
        xvals.append(x)
        tvals.append(t)          

# graphing 
plt.figure(figsize=(15,8))
plt.axhline(y=0, color='k', linestyle='-', linewidth=3) # Bold horizontal line at y=0

# handles sign + or - in the label
sign = lambda x: '+ ' + str(float(x)) if float(x) >= 0 else '- ' + str(abs(float(x)))

DE_label = f"Differental Equation: ${mass} \ddotx {sign(b)} \dotx {sign(K)}x$ \
    \nInitial conditions: x(0) = {x_i_cond} and $\dotx$ (0) = {v_i_cond}"

plt.plot(tvals,xvals,'.-r',label=(DE_label))

plt.title("Damped Harmonic Oscillation: Mass-Spring System \
            \n {} - Position v.s. Time".format(graph_type))
           
plt.xlabel("Time (s)")
plt.ylabel("Positoin (m)")
plt.xlim(tvals[0],tvals[-1])

plt.legend(loc='upper right', shadow=True, fontsize=16, markerscale=4)
plt.grid()
plt.show()