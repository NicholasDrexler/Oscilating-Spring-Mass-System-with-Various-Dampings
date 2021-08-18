# Oscilating-Spring-Mass-System-with-Various-Dampings
For the non-physicists / non-engineers:      
We consider the sum of forces acting on a spring-mass system in the direction of the bob's axial motion.     
These are the friction force, spring force and newtons second law.     
The sum of these forces can be written:  b * (velocty) - K * (position) = (mass)*(acceleration)     

We can then express this as a Linear 2nd Order Differential Equation, where acceleration is the 2nd derivative of position and velocity is the 1st derivative of position.        The system is thus directly affect by b, which is the coeffient of friction (damping coefficient) for the spring, K which is the spring constant which is how stiff the spring is, and the mass of the bob.          

All differential equations require initial conditions for the position and for the velocity. The amount dt, is an iterator that refines detail.   

INSTRUCTIONS: 
Only adjust the quantities in the SETTINGS section of the code.   
Depending on those settings, the code will abutomatically determine if the system is underdamped, overdamped, critically damped, or undamped, and  produce the proper graph with dampening type in the title and the  differential equation in the label. 
