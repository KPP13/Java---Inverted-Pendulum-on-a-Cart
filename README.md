# Java---Inverted-Pendulum-on-a-Cart

The program is simple headless simulation of inverted pendulum on a cart system. 

Program creates single output file with six columns of data: t (simulation time in seconds), x1 (cart position [m]), x2 (cart velocity [m/s]), x3 (pendulum position [rad]), x4 (pedulum velocity [rad/s]), u (control force [N]). They contain calculated time series of state variables that might be plotted in the other program (for example Matlab).

Differential equations of inverted pendulum on a cart system was derived from Euler-Lagrange equations and Rayleigh's dissipation function. The full first order model of the system is shown below:
dx1 = x2
dx2 = (e1*x2 + e2*x4*cos(x3) + e3*x4^2*sin(x3) + e4*sin(2*x3) + e5*u)/(k1 + k2*cos(x3)^2)
dx3 = x3
dx4 = (f1*x2*cos(x3) + f2*x4 + f3*x4^2*sin(2*x3) + f4*sin(x3) + f5*u*cos(x3))/(k1 + k2*cos(x3)^2)

The values of k1, k2, e1,...,e5, f1,...,f5 coefficients may be found in the code.

The inverted pendulum system describes the real inverted pendulum plant. All physical coefficients were found during the identification process.

In the program two kind of regulators were implemented: swing up regulator and LQ regulator to stabilise the system.

The program is just a simulation of physical object. I am going to create a GUI with charts and animation of the plant soon. I will use Spring Framework probably (but I am still learning ;) ).

Output examples are presented in 'output - examples' directory.