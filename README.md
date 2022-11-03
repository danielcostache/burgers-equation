<h1>
    Numerical solver for the Burgers' equation
</h1>

The following repository contains two scripts, written in Python and Julia, with the purpose of computing the numerical solution of the Burgers' Equation in fluid dynamics. 

<h2>
    Theoretical considerations
</h2>

Given a field $u(x,t)$ where $u$ represents the velocity of the fluid with the viscosity $\nu$ at a certain $x$ spatial coordinate and $t$ temporal coordinate, the Burger's equation is as following:
$$ \frac{\partial u}{\partial t} + u\frac{\partial u}{\partial x} = \nu\frac{\partial^2 u}{\partial x^2} $$
In order to find the numerical solution of the equation, the explicit Euler method is used. The value of the function at each step in time is computed using the following formula:
$$ u(x, t+\Delta t)=u(x,t)+\Delta t\left(-u(x,t)\frac{u(x+\Delta x, t)-u(x,t)}{\Delta x}+\nu\frac{u(x+\Delta x,t)-2u(x,t)+u(x-\Delta x,t)}{\Delta x^2}\right) $$
The spatial boundaries are considered absorbing.

<h2>
    Running the code
</h2>

The values for the parameters used in each test are written in the comments at the beginning of the code.

<h3>
    Python
</h3>

The tests have been run on [Python 3.9.7](https://www.python.org/downloads/release/python-397/), requiring the following modules: [matplotlib](https://matplotlib.org/) and [numpy](https://numpy.org/)

<h3>
    Julia
</h3>

The tests have been run on [Julia 1.8.2](https://julialang.org/downloads/), requiring the following package: [Plots](https://docs.juliaplots.org/stable/)
