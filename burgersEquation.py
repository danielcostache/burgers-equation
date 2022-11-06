""" 
    Numerical solver for the Burgers' equation in fluid mechanics,
    describing the kinematics of a viscous fluid.

    The tests have been run on Python 3.9.7, using the following 
    values for the parameters:
        xi = -5
        xf = 5
        Nx = 600
        ti = 0
        tf = 2.5
        Nt = 1000
        nu = 0.015
"""

from matplotlib import pyplot as plt        # Module used for 2D plotting
from mpl_toolkits.mplot3d import Axes3D     # Module used for 3D plotting
import numpy as np                          # Module used for array operations
from matplotlib import animation            # Module used for animation rendering

xi = -5         # Initial spatial coordinate
xf = 5          # Final spatial coordinate
Nx = 600        # Number of spatial steps
ti = 0          # Initial temporal coordinate 
tf = 2.5        # Final temporal coordinate 
Nt = 1000       # Number of temporal steps

x = np.zeros(Nx+1)              # 1D array holding the spatial coordinates
t = np.zeros(Nt+1)              # 1D array holding the temporal coordinates
t[0] = ti
u = np.zeros((Nx+1, Nt+1))      # 2D array holding the velocity of each point at a certain spatial and temporal coordinate

dx = (xf-xi) / Nx       # Spatial step
dx2 = dx ** 2           # Square of the spatial step
dt = (tf-ti) / Nt       # Temporal step

nu = 0.015      # Viscosity of the fluid

for i in range(Nx+1):
    x[i] = xi + (i-1)*dx            # Generating values for the spatial coordinates
    u[i, 0] = np.exp(-x[i]**2)      # Generating values for the velocity at each spatial coordinate at the initial time moment

# Computing the numerical solution of the Burgers' equation
for kt in range(Nt):
    # Initial spatial boundary condition
    u[0, kt+1] = u[0, kt] + dt * (-u[0, kt] * (u[1, kt] - u[0, kt]) / dx + nu * (u[1, kt] - 2*u[0, kt]) / dx2)
    for jx in range(1, Nx):
        u[jx, kt+1] = u[jx, kt] + dt * (-u[jx, kt] * (u[jx+1, kt] - u[jx, kt]) / dx + nu * (u[jx+1, kt] - 2*u[jx, kt] + u[jx-1, kt]) / dx2)
    # Final spatial boundary condition
    u[Nx, kt+1] = u[Nx, kt] + dt * (-u[Nx, kt] * (-u[Nx, kt]) / dx + nu * (- 2 * u[Nx, kt] + u[Nx-1, kt]) / dx2)
    # Incrementing the time values 
    t[kt+1] = t[kt] + dt

# Plotting a surface graph of u(x,t) as a function of the two parameters
fig1 = plt.figure()
ax1 = Axes3D(fig1, auto_add_to_figure=False)
fig1.add_axes(ax1)
ax1.contour3D(t, x, u, 500)
ax1.set_title("Dependency of the speed on time and spatial coordinate")
ax1.set_xlabel("t")
ax1.set_ylabel("x")
ax1.set_zlabel("u(x,t)")
ax1.view_init(25, 0)
plt.savefig("images/python-tests/figure1.png", format="png")

# Creating an animation which pans the camera around the graph
figAni = plt.figure()
axAni = Axes3D(figAni, auto_add_to_figure=False)
figAni.add_axes(axAni)
axAni.set_title("Dependency of the speed on time and spatial coordinate")
axAni.set_xlabel("t")
axAni.set_ylabel("x")
axAni.set_zlabel("u(x,t)")

def init():
    axAni.contour3D(t, x, u, 500)
    return ()

def animate(i):
    axAni.view_init(elev=25, azim=i)
    return ()

anim = animation.FuncAnimation(figAni, animate, init_func=init, frames=360, interval=20, blit=True)
anim.save('images/python-tests/animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

# # Plotting 2D graphs of u(x,t) as a function of x at the initial and final time moment
fig2 = plt.figure()
plt.plot(x, u[:, 1], lw=3, label="Initial time moment")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.title("Dependency of the speed on the spatial coordinate")
plt.legend()
plt.savefig("images/python-tests/figure2.png", format="png")

fig3 = plt.figure()
plt.plot(x, u[:, Nt], lw=3, label="Final time moment")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.title("Dependency of the speed on the spatial coordinate")
plt.legend()
plt.savefig("images/python-tests/figure3.png", format="png")

fig4 = plt.figure()
plt.plot(x, u[:, 1], lw=3, label="Initial time moment")
plt.plot(x, u[:, Nt], lw=3, label="Final time moment")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.title("Dependency of the speed on the spatial coordinate")
plt.legend()
plt.savefig("images/python-tests/figure4.png", format="png")
