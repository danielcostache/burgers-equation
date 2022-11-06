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
frq = 50

x = np.zeros(Nx+1)              # 1D array holding the spatial coordinates
t = np.zeros(Nt+1)              # 1D array holding the temporal coordinates
t[0] = ti
u = np.zeros((Nx+1, Nt+1))      # 2D array holding the velocity of each point at a certain spatial and temporal coordinate
noise = np.zeros((Nx+1, Nt+1))

dx = (xf-xi) / Nx       # Spatial step
dx2 = dx ** 2           # Square of the spatial step
dt = (tf-ti) / Nt       # Temporal step

nu = 0.015      # Viscosity of the fluid

for i in range(Nx+1):
    x[i] = xi + (i-1)*dx            # Generating values for the spatial coordinates
    u[i, 0] = np.exp(-x[i]**2)      # Generating values for the velocity at each spatial coordinate at the initial time moment
    noise[i, :] = np.random.rand()*np.sin(2*np.pi*frq*x[i])/25


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
ax = plt.axes(projection="3d")
ax.contour3D(t, x, (u+noise), 500)
ax.set_title("Dependency of the speed on time and spatial coordinate")
ax.set_xlabel("t")
ax.set_ylabel("x")
ax.set_zlabel("u(x,t)")
ax.view_init(25, 0)
plt.savefig("images/python-noise3d-angle.png", format="png")

fig2 = plt.figure()
ax = plt.axes(projection="3d")
ax.contour3D(t, x, (u+noise), 500)
ax.set_title("Dependency of the speed on time and spatial coordinate")
ax.set_xlabel("t")
ax.set_ylabel("x")
ax.set_zlabel("u(x,t)")
ax.view_init(5, 0)
plt.savefig("images/python-noise3d-facing.png", format="png")

# Creating an animation which pans the camera around the graph
figAni1 = plt.figure()
figAni2 = plt.figure()
axAni1 = Axes3D(figAni1, auto_add_to_figure=False)
axAni2 = Axes3D(figAni2, auto_add_to_figure=False)
figAni1.add_axes(axAni1)
figAni2.add_axes(axAni2)
axAni1.set_title("Dependency of the speed on time and spatial coordinate")
axAni1.set_xlabel("t")
axAni1.set_ylabel("x")
axAni1.set_zlabel("u(x,t)")
axAni2.set_title("Dependency of the speed on time and spatial coordinate")
axAni2.set_xlabel("t")
axAni2.set_ylabel("x")
axAni2.set_zlabel("u(x,t)")

def init():
    axAni1.contour3D(t, x, (u+noise), 500)
    axAni2.contour3D(t, x, (u+noise), 500)
    return ()

def animate(i):
    axAni1.view_init(elev=25, azim=i)
    axAni2.view_init(elev=5, azim=i)
    return ()

anim1 = animation.FuncAnimation(figAni1, animate, init_func=init, frames=360, interval=20, blit=True)
anim1.save('images/python-noise3d-angle-animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
anim2 = animation.FuncAnimation(figAni2, animate, init_func=init, frames=360, interval=20, blit=True)
anim2.save('images/python-noise3d-facing-animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

# Plotting 2D graphs of u(x,t) as a function of x at the initial and final time moment
fig3 = plt.figure()
plt.plot(x, u[:, 0], lw=3, label="Initial time moment")
plt.plot(x, (u[:, 0]+noise[:, 0]), lw=1, label="Noisy signal")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.title("Dependency of the speed on the spatial coordinate")
plt.legend()
plt.savefig("images/python-noise2d-initial.png", format="png")

fig4 = plt.figure()
plt.plot(x, u[:, Nt], lw=3, label="Final time moment")
plt.plot(x, (u[:, Nt]+noise[:, Nt]), lw=1, label="Noisy signal")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.title("Dependency of the speed on the spatial coordinate")
plt.legend()
plt.savefig("images/python-noise2d-final.png", format="png")
