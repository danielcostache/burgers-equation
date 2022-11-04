#= 
    Numerical solver for the Burgers' equation in fluid mechanics,
    describing the kinematics of a viscous fluid.

    The tests have been run on Julia 1.8.2, using the following 
    values for the parameters:
        xi = -5
        xf = 5
        Nx = 600
        ti = 0
        tf = 2.5
        Nt = 1000
        ν = 0.015
=#
using Plots         # Package used for graphical representations 

xi = -5           # Initial spatial coordinate
xf = 5            # Final spatial coordinate
Nx = 600          # Number of spatial steps
ti = 0            # Initial temporal coordinate
tf = 2.5          # Final temporal coordinate
Nt = 1000         # Number of temporal steps
frq = 50

x = zeros(Nx+1)             # 1D array holding the spatial coordinates
t = zeros(Nt+1)             # 1D array holding the temporal coordinates
t[1] = ti
u = zeros(Nx+1, Nt+1)       # 2D array holding the velocity of each point at a certain spatial and temporal coordinate
noise = zeros(Nx+1, Nt+1)

Δx = (xf-xi) / Nx         # Spatial step
Δx2 = Δx^2                # Square of the spatial step
Δt = (tf-ti) / Nt         # Temporal step

ν = 0.015           # Viscosity of the fluid

for i = 1:Nx + 1
    x[i] = xi + (i-1)*Δx        # Generating values for the spatial coordinates
    u[i, 1] = exp(-x[i]^2)      # Generating values for the velocity at each spatial coordinate at the initial time moment
    noise[i, :] .= rand()*sin(2*pi*frq*x[i])/25
end

# Computing the numerical solution of the Burgers' equation
for kt = 1:Nt
    # Initial spatial boundary condition
    u[1, kt+1] = u[1, kt] + Δt * (-u[1, kt] * (u[2, kt] - u[1, kt]) / Δx + ν * (u[2, kt] - 2*u[1, kt]) / Δx2)
    for jx = 2:Nx
        u[jx, kt+1] = u[jx, kt] + Δt * (-u[jx, kt] * (u[jx+1, kt] - u[jx, kt]) / Δx + ν * (u[jx+1, kt] - 2*u[jx, kt] + u[jx-1, kt]) / Δx2)
    end
    # Final spatial boundary condition
    u[Nx+1, kt+1] = u[Nx+1, kt] + Δt * (-u[Nx+1, kt] * (-u[Nx+1, kt]) / Δx + ν * (- 2*u[Nx+1, kt] + u[Nx, kt]) / Δx2)
    # Incrementing the time
    t[kt + 1] = t[kt] + Δt
end

# Plotting a surface graph of u(x,t) as a function of the two parameters
surface(t, x, (u+noise), title="Dependency of the speed on time and spatial coordinate", xaxis="t", yaxis="x", zaxis="u(x,t)", camera = (90, 25))
savefig("images/julia-noise3d-angle.png")

surface(t, x, (u+noise), title="Dependency of the speed on time and spatial coordinate", xaxis="t", yaxis="x", zaxis="u(x,t)", camera = (90, 5))
savefig("images/julia-noise3d-facing.png")

# Plotting 2D lines of u(x,t) as a function of x at the initial and final time moment
plot(x, u[:, 1], lw = 5, title="Dependency of the speed on the spatial coordinate", xaxis="x", yaxis="u(x)", label="Initial time moment")
plot!(x, (u[:, 1]+noise[:, 1]), lw = 2, title="Dependency of the speed on the spatial coordinate", xaxis="x", yaxis="u(x)", label="Initial time moment")
savefig("images/julia-noise2d-initial.png")

plot(x, u[:, Nt+1], lw = 5, title="Dependncy of the speed on the spatial coordinate", xaxis="x", yaxis="u(x)", label="Final time moment", legend=:topleft)
plot!(x, (u[:, Nt+1]+noise[:, Nt+1]), lw = 2, title="Dependncy of the speed on the spatial coordinate", xaxis="x", yaxis="u(x)", label="Final time moment", legend=:topleft)
savefig("images/julia-noise2d-final.png")
