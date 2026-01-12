############################
## Main file for JuLattice #
############################
include("src/Plotter.jl")
include("src/Logger.jl")

using MeshGrid, GLMakie
using .Plotter, .Logger

####################################  Initialize  ####################################
## User Settings
# Domain Settings
length_X = 4;              # m
length_Y = 1;              # m 

# Cylinder Definition
Radius   = 0.1    # m
Position = [1, 0.5] # m

# Fluid Settings 
Fluid_Density = 1000.0;       # kg/m^3
Inflow_Velocity = 0.4;      # m/s
Kinematic_Viscosity = 0.001; # m^2/s 

# Simulation Settings
Simulation_Time = 8000;     # s
delta_x = 0.01;             # discretisation in time and space
τ = 0.65;
Re = (Inflow_Velocity .* Radius)/Kinematic_Viscosity;
Re_Log=floor(Int,Re)

# Plot Requests
Plotvx = true;
Plotvy = true;
Plotvorticity = true;


#### Run Simulation #####
Log_Simulation_Header()
## Compute timestep from relaxation time
# Time step from relaxation time
lattice_speedOfSound = 1 / √3;
delta_t = ((τ - 0.5) * lattice_speedOfSound^2 * delta_x^2) / Kinematic_Viscosity

## Convert user settings to lattice units
# Domain
gridlengthX  = ceil(Int, length_X / delta_x);
gridlengthY  = ceil(Int, length_Y / delta_x);

# Cylinder
cylinder_radius  = Radius/delta_x;
cylinder_position = Position ./ delta_x;

# Fluid
fluiddensity = 100;
lattice_inflow_velocity = Inflow_Velocity * (delta_t / delta_x);
lattice_viscosity = lattice_speedOfSound^2 * (τ -0.5);
#ReynoldsCheck
lattice_Re = (lattice_inflow_velocity .* cylinder_radius)/lattice_viscosity;
lattice_Re_Log=floor(Int,lattice_Re)

#Log 
Log_Discretization_Settings(delta_x, delta_t, lattice_Re_Log)

# Simulation Settings
simulationTime = ceil(Int, Simulation_Time / delta_t);
# Q   = 9;

#Define arrays for each direction D2Q9
f00 = zeros(gridlengthX, gridlengthY) #center
fm0 = zeros(gridlengthX, gridlengthY) #left
f0m = zeros(gridlengthX, gridlengthY) #down
fp0 = zeros(gridlengthX, gridlengthY) #right
f0p = zeros(gridlengthX, gridlengthY) #up
fmm = zeros(gridlengthX, gridlengthY) #left-down
fmp = zeros(gridlengthX, gridlengthY) #left-up
fpp = zeros(gridlengthX, gridlengthY) #right-up
fpm = zeros(gridlengthX, gridlengthY) #right-down

#Define array for each direction after stream+Collision (S)
f00S = zeros(gridlengthX, gridlengthY) #center
fm0S = zeros(gridlengthX, gridlengthY) #left
f0mS = zeros(gridlengthX, gridlengthY) #down
fp0S = zeros(gridlengthX, gridlengthY) #right
f0pS = zeros(gridlengthX, gridlengthY) #up
fmmS = zeros(gridlengthX, gridlengthY) #left-down
fmpS = zeros(gridlengthX, gridlengthY) #left-up
fppS = zeros(gridlengthX, gridlengthY) #right-up
fpmS = zeros(gridlengthX, gridlengthY) #right-down

#Initialise macroscopic variables
rho = ones(gridlengthX, gridlengthY) .* fluiddensity
u = zeros(gridlengthX,gridlengthY) #vx
v = zeros(gridlengthX, gridlengthY) #vy

#define omega
omega = 1.0 / τ

#grid dimensions
rows = gridlengthY
cols = gridlengthX

#Initialise distribution functions 
for x in 1:gridlengthX
    for y in 1:gridlengthY
        #same ux and uy for all cells so no arrays for ux, uy
        ux = lattice_inflow_velocity
        uy = 0.0
        u[x,y] = ux
        v[x,y] = uy

        rho_init = fluiddensity

        f00[x,y] = rho_init * (-2.0 + 3.0*ux*ux) * (-2.0 + 3.0*uy*uy) / 9.0

        fm0[x,y] = rho_init * (1.0 - 3.0*ux + 3.0*ux*ux) * (-2.0 + 3.0*uy*uy) / 18.0
        fp0[x,y] = rho_init * (1.0 + 3.0*ux + 3.0*ux*ux) * (-2.0 + 3.0*uy*uy) / 18.0
        f0m[x,y] = rho_init * (-2.0 + 3.0*ux*ux) * (1.0 - 3.0*uy + 3.0*uy*uy) / 18.0
        f0p[x,y] = rho_init * (-2.0 + 3.0*ux*ux) * (1.0 + 3.0*uy + 3.0*uy*uy) / 18.0

        fmm[x,y] = rho_init * (1.0 - 3.0*ux + 3.0*ux*ux) * (1.0 - 3.0*uy + 3.0*uy*uy) / 36.0
        fmp[x,y] = rho_init * (1.0 - 3.0*ux + 3.0*ux*ux) * (1.0 + 3.0*uy + 3.0*uy*uy) / 36.0
        fpm[x,y] = rho_init * (1.0 + 3.0*ux + 3.0*ux*ux) * (1.0 - 3.0*uy + 3.0*uy*uy) / 36.0
        fpp[x,y] = rho_init * (1.0 + 3.0*ux + 3.0*ux*ux) * (1.0 + 3.0*uy + 3.0*uy*uy) / 36.0

        # #C++ Formeln
        # f00[x][y]=(((-2 + 3*(u[x][y]*u[x][y]))*(-2 + 3*(v[x][y]*v[x][y]))*rho[x][y])/9.);
    
        # fm0[x][y]=(-0.05555555555555555*((-2 + 3*(v[x][y]*v[x][y]))*rho[x][y]*(1 + 3*(u[x][y]*u[x][y]) - 3*u[x][y])));
        # fp0[x][y]=(-0.05555555555555555*((-2 + 3*(v[x][y]*v[x][y]))*rho[x][y]*(1 + 3*(u[x][y]*u[x][y]) + 3*u[x][y])));
        # f0m[x][y]=(-0.05555555555555555*((-2 + 3*(u[x][y]*u[x][y]))*rho[x][y]*(1 + 3*(v[x][y]*v[x][y]) - 3*v[x][y])));
        # f0p[x][y]=(-0.05555555555555555*((-2 + 3*(u[x][y]*u[x][y]))*rho[x][y]*(1 + 3*(v[x][y]*v[x][y]) + 3*v[x][y])));
        
        # fmm[x][y]=((rho[x][y]*(1 - 3*u[x][y] + 3*(u[x][y]*u[x][y]))*(1 - 3*v[x][y] + 3*(v[x][y]*v[x][y])))/36.);
        # fmp[x][y]=((rho[x][y]*(1 + 3*(u[x][y]*u[x][y]) - 3*u[x][y])*(1 + 3*(v[x][y]*v[x][y]) + 3*v[x][y]))/36.);
        # fpm[x][y]=((rho[x][y]*(1 + 3*(u[x][y]*u[x][y]) + 3*u[x][y])*(1 + 3*(v[x][y]*v[x][y]) - 3*v[x][y]))/36.);
        # fpp[x][y]=((rho[x][y]*(1 + 3*(u[x][y]*u[x][y]) + 3*u[x][y])*(1 + 3*(v[x][y]*v[x][y]) + 3*v[x][y]))/36.);
    end
end

#initialise fS-Arrays for the first time
f00S .= f00

fm0S .= fm0
fp0S .= fp0
f0mS .= f0m
f0pS .= f0p

fmmS .= fmm
fmpS .= fmp
fpmS .= fpm
fppS .= fpp


# create grid
gridX, gridY = meshgrid(1:gridlengthX, 1:gridlengthY);
gridX, gridY = gridX', gridY';

# create object indetifier
cylinder = (gridX.-cylinder_position[1]).^2 + (gridY.-cylinder_position[2]).^2 .< cylinder_radius.^2;

# create boundary indetifiers
walls = gridY .== 1 .|| gridY .== gridlengthY;
inlet = gridX .== 1;
outlet = gridX .== gridlengthX;

# # Initialize distributions arrays
# distributions = ones(gridlengthX, gridlengthY, Q) .+ 0.01*rand(gridlengthX, gridlengthY, Q);
# distributions[:,:,4] .+= 2 .* (1 .+ 0.2 .* cos.(2 .* π .*gridX ./ gridlengthX .*4));
# distributions_equilibrium = ones(gridlengthX, gridlengthY, Q);

# # Initialize macroscopic density and scale distribution
# densityGrid = sum(distributions, dims=3);
# distributions .*= fluiddensity ./ densityGrid;

# # Initialize macroscopic velocity arrays
# velocityX   = zeros(gridlengthX, gridlengthY);
# velocityY   = zeros(gridlengthX, gridlengthY);

# # Initialise dotproduct array 
# dotprod_velocities = zeros(gridlengthX, gridlengthY, Q);

if any((Plotvorticity, Plotvx, Plotvy))
    if Plotvorticity==true 
        vorticity, vorticity_obs, text_obj, step_text, fig_vorticity = Create_Plot(gridlengthX, gridlengthY)
        screen1 = GLMakie.Screen()
        GLMakie.display(screen1, fig_vorticity)
    end
    if Plotvx==true 
        velocityX_obs, text_obj_vx, step_text_vx, fig_vx = Create_Plot(gridlengthX, gridlengthY, u, "X")
        screen2 = GLMakie.Screen(; position = (600, 0))
        GLMakie.display(screen2, fig_vx)
    end
    if Plotvy==true 
        velocityY_obs, text_obj_vy, step_text_vy, fig_vy = Create_Plot(gridlengthX, gridlengthY, v, "Y")
        screen3 = GLMakie.Screen()
        display(screen3, fig_vy)
    end
end
####################################  Initialize  ####################################


####################################  Sim-Loop ####################################
println("#################################")
println("Starting Simulation:")
# Run Simulation Loop
for i in 1:simulationTime

    ###### NEW STABILIZATION #####
    for x in 2:cols-1
        for y in 2:rows-1
            # Get Macroscopic values
            rho[x, y] = f00[x, y] + (((fmm[x, y] + fpp[x, y]) + (fmp[x, y] + fpm[x, y])) + 
                                    ((fm0[x, y] + fp0[x, y]) + (f0p[x, y] + f0m[x, y])))
            
            u[x, y] = (((-fmm[x, y] + fpp[x, y]) + (-fmp[x, y] + fpm[x, y])) + 
                    ((-fm0[x, y] + fp0[x, y]))) / rho[x, y]
            
            v[x, y] = (((-fmm[x, y] + fpp[x, y]) + (fmp[x, y] - fpm[x, y])) + 
                    ((f0p[x, y] - f0m[x, y]))) / rho[x, y]
            
            # Push and Collision scheme
            # Collision is computed and initialised as the corresponding cell after Push
            fmmS[x-1, y-1] = fmm[x, y] + omega * ((rho[x, y] * (1 - 3*u[x, y] + 3*(u[x, y]*u[x, y])) * 
                                                (1 - 3*v[x, y] + 3*(v[x, y]*v[x, y]))) / 36.0 - fmm[x, y])
            
            f0mS[x, y-1] = f0m[x, y] + omega * (-0.05555555555555555 * ((-2 + 3*(u[x, y]*u[x, y])) * 
                                                rho[x, y] * (1 + 3*(v[x, y]*v[x, y]) - 3*v[x, y])) - f0m[x, y])
            
            fpmS[x+1, y-1] = fpm[x, y] + omega * ((rho[x, y] * (1 + 3*(u[x, y]*u[x, y]) + 3*u[x, y]) * 
                                                (1 + 3*(v[x, y]*v[x, y]) - 3*v[x, y])) / 36.0 - fpm[x, y])
            
            fm0S[x-1, y] = fm0[x, y] + omega * (-0.05555555555555555 * ((-2 + 3*(v[x, y]*v[x, y])) * 
                                                rho[x, y] * (1 + 3*(u[x, y]*u[x, y]) - 3*u[x, y])) - fm0[x, y])
            
            f00S[x, y] = f00[x, y] + omega * (((-2 + 3*(u[x, y]*u[x, y])) * (-2 + 3*(v[x, y]*v[x, y])) * 
                                            rho[x, y]) / 9.0 - f00[x, y])
            
            fp0S[x+1, y] = fp0[x, y] + omega * (-0.05555555555555555 * ((-2 + 3*(v[x, y]*v[x, y])) * 
                                                rho[x, y] * (1 + 3*(u[x, y]*u[x, y]) + 3*u[x, y])) - fp0[x, y])
            
            fmpS[x-1, y+1] = fmp[x, y] + omega * ((rho[x, y] * (1 + 3*(u[x, y]*u[x, y]) - 3*u[x, y]) * 
                                                (1 + 3*(v[x, y]*v[x, y]) + 3*v[x, y])) / 36.0 - fmp[x, y])
            
            f0pS[x, y+1] = f0p[x, y] + omega * (-0.05555555555555555 * ((-2 + 3*(u[x, y]*u[x, y])) * 
                                                rho[x, y] * (1 + 3*(v[x, y]*v[x, y]) + 3*v[x, y])) - f0p[x, y])
            
            fppS[x+1, y+1] = fpp[x, y] + omega * ((rho[x, y] * (1 + 3*(u[x, y]*u[x, y]) + 3*u[x, y]) * 
                                                (1 + 3*(v[x, y]*v[x, y]) + 3*v[x, y])) / 36.0 - fpp[x, y])
        end
    end

    #Swap: copy new distributions to array
    f00 .= f00S
    fm0 .= fm0S
    f0m .= f0mS
    fp0 .= fp0S
    f0p .= f0pS
    fmm .= fmmS
    fmp .= fmpS
    fpp .= fppS
    fpm .= fpmS
    
    ##### Boundary Conditions #####
    #Bounceback walls
    for x in 1:cols
        #bottom wall (y=1)
        f0p[x, 1] = f0m[x, 1]   #top = bottom
        fpp[x, 1] = fmm[x, 1]   #right-top = left-bottom
        fmp[x, 1] = fpm[x, 1]   #left-top = right-bottom

        #top wall(y=rows)
        f0m[x, rows] = f0p[x, rows] #bottom = top
        fmm[x, rows] = fpp[x, rows] #left-bottom = right-top
        fpm[x, rows] = fmp[x, rows] #right-bottom = left-top
    end

    #Bounceback cylinder
    for x in 1:cols
        for y in 1:rows
            if cylinder[x,y]
                fp0[x,y], fm0[x,y] = fm0[x,y], fp0[x,y] #horizontal getauscht
                f0p[x,y], f0m[x,y] = f0m[x,y], f0p[x,y] #vertikal getauscht
                fpp[x,y], fmm[x,y] = fmm[x,y], fpp[x,y] #diagonal getauscht rechtsoben <-> linksunten
                fpm[x,y], fmp[x,y] = fmp[x,y], fpm[x,y] #diagonal getauscht rechtsunten <-> linksoben
            end
        end
    end
    ##### Boundary Conditions #####

    ###### NEW STABILIZATION #####




    # # Get Macroscopic values
    # global densityGrid = sum(distributions, dims=3);
    # velocityX .= (1 ./ densityGrid) .* sum(distributions.*velocity_vector_x, dims=3); 
    # velocityY .= (1 ./ densityGrid) .* sum(distributions.*velocity_vector_y, dims=3); 

    # ## Apply Collision
    # # Compute equilibrium state
    # dotprod_velocities .= (velocity_vector_x .* velocityX) .+ (velocity_vector_y .* velocityY);
    # distributions_equilibrium .= weights .* densityGrid .*(1 .+ 3 .*dotprod_velocities .+ 4.5 .*dotprod_velocities.^2 .- 1.5 .*(velocityX.^2 .+ velocityY.^2));
    # # Relax towards equilibrium
    # distributions .+= -(1/τ) .* (distributions .- distributions_equilibrium);

    # # Stream 
    # for j in 1:Q
    #     distributions[:,:,j] = circshift(distributions[:,:,j], (velocity_vector_x[j], velocity_vector_y[j]))
    # end

    # ## Apply Boundary conditions
    # #Inlet velocity bc (unknown: f_1, f_8, f_9)
    # densityGrid[inlet, :] .= (sum(distributions[inlet, [1,3,5]], dims=2).+ 2 .*sum(distributions[inlet, [2,6,7]], dims=2)) ./ (1-lattice_inflow_velocity)
    # distributions[inlet, 4] .= distributions[inlet, 2] .+ (2/3 .* densityGrid[inlet,:] .* lattice_inflow_velocity)
    # distributions[inlet, 8] .= distributions[inlet, 6] .+ (1/6 .* densityGrid[inlet,:] .* lattice_inflow_velocity) .- (1/2 .* (distributions[inlet, 3] .- distributions[inlet, 5]))
    # distributions[inlet, 9] .= distributions[inlet, 7] .+ (1/6 .* densityGrid[inlet,:] .* lattice_inflow_velocity) .+ (1/2 .* (distributions[inlet, 3] .- distributions[inlet, 5]))


    # #Outlet zero gradient bc
    # distributions[outlet, [4, 8, 9]] .= distributions[gridlengthX-1, :, [4, 8, 9]]

    # #No Slip Walls
    # distributions[walls, 1:Q] .= distributions[walls, [1,4,5,2,3,8,9,6,7]];

    # # Apply object boundary condition
    # distributions[cylinder, 1:Q] .= distributions[cylinder, [1,4,5,2,3,8,9,6,7]];

        # Plot of the field
        if ((i % 10 == 0)) || (i == simulationTime)
            # Set velocities inside the cylinder to zero
            u[cylinder] .= NaN
            v[cylinder] .= NaN

            # velocityX[cylinder] .= NaN
            # velocityY[cylinder] .= NaN

            # Compute vorticity
            fill!(vorticity, 0.0)
            dv_dx = circshift(v, (-1, 0)) .- circshift(v, (1, 0))
            du_dy = circshift(u, (0, -1)) .- circshift(u, (0, 1))
            vorticity .= dv_dx .- du_dy
            vorticity[inlet] .= 0.0
            vorticity[outlet] .= 0.0

            # Mask the cylinder region
            vorticity[cylinder] .= NaN

            if ((i % 100 == 0)) || (i == simulationTime)
                Log_Simulation_Runtime(i, simulationTime)
            end
            # Update the observables
            if Plotvorticity==true 
                vorticity_obs[] = copy(vorticity) 
                step_text[] = "Time step: $i, $(floor(Int, i*delta_t))s"
            end
            if Plotvx==true 
                velocityX_obs[] = copy(u) 
                step_text_vx[] = "Time step: $i, $(floor(Int, i*delta_t))s"
            end
            if Plotvy==true 
                velocityY_obs[] = copy(v) 
                step_text_vy[] = "Time step: $i, $(floor(Int, i*delta_t))s"
            end

            yield()
            sleep(0.05)
        end

end

Log_Simulation_Tail()