############################
## Main file for JuLattice #
############################
include("src/Plotter.jl")
include("src/Logger.jl")

using Revise, MeshGrid, GLMakie
using .Plotter, .Logger, .SimulationCore

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
Q   = 9;

velocity_vector = [     [0, -1, 0, 1, 0, -1, -1, 1, 1],
                        [0, 0, 1, 0, -1, -1, 1, 1, -1]];

velocity_vector_x = reshape(velocity_vector[1,],1,1,Q)
velocity_vector_y = reshape(velocity_vector[2,],1,1,Q)

weights =   [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
weights =   reshape(weights,1,1,Q);

# create grid
gridX, gridY = meshgrid(1:gridlengthX, 1:gridlengthY);
gridX, gridY = gridX', gridY';

# create object indetifier
cylinder = (gridX.-cylinder_position[1]).^2 + (gridY.-cylinder_position[2]).^2 .< cylinder_radius.^2;

# create boundary indetifiers
walls = gridY .== 1 .|| gridY .== gridlengthY;
inlet = gridX .== 1;
outlet = gridX .== gridlengthX;

# Initialize distributions arrays
distributions = ones(gridlengthX, gridlengthY, Q) .+ 0.01*rand(gridlengthX, gridlengthY, Q);
distributions[:,:,4] .+= 2 .* (1 .+ 0.2 .* cos.(2 .* π .*gridX ./ gridlengthX .*4));
distributions_equilibrium = ones(gridlengthX, gridlengthY, Q);

# Initialize macroscopic density and scale distribution
densityGrid = sum(distributions, dims=3);
distributions .*= fluiddensity ./ densityGrid;

# Initialize macroscopic velocity arrays
velocityX   = zeros(gridlengthX, gridlengthY);
velocityY   = zeros(gridlengthX, gridlengthY);

# Initialise dotproduct array 
dotprod_velocities = zeros(gridlengthX, gridlengthY, Q);

if any((Plotvorticity, Plotvx, Plotvy))
    if Plotvorticity==true 
        vorticity, vorticity_obs, text_obj, step_text, fig_vorticity = Create_Plot(gridlengthX, gridlengthY)
        screen1 = GLMakie.Screen()
        GLMakie.display(screen1, fig_vorticity)
    end
    if Plotvx==true 
        velocityX_obs, text_obj_vx, step_text_vx, fig_vx = Create_Plot(gridlengthX, gridlengthY, velocityX, "X")
        screen2 = GLMakie.Screen(; position = (600, 0))
        GLMakie.display(screen2, fig_vx)
    end
    if Plotvy==true 
        velocityY_obs, text_obj_vy, step_text_vy, fig_vy = Create_Plot(gridlengthX, gridlengthY, velocityY, "Y")
        screen3 = GLMakie.Screen()
        display(screen3, fig_vy)
    end
end

println("#################################")
println("Starting Simulation:")
# Run Simulation Loop
for i in 1:simulationTime

    # Get Macroscopic values
    global densityGrid = sum(distributions, dims=3);
    velocityX .= (1 ./ densityGrid) .* sum(distributions.*velocity_vector_x, dims=3); 
    velocityY .= (1 ./ densityGrid) .* sum(distributions.*velocity_vector_y, dims=3); 

    ## Apply Collision
    # Compute equilibrium state
    dotprod_velocities .= (velocity_vector_x .* velocityX) .+ (velocity_vector_y .* velocityY);
    distributions_equilibrium .= weights .* densityGrid .*(1 .+ 3 .*dotprod_velocities .+ 4.5 .*dotprod_velocities.^2 .- 1.5 .*(velocityX.^2 .+ velocityY.^2));
    # Relax towards equilibrium
    distributions .+= -(1/τ) .* (distributions .- distributions_equilibrium);

    # Stream 
    for j in 1:Q
        distributions[:,:,j] = circshift(distributions[:,:,j], (velocity_vector_x[j], velocity_vector_y[j]))
    end

    ## Apply Boundary conditions
    #Inlet velocity bc (unknown: f_1, f_8, f_9)
    densityGrid[inlet, :] .= (sum(distributions[inlet, [1,3,5]], dims=2).+ 2 .*sum(distributions[inlet, [2,6,7]], dims=2)) ./ (1-lattice_inflow_velocity)
    distributions[inlet, 4] .= distributions[inlet, 2] .+ (2/3 .* densityGrid[inlet,:] .* lattice_inflow_velocity)
    distributions[inlet, 8] .= distributions[inlet, 6] .+ (1/6 .* densityGrid[inlet,:] .* lattice_inflow_velocity) .- (1/2 .* (distributions[inlet, 3] .- distributions[inlet, 5]))
    distributions[inlet, 9] .= distributions[inlet, 7] .+ (1/6 .* densityGrid[inlet,:] .* lattice_inflow_velocity) .+ (1/2 .* (distributions[inlet, 3] .- distributions[inlet, 5]))


    #Outlet zero gradient bc
    distributions[outlet, [4, 8, 9]] .= distributions[gridlengthX-1, :, [4, 8, 9]]

    #No Slip Walls
    distributions[walls, 1:Q] .= distributions[walls, [1,4,5,2,3,8,9,6,7]];

    # Apply object boundary condition
    distributions[cylinder, 1:Q] .= distributions[cylinder, [1,4,5,2,3,8,9,6,7]];

        # Plot of the field
        if ((i % 10 == 0)) || (i == simulationTime)
            # Set velocities inside the cylinder to zero
            velocityX[cylinder] .= NaN
            velocityY[cylinder] .= NaN

            # Compute vorticity
            fill!(vorticity, 0.0)
            dv_dx = circshift(velocityY, (-1, 0)) .- circshift(velocityY, (1, 0))
            du_dy = circshift(velocityX, (0, -1)) .- circshift(velocityX, (0, 1))
            vorticity .= dv_dx .- du_dy
            vorticity[inlet] .= 0.0
            vorticity[outlet] .= 0.0

            # Mask the cylinder region
            vorticity[cylinder] .= NaN

            if ((i % 100 == 0)) || (i == simulationTime)
                Simulation_RuntimeLog(i, simulationTime)
            end
            # Update the observables
            if Plotvorticity==true 
                vorticity_obs[] = copy(vorticity) 
                step_text[] = "Time step: $i, $(floor(Int, i*delta_t))s"
            end
            if Plotvx==true 
                velocityX_obs[] = copy(velocityX) 
                step_text_vx[] = "Time step: $i, $(floor(Int, i*delta_t))s"
            end
            if Plotvy==true 
                velocityY_obs[] = copy(velocityY) 
                step_text_vy[] = "Time step: $i, $(floor(Int, i*delta_t))s"
            end

            yield()
            sleep(0.05)
        end

end

Log_Simulation_Tail()