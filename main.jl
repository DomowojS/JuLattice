############################
## Main file for JuLattice #
############################
using Revise, MeshGrid, GLMakie

# User Settings
gridlengthX  = 400;
gridlengthY  = 100;

fluiddensity   = 100;
inflow_velocity= 0.1;     #Lattice unit

simulationTime = 8000;

cylinder_radius  = 10;
cylinder_position = [gridlengthX/4, gridlengthY/2]

## Simulation Settings
Q   = 9;
τ   = 1/√3 

lattice_velocity_unit = [   [0, -1, 0, 1, 0, -1, -1, 1, 1],
                            [0, 0, 1, 0, -1, -1, 1, 1, -1]];

lattice_velx = reshape(lattice_velocity_unit[1,],1,1,Q)
lattice_vely = reshape(lattice_velocity_unit[2,],1,1,Q)

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

# Initialize Plot arrays
vorticity = zeros(gridlengthX, gridlengthY);
vorticity_obs = Observable(vorticity);

    ## Set up the figure and axis       
        # Set up the figure and axis with explicit sizing
        fig = Figure(size = (1000, 600))
        ax = Axis(fig[1, 1], aspect = DataAspect(), title = "LBM Simulation")

        # Display the vorticity field using a heatmap with dynamic color range
        hm = heatmap!(ax, 1:gridlengthX, 1:gridlengthY, vorticity_obs, 
                    colormap = :curl, 
                    nan_color = :black,
                    colorrange = (-0.3, 0.3))
        Colorbar(fig[1, 2], hm, label = "Vorticity")

        # Set axis limits explicitly
        xlims!(ax, 1, gridlengthX)
        ylims!(ax, 1, gridlengthY)

        # Create a text element for time step display
        step_text = Observable("Time step: 0")
        text_obj = text!(ax, 2, gridlengthY-8, text = step_text, 
                color = :black, fontsize = 14)

        # Make sure figure is displayed explicitly
        display(fig)

# Run Simulation Loop
for i in 1:simulationTime

    # Get Macroscopic values
    global densityGrid = sum(distributions, dims=3);
    velocityX .= (1 ./ densityGrid) .* sum(distributions.*lattice_velx, dims=3); 
    velocityY .= (1 ./ densityGrid) .* sum(distributions.*lattice_vely, dims=3); 

    ## Apply Collision
    # Compute equilibrium state
    dotprod_velocities .= (lattice_velx .* velocityX) .+ (lattice_vely .* velocityY);
    distributions_equilibrium .= weights .* densityGrid .*(1 .+ 3 .*dotprod_velocities .+ 9/2 .*dotprod_velocities.^2 .- 3/2 .*(velocityX.^2 .+ velocityY.^2));
    # Relax towards equilibrium
    distributions .+= -(1/τ) .* (distributions .- distributions_equilibrium);

    # Stream 
    for j in 1:Q
        distributions[:,:,j] = circshift(distributions[:,:,j], (lattice_velx[j], lattice_vely[j]))
    end

    ## Apply Boundary conditions
    #Inlet velocity bc (unknown: f_1, f_8, f_9)
    densityGrid[inlet, :] .= (sum(distributions[inlet, [1,3,5]], dims=2).+ 2 .*sum(distributions[inlet, [2,6,7]], dims=2)) ./ (1-inflow_velocity)
    distributions[inlet, 4] .= distributions[inlet, 2] .+ (2/3 .* densityGrid[inlet,:] .* inflow_velocity)
    distributions[inlet, 8] .= distributions[inlet, 6] .+ (1/6 .* densityGrid[inlet,:] .* inflow_velocity) .- (1/2 .* (distributions[inlet, 3] .- distributions[inlet, 5]))
    distributions[inlet, 9] .= distributions[inlet, 7] .+ (1/6 .* densityGrid[inlet,:] .* inflow_velocity) .+ (1/2 .* (distributions[inlet, 3] .- distributions[inlet, 5]))


    #Outlet zero gradient bc
    distributions[outlet, [4, 8, 9]] .= distributions[gridlengthX-1, :, [4, 8, 9]]

    #No Slip Walls
    distributions[walls, 1:Q] .= distributions[walls, [1,4,5,2,3,8,9,6,7]];

    # Apply object boundary condition
    distributions[cylinder, 1:Q] .= distributions[cylinder, [1,4,5,2,3,8,9,6,7]];

        # Plot of the field
        if ((i % 10 == 0)) || (i == simulationTime)
            # Set velocities inside the cylinder to zero
            velocityX[cylinder] .= 0.0
            velocityY[cylinder] .= 0.0

            # Compute vorticity
            fill!(vorticity, 0.0)
            dv_dx = circshift(velocityY, (-1, 0)) .- circshift(velocityY, (1, 0))
            du_dy = circshift(velocityX, (0, -1)) .- circshift(velocityX, (0, 1))
            vorticity .= dv_dx .- du_dy
            vorticity[inlet] .= 0.0
            vorticity[outlet] .= 0.0

            # Mask the cylinder region
            vorticity[cylinder] .= NaN

            # Print some statistics to monitor the simulation
            valid_values = filter(!isnan, vorticity)
            if !isempty(valid_values)
                println("Step $i - Min: $(minimum(valid_values)), Max: $(maximum(valid_values))")
            end

            # Force the figure to update
            # Update the observable
            vorticity_obs[] = copy(vorticity)
            # Update time step text
            step_text[] = "Time step: $i"
            # Force a draw
            GLMakie.display(fig)
            yield()
            sleep(0.05)
        end

end

println("Everything Set")