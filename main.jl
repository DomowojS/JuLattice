############################
## Main file for JuLattice #
############################
using Revise, MeshGrid

# User Settings
gridlengthX  = 150;
gridlengthY  = 50;

fluiddensity   = 100;

simulationTime = 1000;

cylinder_radius  = 10;
cylinder_position = [gridlengthX/3, gridlengthY/2]

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

# Run Simulation Loop
for i in 1:simulationTime

    # Get Macroscopic values
    global densityGrid = sum(distributions, dims=3);
    velocityX .= (1 ./ densityGrid) .* sum(distributions.*lattice_velx, dims=3); 
    velocityY .= (1 ./ densityGrid) .* sum(distributions.*lattice_vely, dims=3); 

    ## Apply Collision
    # Compute equilibrium state
    dotprod_velocities .= (lattice_velx .* velocityX) .+ (lattice_vely .* velocityY);
    distributions_equilibrium .= weights .* densityGrid .*(1 .+ 3 .*dotprod_velocities .+ 9/2 .*dotprod_velocities.^2 .- 3/2 .*(velocityX.^2 + velocityY.^2));
    # Relax towards equilibrium
    distributions .+= -(1/τ) .* (distributions .- distributions_equilibrium);

    # Stream 
    for j in 1:Q
        distributions[:,:,j] = circshift(distributions[:,:,j], (lattice_velx[j], lattice_vely[j]))
    end

    # Apply Boundary conditions

    # Apply object boundary condition
    distributions[cylinder, 1:Q] .= distributions[cylinder, [1,4,5,2,3,8,9,6,7]];

end

println("Everything Set")