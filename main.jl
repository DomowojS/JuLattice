############################
## Main file for JuLattice #
############################
using Revise


# User Settings
gridlengthX  = 150;
gridlengthY  = 50;

fluiddensity      = 100;

simulationTime = 1000;

# Simulation Settings
Q   = 9;
τ   = 1/√3 

lattice_velocity_unit = [   [0, -1, 0, 1, 0, -1, -1, 1, 1],
                            [0, 0, 1, 0, -1, -1, 1, 1, -1]
                        ];

weights =   [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

# Initialize distributions arrays
distributions = ones(gridlengthX, gridlengthY, Q) .+ 0.01*rand(gridlengthX, gridlengthY, Q);
distributions_equilibrium = ones(gridlengthX, gridlengthY, Q);
distributions_streamed = ones(gridlengthX, gridlengthY, Q);

# Initialize macroscopic density and scale distribution
densityGrid = sum(distributions, dims=3);
distributions[:,:,1:Q] .*= fluiddensity ./ densityGrid[:,:,1];

# Initialize macroscopic velocity arrays
velocityX   = zeros(gridlengthX, gridlengthY, Q);
velocityY   = zeros(gridlengthX, gridlengthY, Q);

# Initialize macroscopic velocity
velocityX[:,:,1:Q] .= 1 ./ densityGrid[:,:,1] .* (sum(distributions[:,:,1:Q].*reshape(lattice_velocity_unit[1,],1,1,Q), dims=3)); 
velocityY[:,:,1:Q] .= 1 ./ densityGrid[:,:,1] .* (sum(distributions[:,:,1:Q].*reshape(lattice_velocity_unit[2,],1,1,Q), dims=3)); 

# Initialize macroscopic entities arrays

println("Everything set")