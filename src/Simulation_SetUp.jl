module Simulation_SetUp
using MeshGrid
export Set_Simulation_Params, Create_Grid, Initialize_Distributions, Simulation_Params

#Structs
    struct Simulation_Params
        delta_x::Real;
        delta_t::Real;
        time_steps::Int64;
        gridlengthX::Int64;
        gridlengthY::Int64;
        Q::Int64;
        velocity_vector_x::Array{Int64, 3};
        velocity_vector_y::Array{Int64, 3};
        weights::Array{Float64, 3};
    end#simulation_params

    struct Fluid
        Kinematic_Viscosity::Real
        Fluid_Density::Real
        Inflow_Velocity::Real
        fluiddensity::Int64;    
        lattice_inflow_velocity::Float64;
        lattice_viscosity::Float64;
    end#fluid

    struct Grid
        gridX::Matrix{Int64}
        gridY::Matrix{Int64}
        Left_Boundary::BitMatrix
        Right_Boundary::BitMatrix
        Top_Boundary::BitMatrix
        Bottom_Boundary::BitMatrix
    end#Grid

    struct Mutable_Grid
        distributions::Array{Float64, 3}
        distributions_equilibrium::Array{Float64, 3}
        densityGrid::Array{Float64, 3}
        velocityX::Matrix{Float64}
        velocityY::Matrix{Float64}
        dotprod_velocity::Array{Float64, 3}
    end#Mutable_Grid

#Functions
    function Set_Simulation_Params(config)
        #Fixed values and read from config
        delta_x = config.delta_x 
        τ       = config.τ
        lattice_speedOfSound = 1 / √3;
        fluiddensity         = 100;

        delta_t = Compute_Time_Step(config.delta_x, τ, config.Kinematic_Viscosity, lattice_speedOfSound)
        lattice_inflow_velocity, lattice_viscosity, time_steps, gridlengthX, gridlengthY = 
            Convert_to_Lattice_Units(delta_t, delta_x, config.Inflow_Velocity, lattice_speedOfSound, τ, 
                                        config.Simulation_Time, config.length_X, config.length_Y)

        # Hardcoded simulation parameters for D2Q9 LBM
        Q   = 9;
        velocity_vector = [     [0, -1, 0, 1, 0, -1, -1, 1, 1],
                                [0, 0, 1, 0, -1, -1, 1, 1, -1]];
        velocity_vector_x = reshape(velocity_vector[1,],1,1,Q)
        velocity_vector_y = reshape(velocity_vector[2,],1,1,Q)
        weights =   [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
        weights =   reshape(weights,1,1,Q);


        simulation = Simulation_Params(delta_x, delta_t, time_steps, gridlengthX, gridlengthY, Q, velocity_vector_x, velocity_vector_y, weights);
        fluid      = Fluid(config.Kinematic_Viscosity, config.Fluid_Density, config.Inflow_Velocity, fluiddensity, lattice_inflow_velocity, lattice_viscosity);
        
        return simulation, fluid
    end#Set_Simulation_Params

    function Compute_Time_Step(delta_x, τ, Kinematic_Viscosity, lattice_speedOfSound)
        delta_t = ((τ - 0.5) * lattice_speedOfSound^2 * delta_x^2) / Kinematic_Viscosity
        return delta_t
    end#Compute_TimeStep

    function Convert_to_Lattice_Units(delta_t, delta_x, Inflow_Velocity, lattice_speedOfSound, τ, simulation_time, length_X, length_Y)
        lattice_inflow_velocity = Inflow_Velocity * (delta_t / delta_x);
        lattice_viscosity   = lattice_speedOfSound^2 * (τ -0.5);
        time_steps     = ceil(Int, simulation_time / delta_t);
        gridlengthX  = ceil(Int, length_X / delta_x);
        gridlengthY  = ceil(Int, length_Y / delta_x);
        return lattice_inflow_velocity, lattice_viscosity, time_steps, gridlengthX, gridlengthY
    end#Convert_to_Lattice_Units

    function Create_Grid(simulation::Simulation_Params)
        gridlengthX = simulation.gridlengthX
        gridlengthY = simulation.gridlengthY

        ##Identifiers
        # create grid
        gridX, gridY = meshgrid(1:gridlengthX, 1:gridlengthY);
        gridX, gridY = gridX', gridY';

        # create boundary indetifiers
        Left_Boundary    = gridX .== 1;
        Right_Boundary   = gridX .== gridlengthX;
        Top_Boundary     = gridY .== gridlengthY;
        Bottom_Boundary  = gridY .== 1;

        #Immutable grid - only identifiers
        grid = Grid(gridX, gridY, Left_Boundary, Right_Boundary, Top_Boundary, Bottom_Boundary)

        return grid
    end#Create_Grid

    function Initialize_Distributions(simulation::Simulation_Params, fluid::Fluid, grid::Grid)
        gridlengthX = simulation.gridlengthX
        gridlengthY = simulation.gridlengthY
        Q = simulation.Q
        fluiddensity = fluid.fluiddensity
        gridX = grid.gridX
        
        # Initialize distributions arrays
        distributions = ones(gridlengthX, gridlengthY, Q) .+ 0.01*rand(gridlengthX, gridlengthY, Q);
        distributions[:,:,4] .+= 2 .* (1 .+ 0.2 .* cos.(2 .* π .*gridX ./ gridlengthX .*4));
        distributions_equilibrium = ones(gridlengthX, gridlengthY, Q);

        # Initialize density and scale distribution
        densityGrid = sum(distributions, dims=3);
        distributions .*= fluiddensity ./ densityGrid;

        # Initialize macroscopic velocity arrays
        velocityX   = zeros(gridlengthX, gridlengthY);
        velocityY   = zeros(gridlengthX, gridlengthY);

        # Initialise dotproduct array 
        dotprod_velocities = zeros(gridlengthX, gridlengthY, Q);

        mutable_grid = Mutable_Grid(distributions, distributions_equilibrium, densityGrid, velocityX, velocityY, dotprod_velocities)

        return mutable_grid

    end#Initialize_Distributions

    function Create_Object(config, simulation::Simulation_Params, fluid::Fluid)
        return object
    end

        struct Cylinder_Object
            lattice_Re
        end#cylinder_object

        struct Rectangle_Object
        end#cylinder_object
        
    end#Create_Object



end#Simulation_SetUp