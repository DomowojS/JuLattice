module Simulation_SetUp
using MeshGrid
export Set_Simulation_Params, Create_Grid, Initialize_Distributions, Simulation_Params, Create_Object

#Types
abstract type Geometry end

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

    Base.@kwdef struct Cylinder <: Geometry
        Radius::Real
        Position::Vector{Float64}
        lattice_radius::Union{Real, Nothing} = nothing
        lattice_position::Union{Vector{Float64}, Nothing} = nothing
        mask::Union{BitMatrix, Nothing} = nothing
        lattice_Reynolds::Union{Float64, Nothing} = nothing
    end#Rectangle

    Base.@kwdef struct Rectangle <: Geometry
        Width::Real
        Height::Real
        Angle::Real
        Position::Vector{Float64}
        lattice_width::Union{Real, Nothing} = nothing
        lattice_height::Union{Real, Nothing} = nothing
        lattice_position::Union{Vector{Float64}, Nothing} = nothing
        mask::Union{BitMatrix, Nothing} = nothing
        lattice_Reynolds::Union{Float64, Nothing} = nothing
    end#Rectangle

    struct none <: Geometry
        mask::BitMatrix
        lattice_Reynolds::Float64
    end#Rectangle

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


    #Object definitions
    function Create_Object(config, simulation::Simulation_Params, fluid::Fluid, grid::Grid)
        if isdefined(config, :Object_Type)
            type = config.Object_Type
                if type == "Cylinder"
                    object = Cylinder(Radius=config.Object_Radius, Position=config.Object_Position)
                elseif type == "Rectangle"
                    object = Rectangle(Width=config.Object_Width, Height=config.Object_Height, Position=config.Object_Position)  
                else 
                    error("Unsupported geometry type: $type")
                end
            lattice_dimensions = Convert_to_Lattice_Units(simulation.delta_x, object)
        else
            lattice_dimensions = nothing
            info("No Object in flow specified")
        end

        mask = Set_Object_Mask(grid, lattice_dimensions)
        lattice_Reynolds = Compute_Reynolds_Number(grid, fluid, lattice_dimensions)
        
        if type == "Cylinder"
            return Cylinder(config.Object_Radius, config.Object_Position, lattice_dimensions.Radius, lattice_dimensions.Position, mask, lattice_Reynolds)
        elseif type == "Rectangle"
            return Rectangle(config.Object_Width, config.Object_Height, config.Object_Angle, lattice_dimensions.Width, lattice_dimensions.Height, lattice_dimensions.Position, lattice_position, mask, lattice_Reynolds)
        else
            return none(mask, lattice_Reynolds)
        end                
    end#Create_Object

    function Convert_to_Lattice_Units(delta_x::Real, object::Cylinder)
        lattice_radius  = object.Radius/delta_x;
        lattice_position = object.Position ./ delta_x;
        return (Type="Cylinder", Radius=lattice_radius, Position=lattice_position)
    end#Convert_to_Lattice_Units

    function Convert_to_Lattice_Units(delta_x::Real, object::Rectangle)
        lattice_width  = object.width/delta_x;
        lattice_height = object.height/delta_x;
        lattice_position = object.Position ./ delta_x;
        return (Type="Rectangle", Width=lattice_width, Height=lattice_height, Angle=object.Angle, Position=lattice_position)
    end#Convert_to_Lattice_Units
    
    function Set_Object_Mask(grid::Grid, lattice_dimensions::NamedTuple{(:Type, :Radius, :Position)})
        gridX = grid.gridX
        gridY = grid.gridY
        Radius     = lattice_dimensions.Radius
        Position_X = lattice_dimensions.Position[1]
        Position_Y = lattice_dimensions.Position[2]
        
        return (gridX.-Position_X).^2 + (gridY.-Position_Y).^2 .< Radius.^2;
    end#Set_Object_Mask

    function Set_Object_Mask(grid::Grid, lattice_dimensions::NamedTuple{(:Type, :Width, :Height, :Angle, :Position)})
        gridX = grid.gridX
        gridY = grid.gridY        
        Width     = lattice_dimensions.Width
        Height    = lattice_dimensions.Height
        Angle     = deg2rad(lattice_dimensions.Angle)
        Position_X = lattice_dimensions.Position[1]
        Position_Y = lattice_dimensions.Position[2]

        
        # Translate grid points to the rectangle center
        gridX_rel = gridX .- Position_X
        gridY_rel = gridY .- Position_Y

        # Rotate the grid points around the center (Position_X, Position_Y)
        rotatedX = gridX_rel .* cos(Angle) .+ gridY_rel .* sin(Angle)
        rotatedY = -gridX_rel .* sin(Angle) .+ gridY_rel .* cos(Angle)

        # Now check if the points are within the bounds of the rotated rectangle
        return abs.(rotatedX) .<= Width / 2 .&& abs.(rotatedY) .<= Height / 2

    end#Set_Object_Mask

    function Set_Object_Mask(grid::Grid, lattice_dimensions::Nothing)
        gridX = grid.gridX
        return falses(gridX)
    end#Set_Object_Mask

    function Compute_Reynolds_Number(grid::Grid, fluid::Fluid, lattice_dimensions::NamedTuple{(:Type, :Radius, :Position)})      
        lattice_inflow_velocity = fluid.lattice_inflow_velocity
        lattice_viscosity = fluid.lattice_viscosity
        Characteristic_Length = lattice_dimensions.Radius

        lattice_Reynolds = (lattice_inflow_velocity .* Characteristic_Length)/lattice_viscosity

        return lattice_Reynolds
    end#Compute_Reynolds_Number

    function Compute_Reynolds_Number(grid::Grid, fluid::Fluid, lattice_dimensions::NamedTuple{(:Type, :Width, :Height, :Angle, :Position)})
        gridY = grid.gridY

        lattice_inflow_velocity = fluid.lattice_inflow_velocity
        lattice_viscosity = fluid.lattice_viscosity
        Characteristic_Length = maximum(gridY[mask]) - minimum(gridY[mask])

        lattice_Reynolds = (lattice_inflow_velocity .* Characteristic_Length)/lattice_viscosity

        return lattice_Reynolds
    end#Compute_Reynolds_Number

    function Compute_Reynolds_Number(grid::Grid, fluid::Fluid, lattice_dimensions::Nothing)
        gridY = grid.gridY

        lattice_inflow_velocity = fluid.lattice_inflow_velocity
        lattice_viscosity = fluid.lattice_viscosity
        Characteristic_Length = maximum(gridY) - minimum(gridY)

        lattice_Reynolds = (lattice_inflow_velocity .* Characteristic_Length)/lattice_viscosity

        return lattice_Reynolds
    end#Compute_Reynolds_Number


end#Simulation_SetUp