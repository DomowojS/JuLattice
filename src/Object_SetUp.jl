module Object_SetUp
using ..ConfigReader, ..Simulation_SetUp
export Create_Object, Geometry, Cylinder, Rectangle, none

    #Types
    abstract type Geometry end

    #Structs
    Base.@kwdef struct Cylinder <: Geometry
        Type::String
        Radius::Real
        Position::Vector{Float64}
        lattice_radius::Union{Real, Nothing} = nothing
        lattice_position::Union{Vector{Float64}, Nothing} = nothing
        mask::Union{BitMatrix, Nothing} = nothing
        lattice_Reynolds::Union{Float64, Nothing} = nothing
    end#Cylinder

    Base.@kwdef struct Rectangle <: Geometry
        Type::String
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
        Type::String
        mask::BitMatrix
        lattice_Reynolds::Float64
    end#Rectangle


    #Object definitions
    function Create_Object(config::Config, simulation::Simulation_Params, fluid::Fluid, grid::Grid)
        if isdefined(config, :Object_Type)
            type = config.Object_Type
                if type == "Cylinder"
                    object = Cylinder(Type="Cylinder", Radius=config.Object_Radius, Position=config.Object_Position)
                elseif type == "Rectangle"
                    object = Rectangle(Type="Rectangle", Width=config.Object_Width, Height=config.Object_Height, Position=config.Object_Position)  
                else 
                    error("Unsupported geometry type: $type")
                end
            lattice_dimensions = Convert_to_Lattice_Units(simulation.delta_x, object)
        else
            lattice_dimensions = nothing
            @info "Nothing in flow skipping object set-up"
        end

        mask = Set_Object_Mask(grid, lattice_dimensions)
        lattice_Reynolds = Compute_Reynolds_Number(grid, fluid, lattice_dimensions)
        
        if type == "Cylinder"
            return Cylinder("Cylinder", config.Object_Radius, config.Object_Position, lattice_dimensions.Radius, lattice_dimensions.Position, mask, lattice_Reynolds)
        elseif type == "Rectangle"
            return Rectangle("Rectangle", config.Object_Width, config.Object_Height, config.Object_Angle, lattice_dimensions.Width, lattice_dimensions.Height, lattice_dimensions.Position, lattice_position, mask, lattice_Reynolds)
        else
            return none("none", mask, lattice_Reynolds)
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
end#Object_SetUp