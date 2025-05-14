module Object_SetUp
using Images, FileIO
using ..ConfigReader, ..Simulation_SetUp
export Create_Object, Geometry, Cylinder, Rectangle, Custom, none

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

    Base.@kwdef struct Custom <: Geometry
        Type::String
        Width::Real
        Height::Real
        Angle::Real
        Position::Vector{Float64}
        Path_to_png::String
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
        if config.Object_Type in ["Cylinder", "Rectangle", "Custom"]
            type = config.Object_Type
                if type == "Cylinder"
                    object = Cylinder(Type="Cylinder", Radius=config.Object_Radius, Position=config.Object_Position)
                elseif type == "Rectangle"
                    object = Rectangle(Type="Rectangle", Width=config.Object_Width, Height=config.Object_Height, Angle=config.Object_Angle, Position=config.Object_Position)  
                elseif type == "Custom"
                    object = Custom(Type="Custom", Width=config.Object_Width, Height=config.Object_Height, Angle=config.Object_Angle, Position=config.Object_Position, Path_to_png=config.Object_path_to_png)  
                else 
                    error("Unsupported geometry type: $type")
                end
            lattice_dimensions = Convert_to_Lattice_Units(simulation.delta_x, object)
        else
            lattice_dimensions = nothing
        end

        mask = Set_Object_Mask(grid, lattice_dimensions)
        lattice_Reynolds = Compute_Reynolds_Number(grid, fluid, lattice_dimensions, mask)
        
        if config.Object_Type == "Cylinder"
            return Cylinder("Cylinder", config.Object_Radius, config.Object_Position, lattice_dimensions.Radius, lattice_dimensions.Position, mask, lattice_Reynolds)
        elseif config.Object_Type == "Rectangle"
            return Rectangle("Rectangle", config.Object_Width, config.Object_Height, config.Object_Angle, config.Object_Position, lattice_dimensions.Width, lattice_dimensions.Height, lattice_dimensions.Position, mask, lattice_Reynolds)
        elseif type == "Custom"
            return Custom("Custom", config.Object_Width, config.Object_Height, config.Object_Angle, config.Object_Position, config.Object_path_to_png , lattice_dimensions.Width, lattice_dimensions.Height, lattice_dimensions.Position, mask, lattice_Reynolds)  
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
        lattice_width  = object.Width/delta_x;
        lattice_height = object.Height/delta_x;
        lattice_position = object.Position ./ delta_x;
        return (Type="Rectangle", Width=lattice_width, Height=lattice_height, Angle=object.Angle, Position=lattice_position)
    end#Convert_to_Lattice_Units

    function Convert_to_Lattice_Units(delta_x::Real, object::Custom)
        lattice_width  = object.Width/delta_x;
        lattice_height = object.Height/delta_x;
        lattice_position = object.Position ./ delta_x;
        return (Type="Custom", Width=lattice_width, Height=lattice_height, Angle=object.Angle, Position=lattice_position, Path=object.Path_to_png)
    end#Convert_to_Lattice_Units

    function Set_Object_Mask(grid::Grid, lattice_dimensions::NamedTuple{(:Type, :Radius, :Position)})
        gridX = grid.gridX
        gridY = grid.gridY
        Radius     = lattice_dimensions.Radius
        Position_X = lattice_dimensions.Position[1]
        Position_Y = lattice_dimensions.Position[2]
        
        return (gridX.-Position_X).^2 + (gridY.-Position_Y).^2 .< Radius.^2;
    end#Set_Object_Mask

    function Set_Object_Mask(grid::Grid, lattice_dimensions::NamedTuple{(:Type, :Width, :Height, :Angle, :Position, :Path)})

        img = load(lattice_dimensions.Path)
        custom_object = BitMatrix(Gray.(img) .> 0.5)

        mask = place_custom_object_in_grid(grid, custom_object, lattice_dimensions)

        return mask

    end#Set_Object_Mask

    function Set_Object_Mask(grid::Grid, lattice_dimensions::Nothing)
        gridX = grid.gridX
        return falses(size(gridX))
    end#Set_Object_Mask

    function place_custom_object_in_grid(grid::Grid, custom_object::BitMatrix, lattice_dimensions::NamedTuple{(:Type, :Width, :Height, :Angle, :Position, :Path)})
        gridX = grid.gridX
        gridY = grid.gridY        
        Width     = lattice_dimensions.Width
        Height    = lattice_dimensions.Height
        Angle     = deg2rad(lattice_dimensions.Angle)
        Position_X = lattice_dimensions.Position[1]
        Position_Y = lattice_dimensions.Position[2]    
        
        mask = falses(size(gridX))
            
        # Convert custom_object to grayscale image for transformations
        custom_object_img = Gray.(custom_object)
        
        # Calculate grid cell size
        dx = abs(gridX[min(2,end),1] - gridX[1,1])
        dy = abs(gridY[1,min(2,end)] - gridY[1,1])
        
        # Scale custom_object to grid resolution
        cells_width = round(Int, Width / dx)
        cells_height = round(Int, Height / dy)
        scaled_custom_object = imresize(custom_object_img, (cells_height, cells_width))
        
        # Rotate image (centered), pad with white (1.0)
        rotated_custom_object = imrotate(scaled_custom_object, -Angle, axes(scaled_custom_object), fillvalue=1.0)
        rot_height, rot_width = size(rotated_custom_object)

        # Loop over grid
        for i in 1:size(gridX, 1)
            for j in 1:size(gridX, 2)
                x, y = gridX[i, j], gridY[i, j]

                # Translate to object-centered coordinates
                rel_x = x - Position_X
                rel_y = y - Position_Y

                # Inverse rotation (to map from grid space back to object space)
                inv_rot_x =  rel_x * cos(Angle) + rel_y * sin(Angle)
                inv_rot_y = -rel_x * sin(Angle) + rel_y * cos(Angle)

                # Flip y to account for image coordinates (top-down)
                pixel_x = round(Int, rot_width / 2 + inv_rot_x / dx)
                pixel_y = round(Int, rot_height / 2 - inv_rot_y / dy)  # minus here flips the y-axis

                # Check bounds
                if 1 <= pixel_x <= rot_width && 1 <= pixel_y <= rot_height
                    mask[i, j] = rotated_custom_object[pixel_y, pixel_x] < 0.5
                end
            end
        end

        return mask
    end


    function Compute_Reynolds_Number(grid::Grid, fluid::Fluid, lattice_dimensions::NamedTuple{(:Type, :Radius, :Position)}, mask::BitMatrix)      
        lattice_inflow_velocity = fluid.lattice_inflow_velocity
        lattice_viscosity = fluid.lattice_viscosity
        Characteristic_Length = lattice_dimensions.Radius

        lattice_Reynolds = (lattice_inflow_velocity .* Characteristic_Length)/lattice_viscosity

        return lattice_Reynolds
    end#Compute_Reynolds_Number

    function Compute_Reynolds_Number(grid::Grid, fluid::Fluid, lattice_dimensions::NamedTuple{(:Type, :Width, :Height, :Angle, :Position)}, mask::BitMatrix)
        gridY = grid.gridY

        lattice_inflow_velocity = fluid.lattice_inflow_velocity
        lattice_viscosity = fluid.lattice_viscosity
        Characteristic_Length = maximum(gridY[mask]) - minimum(gridY[mask])

        lattice_Reynolds = (lattice_inflow_velocity .* Characteristic_Length)/lattice_viscosity

        return lattice_Reynolds
    end#Compute_Reynolds_Number

    function Compute_Reynolds_Number(grid::Grid, fluid::Fluid, lattice_dimensions::NamedTuple{(:Type, :Width, :Height, :Angle, :Position, :Path)}, mask::BitMatrix)
        gridY = grid.gridY

        lattice_inflow_velocity = fluid.lattice_inflow_velocity
        lattice_viscosity = fluid.lattice_viscosity
        Characteristic_Length = maximum(gridY[mask]) - minimum(gridY[mask])

        lattice_Reynolds = (lattice_inflow_velocity .* Characteristic_Length)/lattice_viscosity

        return lattice_Reynolds
    end#Compute_Reynolds_Number

    function Compute_Reynolds_Number(grid::Grid, fluid::Fluid, lattice_dimensions::Nothing, mask::BitMatrix)
        gridY = grid.gridY

        lattice_inflow_velocity = fluid.lattice_inflow_velocity
        lattice_viscosity = fluid.lattice_viscosity
        Characteristic_Length = maximum(gridY) - minimum(gridY)

        lattice_Reynolds = (lattice_inflow_velocity .* Characteristic_Length)/lattice_viscosity

        return lattice_Reynolds
    end#Compute_Reynolds_Number
end#Object_SetUp