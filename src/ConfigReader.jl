module ConfigReader
using JSON
export Config, Create_Config_From_JSON

    #Type
    abstract type ConfigGroup end

    #Struct
    struct Config
        length_X::Real; 
        length_Y::Real; 
        Object_Type::String; 
        Object_Radius::Real;
        Object_Width::Real;
        Object_Height::Real;
        Object_Angle::Real; 
        Object_Position::Vector{Real}; 
        Fluid_Density::Real; Inflow_Velocity::Real; 
        Kinematic_Viscosity::Real; 
        Simulation_Time::Real; 
        delta_x::Real; 
        τ::Real; 
        Left_BC::String; 
        Right_BC::String;
        Top_BC::String;
        Bottom_BC::String;
        Left_BC_Velocity::Real;
        Right_BC_Velocity::Real;
        Top_BC_Velocity::Real;
        Bottom_BC_Velocity::Real;
        Plotvx::Bool; 
        Plotvy::Bool; 
        Plotvorticity::Bool;
    end#config

    struct Config_Domain <: ConfigGroup
        data::Dict{String, Any}
    end

    struct Config_Object <: ConfigGroup
        data::Dict{String, Any}
    end

    struct Config_Fluid <: ConfigGroup
        data::Dict{String, Any}
    end

    struct Config_Simulation <: ConfigGroup
        data::Dict{String, Any}
    end

    struct Config_BC <: ConfigGroup
        data::Dict{String, Any}
    end

    struct Config_Plot <: ConfigGroup
        data::Dict{String, Any}
    end

    function Create_Config_From_JSON(json_path::String)
        config_dictionary = Read_Config_From_JSON(json_path)

        config_struct     = Create_Struct_From_Dict(config_dictionary)

        return config_struct

    end#Create_Config_From_JSON
    
    function Read_Config_From_JSON(json_path::String)
        
        #Open Json file seperately - to prevent saving issues
        config_data = open(json_path, "r") do file
            JSON.parse(read(file, String))
        end

        config_objects = Dict{Symbol, Any}()
    
        # Required groups
        required_configs = [
            (:domain, Config_Domain, "domain"),
            (:fluid, Config_Fluid, "fluid"),
            (:simulation, Config_Simulation, "simulation")
        ]
    
        for (key, constructor, name) in required_configs
            if haskey(config_data, name)
                config_objects[key] = constructor(config_data[name])
                Verify_Settings(config_objects[key])
            else
                error("$name not defined as own group in JSON file")
            end
        end
    
        # Optional groups
        optional_configs = [
            (:object, Config_Object, "Object"),
            (:boundary_conditions, Config_BC, "Boundary Conditions"),
            (:plotting, Config_Plot, "Plotting")
        ]
    
        for (key, constructor, name) in optional_configs
            string_key = String(key)
            if haskey(config_data, string_key)
                config_objects[key] = constructor(config_data[string_key])
                Verify_Settings(config_objects[key])
            else
                @info "No $name defined in JSON file"
            end
        end
    
        return config_objects
    end#Read_Config_From_JSON

    function Create_Struct_From_Dict(config_dictionary::Dict{Symbol, Any})
        
        flattened_dictionary = Flatten_Dictionary(config_dictionary)
        
        config_struct = Struct_From_Dictionary(flattened_dictionary)

        return config_struct

    end#Create_Struct_From_Dict

    function Flatten_Dictionary(config_dictionary::Dict{Symbol, Any})
        defaults = Dict{String, Any}()

        defaults = Dict(
            :length_X => 1.0,
            :length_Y => 1.0,
            :Object_Type => "None",
            :Object_Radius => 0.0,
            :Object_Width => 0.0,
            :Object_Height => 0.0,
            :Object_Angle => 0.0,
            :Object_Position => [0.0, 0.0],
            :Fluid_Density => 1000.0,
            :Inflow_Velocity => 0.0,
            :Kinematic_Viscosity => 0.001,
            :Simulation_Time => 1000,
            :delta_x => 0.01,
            :τ => 0.65,
            :Left_BC => "NoSlip",
            :Right_BC => "NoSlip",
            :Top_BC => "NoSlip",
            :Bottom_BC => "NoSlip",
            :Left_BC_Velocity => 0.0,
            :Right_BC_Velocity => 0.0,
            :Top_BC_Velocity => 0.0,
            :Bottom_BC_Velocity => 0.0,
            :Plotvx => false,
            :Plotvy => false,
            :Plotvorticity => false
        )

        flattened_dictionary=copy(defaults)

        mappings = Dict(
            :domain => [
                ("length_X", :length_X),
                ("length_Y", :length_Y)
            ],
            :object => [
                ("type", :Object_Type),
                ("radius", :Object_Radius),
                ("width", :Object_Width),
                ("height", :Object_Height),
                ("rotation", :Object_Angle),
                ("position", :Object_Position)
            ],
            :fluid => [
                ("density", :Fluid_Density),
                ("inflow_velocity", :Inflow_Velocity),
                ("kinematic_viscosity", :Kinematic_Viscosity)
            ],
            :simulation => [
                ("time", :Simulation_Time),
                ("delta_x", :delta_x),
                ("τ", :τ)
            ],
            :plotting => [
                ("velocity_x", :Plotvx),
                ("velocity_y", :Plotvy),
                ("vorticity", :Plotvorticity)
            ]
        )
        
        # Process mappings for each section
        for (section, keymap) in mappings
            if haskey(config_dictionary, section)
                section_data = config_dictionary[section].data
                for (source_key, target_key) in keymap
                    if haskey(section_data, source_key)
                        flattened_dictionary[target_key] = section_data[source_key]
                    end
                end
            end
        end
        
        # Special handling for boundary conditions because of its nested structure
        if haskey(config_dictionary, :boundary_conditions)
            bc_data = config_dictionary[:boundary_conditions].data
            
            # Process each boundary (left, right, top, bottom)
            for (boundary, value) in bc_data
                # Capitalize first letter of boundary name (e.g., "left" -> "Left")
                boundary_prefix = uppercase(first(boundary)) * lowercase(boundary[2:end])
                
                # Handle boundary type
                if haskey(value, "type")
                    bc_key = Symbol(boundary_prefix * "_BC")
                    flattened_dictionary[bc_key] = value["type"]
                end
                
                # Handle velocity if present
                if haskey(value, "velocity")
                    velocity_key = Symbol(boundary_prefix * "_BC_Velocity")
                    flattened_dictionary[velocity_key] = value["velocity"]
                end
            end
        end
        
        return flattened_dictionary

    end#Flatten_Dictionary

    function Struct_From_Dictionary(flattened_dictionary::Dict{Symbol, Any})
        dict = flattened_dictionary

        config = Config(
            dict[:length_X],
            dict[:length_Y],
            dict[:Object_Type],
            dict[:Object_Radius],
            dict[:Object_Width],
            dict[:Object_Height],
            dict[:Object_Angle],
            dict[:Object_Position],
            dict[:Fluid_Density],
            dict[:Inflow_Velocity],
            dict[:Kinematic_Viscosity],
            dict[:Simulation_Time],
            dict[:delta_x],
            dict[:τ],
            dict[:Left_BC],
            dict[:Right_BC],
            dict[:Top_BC],
            dict[:Bottom_BC],
            dict[:Left_BC_Velocity],
            dict[:Right_BC_Velocity],
            dict[:Top_BC_Velocity],
            dict[:Bottom_BC_Velocity],
            dict[:Plotvx],
            dict[:Plotvy],
            dict[:Plotvorticity]
        )
        
        return config

    end#Struct_From_Dictionary

    function Verify_Settings(domain::Config_Domain)
        for key in ["length_X", "length_Y"]
            if !haskey(domain.data, key)
                error("$key of the domain not defined in JSON config.")
            end

            val = domain.data[key]

            if !(isa(val, Real) && val > 0)
                error("Value for key '$key' must be a real number greater than 0, got: $val")
            end
        end
    end#Verify_Settings

    function Verify_Settings(object::Config_Object)
        data = object.data

        if !haskey(data, "type")
            error("Object must define a 'type' field.")
        end

        obj_type = data["type"]

        if obj_type == "Cylinder"
            required_keys = ["radius", "position"]
            for key in required_keys
                if !haskey(data, key)
                    error("Cylinder object must define '$key'")
                end
            end

            radius = data["radius"]
            position = data["position"]

            if !(isa(radius, Real) && radius > 0)
                error("'radius' must be a positive real number.")
            end

            if !(isa(position, AbstractVector) && length(position) == 2 &&
                all(x -> isa(x, Real) && x > 0, position))
                error("'position' must be a vector of two positive real numbers.")
            end

        elseif obj_type == "Rectangle"
            required_keys = ["width", "height", "position", "rotation"]
            for key in required_keys
                if !haskey(data, key)
                    error("Rectangle object must define '$key'")
                end
            end

            width    = data["width"]
            height   = data["height"]
            rotation = data["rotation"]
            position = data["position"]

            if !(isa(width, Real) && width > 0)
                error("'width' must be a positive real number.")
            end
            if !(isa(height, Real) && height > 0)
                error("'height' must be a positive real number.")
            end
            if !(isa(rotation, Real))
                error("'rotation' must be a real number.")
            end
            if !(isa(position, AbstractVector) && length(position) == 2 &&
                all(x -> isa(x, Real) && x > 0, position))
                error("'position' must be a vector of two positive real numbers.")
            end

        else
            @warn "Unsupported object type '$obj_type'. Only 'Cylinder' and 'Rectangle' are allowed. Object will be omitted."
        end
    end#Verify_Settings

    function Verify_Settings(fluid::Config_Fluid)
        required_keys = ["density", "kinematic_viscosity"]
        for key in required_keys
            if !haskey(fluid.data, key)
                error("Fluid object must define '$key'")
            end    
            
            val = fluid.data[key]

            if !(isa(val, Real) && val > 0)
            error("'$key' must be a positive real number.")
            end
        end

        if !haskey(fluid.data, "inflow_velocity")
            @info "No inflow velocity defined."
        end
    end#Verify_Settings

    function Verify_Settings(simulation::Config_Simulation)
        required_keys = ["time", "delta_x", "τ"]
        for key in required_keys
            if !haskey(simulation.data, key)
                error("Simulation object must define '$key'")
            end    
            
            val = simulation.data[key]

            if !(isa(val, Real) && val > 0)
            error("'$key' must be a positive real number.")
            end
        end

        if simulation.data["τ"] < 0.6 || simulation.data["τ"] > 1
            @warn "Recommended settings for relaxation time is 0.6 < τ < 1.0."
        end

    end#Verify_Settings

    function Verify_Settings(boundary_conditions::Config_BC)
        for key in ["left", "right", "top", "bottom"]
            if !haskey(boundary_conditions.data, key)
                @info "No boundary condition defined for $key side of the domain. Default: Periodic BC."
            else
                if boundary_conditions.data[key]["type"] == "Velocity" && !haskey(boundary_conditions.data[key], "velocity") 
                    error("No velocity defined for velocity boundary condition for $key side of the domain!")
                end
    
                allowed_types = ["Velocity", "ZeroGradient", "NoSlip"]
                if boundary_conditions.data[key]["type"] in allowed_types
                    # Valid type
                else
                    @warn "Invalid boundary condition type: $(boundary_conditions["data"][key]["type"]) for $key side of the domain.'"
                    @warn "Default, periodic BC will be used instead."
                end
            end
        end
    end#Verify_Settings

    function Verify_Settings(plot::Config_Plot)
        plotcounter = 0

        for key in ["velocity_x", "velocity_y", "vorticity"]
            if haskey(plot.data, key) && plot.data[key] == true
                @info "$key will be plotted" 
            end
            plotcounter += 1
        end

        if plotcounter == 0 
            @info "No plots requested"
        end
    end#Verify_Settings

    function ReadConfig(length_X::Real, length_Y::Real, Object_Type::String, Object_Radius::Real, Position::Vector{Float64}, Fluid_Density::Real, Inflow_Velocity::Real, Kinematic_Viscosity::Real, Simulation_Time::Real, delta_x::Real, τ::Real, Left_BC::String, Right_BC::String, Top_BC::String, Bottom_BC::String, Left_BC_Velocity::Real, Right_BC_Velocity::Real, Top_BC_Velocity::Real, Bottom_BC_Velocity::Real, Plotvx::Bool, Plotvy::Bool, Plotvorticity::Bool)
        config = Config(length_X::Real, length_Y::Real, Object_Type::String, Object_Radius::Real, Position::Vector{Float64}, Fluid_Density::Real, Inflow_Velocity::Real, Kinematic_Viscosity::Real, Simulation_Time::Real, delta_x::Real, τ::Real, Left_BC::String, Right_BC::String, Top_BC::String, Bottom_BC::String, Left_BC_Velocity::Real, Right_BC_Velocity::Real, Top_BC_Velocity::Real, Bottom_BC_Velocity::Real, Plotvx::Bool, Plotvy::Bool, Plotvorticity::Bool)
    return config

    end#ReadConfig

end#ConfigReader