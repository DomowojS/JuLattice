module ConfigReader
using JSON
export ReadConfig, Verification, Config, Read_Config_From_JSON

#Type
abstract type ConfigGroup end

#Struct
struct Config
    length_X::Real; 
    length_Y::Real; 
    Object_Type::String; 
    Object_Radius::Real; 
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


    function Read_Config_From_JSON(json_path::String)
        config_data = JSON.parsefile(json_path)

        # Required groups
        for (key, constructor, name) in [
            ("domain",      Config_Domain,      "Domain"),
            ("fluid",       Config_Fluid,       "Fluid"),
            ("simulation",  Config_Simulation,  "Simulation")
            ]
            if haskey(config_data, key)
                @eval $(Symbol(key)) = $constructor(config_data[$key])
                Verify_Settings(@eval $(Symbol(s)))
            else
                error("$name not defined as own group '$key' in JSON file")
            end
        end

        # Optional groups
        for (key, constructor, name) in [
            ("object",             Config_Object, "Object"),
            ("boundary_conditions", Config_BC,    "Boundary Conditions"),
            ("plotting",           Config_Plot,   "Plotting")
            ]   
            if haskey(config_data, key)
                @eval $(Symbol(key)) = $constructor(config_data[$key])
                Verify_Settings(@eval $(Symbol(s)))
            else
                @info "No $name defined in JSON file"
            end
        end

        

    end#Read_Config_From_JSON

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
            required_keys = ["Radius", "Position"]
            for key in required_keys
                if !haskey(data, key)
                    error("Cylinder object must define '$key'")
                end
            end

            radius = data["Radius"]
            position = data["Position"]

            if !(isa(radius, Real) && radius > 0)
                error("'Radius' must be a positive real number.")
            end

            if !(isa(position, AbstractVector) && length(position) == 2 &&
                all(x -> isa(x, Real) && x > 0, position))
                error("'Position' must be a vector of two positive real numbers.")
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
                error("'Position' must be a vector of two positive real numbers.")
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
            
            if !(isa(key, Real) && radius > 0)
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
            
            if !(isa(key, Real) && radius > 0)
            error("'$key' must be a positive real number.")
            end
        end

        if simulation.data["τ"] < 0.6 || simulation.data[τ] > 1
            @warn "Recommended settings for relaxation time is 0.6 < τ < 1.0."
        end

    end#Verify_Settings

function Verify_Settings(boundary_conditions::Config_BC)
    for key in ["left", "right", "Top", "Bottom"]
        if !haskey(boundary_conditions.data, key)
            @info "No boundary condition defined for $key side of the domain. Default: Periodic BC."
        end

        if boundary_conditions.data[key]["type"] == "Velocity" && !haskey(boundary_conditions.data[key], "velocity") 
            error("No velocity defined for velocity boundary condition for $key side of the domain!")
        end

        allowed_types = ["Velocity", "ZeroGradient", "NoSlip"]
        if boundary_conditions["data"][key]["type"] in allowed_types
            # Valid type
        else
            @warn "Invalid boundary condition type: $(boundary_conditions["data"][key]["type"]) for $key side of the domain.'"
            @warn "Default, periodic BC will be used instead."
        end
    end
end#Verify_Settings

    function Verify_Settings(plot::Config_Plot)
        plotcounter += 0

        for key in ["Velocity_x", "Velocity_y", "Vorticity"]
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

    function Verification(length_X::Real, length_Y::Real, Object_Type::String, Object_Radius::Real, Position::Vector{Float64}, Fluid_Density::Real, Inflow_Velocity::Real, Kinematic_Viscosity::Real, Simulation_Time::Real, delta_x::Real, τ::Real, Left_BC::String, Right_BC::String, Top_BC::String, Bottom_BC::String, Left_BC_Velocity::Real, Right_BC_Velocity::Real, Top_BC_Velocity::Real, Bottom_BC_Velocity::Real, Plotvx::Bool, Plotvy::Bool, Plotvorticity::Bool)
        #Object
            #Type Allowed? Rectangular/ Circular
            #Right Settings for type? Radius Position ...
            #Position within domain?
                #-> Error if no, warning if close to domain
            #Compute Reynolds number
                #If above 100 warning 
            #If no object --> Set object definition to 0 and print info
        #τ --> If not between 0.6 and 2 warning. If above 1.0 warning.
        #BC --> Copy from main file
    end#Verification


end#ConfigReader