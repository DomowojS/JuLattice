module Boundary_Conditions_SetUp
using ..ConfigReader, ..Simulation_SetUp, ..Object_SetUp
export Create_Boundary_Conditions, NoSlip, Velocity, ZeroGradient

    #Types
    abstract type BoundaryCondition end

    #Structs
    struct NoSlip <: BoundaryCondition
        Type::String
    end

    struct Velocity <: BoundaryCondition
        Type::String
        Physical_Velocity::Real
        Lattice_Velocity::Real
    end

    struct ZeroGradient <: BoundaryCondition
        Type::String
    end

    #Functions
    function Create_Boundary_Conditions(config::Config, fluid::Fluid, object::Geometry)
        # Merge all conditions to one tuple
        boundary_conditions = NamedTuple()

        boundary_conditions = BC_for_DomainWalls(boundary_conditions, config, fluid)
        boundary_conditions = BC_for_Object(boundary_conditions, object)
        return boundary_conditions
    end#Create_Boundary_Conditions

    function BC_for_DomainWalls(boundary_conditions, config::Config, fluid::Fluid)
        i=0
        for key in [:Left_BC, :Right_BC, :Top_BC, :Bottom_BC]
            i += 1
            if isdefined(config, key)
                varname = Symbol(replace(String(key), "_BC" => ""))
                type = getfield(config, key)
                
                # Create the appropriate struct based on the config
                instance = if type == "Velocity"
                    Velocity("Velocity", fluid.Inflow_Velocity, fluid.lattice_inflow_velocity)
                elseif type == "NoSlip"
                    NoSlip("NoSlip")
                elseif type == "ZeroGradient"
                    ZeroGradient("ZeroGradient")
                else
                    error("Unsupported $(key)_BC type:  $type")
                end
        
                boundary_conditions = merge(boundary_conditions, NamedTuple{(varname,)}((instance,)))
            end
        end

        if i>1
            return boundary_conditions
        end
    end#BC_for_DomainWalls

    function BC_for_Object(boundary_conditions, object::Geometry)
        if object.Type != "none"
            instance = NoSlip("NoSlip")
            boundary_conditions = merge(boundary_conditions, NamedTuple{(:Object,)}((instance,)))
            return boundary_conditions
        end
    end#BC_for_Object

end#Boundary_Conditions_SetUp