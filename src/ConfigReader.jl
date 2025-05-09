module ConfigReader
export ReadConfig, Verification, Config

#Struct
struct Config
    length_X::Real; 
    length_Y::Real; 
    Object_Type::String; 
    Object_Radius::Real; 
    Object_Position::Vector{Float64}; 
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