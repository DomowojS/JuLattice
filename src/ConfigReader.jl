module ConfigReader
export ReadConfig

    function ReadConfig(length_X::Real, length_Y::Real, Object_Type::String, Object_Radius::Real, Position::Vector{Float64}, Fluid_Density::Real, Inflow_Velocity::Real, Kinematic_Viscosity::Real, Simulation_Time::Real, delta_x::Real, τ::Real, Left_BC::String, Right_BC::String, Top_BC::String, Bottom_BC::String, Left_BC_Velocity::Real, Right_BC_Velocity::Real, Top_BC_Velocity::Real, Bottom_BC_Velocity::Real, Plotvx::Bool, Plotvy::Bool, Plotvorticity::Bool)
        config = Config(length_X::Real, length_Y::Real, Object_Type::String, Object_Radius::Real, Position::Vector{Float64}, Fluid_Density::Real, Inflow_Velocity::Real, Kinematic_Viscosity::Real, Simulation_Time::Real, delta_x::Real, τ::Real, Left_BC::String, Right_BC::String, Top_BC::String, Bottom_BC::String, Left_BC_Velocity::Real, Right_BC_Velocity::Real, Top_BC_Velocity::Real, Bottom_BC_Velocity::Real, Plotvx::Bool, Plotvy::Bool, Plotvorticity::Bool)
    return config

    end#ReadConfig

    struct Config
        length_X::Real; 
        length_Y::Real; 
        Object_Type::String; 
        Object_Radius::Real; 
        Position::Vector{Float64}; 
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

end#ConfigReader