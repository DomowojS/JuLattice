module Simulation_SetUp
export Set_Simulation_Params

    function Set_Simulation_Params(config::Config)
        #Fixed values and read from config
        delta_x = config.delta_x 
        τ       = config.τ
        lattice_speedOfSound = 1 / √3;
        fluiddensity         = 100;

        delta_t = Compute_Time_Step(config.delta_x, τ, config.Kinematic_Viscosity, lattice_speedOfSound)
        lattice_inflow_velocity, lattice_viscosity, time_steps, gridlengthX, gridlengthY = 
            Convert_to_Lattice_Units(delta_t, delta_x, config.Inflow_Velocity, lattice_speedOfSound, τ, 
                                        config.Simulation_time, config.length_X, config.length_Y)

        # Hardcoded simulation parameters for D2Q9 LBM
        Q   = 9;
        velocity_vector = [     [0, -1, 0, 1, 0, -1, -1, 1, 1],
                                [0, 0, 1, 0, -1, -1, 1, 1, -1]];
        velocity_vector_x = reshape(velocity_vector[1,],1,1,Q)
        velocity_vector_y = reshape(velocity_vector[2,],1,1,Q)
        weights =   [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
        weights =   reshape(weights,1,1,Q);


        simulation = Simulation_Params(delta_x, delta_t, time_steps, gridlengthX, gridlengthY, velocity_vector_x, velocity_vector_y, weights);
        fluid      = Fluid(config.Kinematic_Viscosity, config.Fluid_Density, config.Inflow_Velocity, fluiddensity, lattice_inflow_velocity, lattice_viscosity);
        
        if 

    
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








    struct Simulation_Params
        delta_x::Real;
        delta_t::Real;
        time_steps::Int64;
        gridlengthX::Int64;
        gridlengthY::Int64;
        Q::Int64;
        velocity_vector_x::Array{Int64, 3};
        velocity_vector_y::Array{Int64, 3};
        weights::Array{Int64, 3};


    end#simulation_params

    struct Fluid
        Kinematic_Viscosity::Real
        Fluid_Density::Real
        Inflow_Velocity::Real
        fluiddensity::Int64;    
        lattice_inflow_velocity::Float64;
        lattice_viscosity::Float64;
    end#fluid









    function Create_Object(Type::String, Radius::Real, Position::Vector{Float64})
        lattice_Re
        struct Cylinder_Object
            lattice_Re
        end#cylinder_object

    end#Create_Object

    function Create_Object(Type::String, Width::Real, Height::Real, Angle::Real, Position::Vector{Float64})
        
        struct Rectangle_Object
        end#cylinder_object
        
    end#Create_Object




end#Simulation_SetUp