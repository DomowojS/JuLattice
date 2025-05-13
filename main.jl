module JuLattice
    ############################
    ## Main file for JuLattice #
    ############################
    include("src/Logger.jl")
    include("src/ConfigReader.jl")
    include("src/Simulation_SetUp.jl")
    include("src/Object_SetUp.jl")
    include("src/Boundary_Conditions_SetUp.jl")    
    include("src/Plotter.jl")
    include("src/SimulationCore.jl")


    using Revise
    using .Logger, .ConfigReader, .Simulation_SetUp, .Object_SetUp, .Boundary_Conditions_SetUp, .Plotter, .SimulationCore

    function run(json_path::String)

        config = Read_Config_From_JSON(json_path)

        ## User Settings
        # Domain Settings
        length_X = 4;              # m
        length_Y = 1;              # m 

        # Object definition
        Object = "Cylinder";
        Radius   = 0.1    # m
        Position = [1, 0.5] # m

        # Fluid Settings 
        Fluid_Density = 1000.0;       # kg/m^3
        Inflow_Velocity = 0.4;      # m/s
        Kinematic_Viscosity = 0.001; # m^2/s 

        # Simulation Settings
        Simulation_Time = 8000;     # s
        delta_x = 0.01;             # discretisation in time and space
        τ = 0.65;

        # Boundary conditions
        Left_BC     = "Velocity"
        Left_BC_Velocity = Inflow_Velocity

        Right_BC    = "ZeroGradient"
        Top_BC      = "NoSlip"
        Bottom_BC   = "NoSlip"

        # Plot Requests
        Plotvx = false;
        Plotvy = false;
        Plotvorticity = true;

        config = ReadConfig(length_X, length_Y, Object, Radius, Position, Fluid_Density, Inflow_Velocity, 
                            Kinematic_Viscosity, Simulation_Time, delta_x, τ, 
                            Left_BC, Right_BC, Top_BC, Bottom_BC, Left_BC_Velocity, Right_BC_Velocity, Top_BC_Velocity, Bottom_BC_Velocity, 
                            Plotvx, Plotvy, Plotvorticity)

        Log_Simulation_Header()

        simulation, fluid   = Set_Simulation_Params(config)

        grid                = Create_Grid(simulation)
        mutable_grid        = Initialize_Distributions(simulation, fluid, grid)

        object              = Create_Object(config, simulation, fluid, grid)

        boundary_conditions = Create_Boundary_Conditions(config, fluid, object)

        Log_Discretization_Settings(simulation.delta_x, simulation.delta_t, object.lattice_Reynolds, grid.gridX)

        plot_data = Set_Up_Figures(config, simulation, mutable_grid)

        Log_Simuation_Start()

        Run_Simulation!(config, simulation, fluid, grid, mutable_grid, object, boundary_conditions, plot_data)

        Log_Simulation_Tail()
    end#run
end#JuLattice