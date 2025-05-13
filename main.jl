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

    using .Logger, .ConfigReader, .Simulation_SetUp, .Object_SetUp, .Boundary_Conditions_SetUp, .Plotter, .SimulationCore

    function run(json_path::String)

        Log_Simulation_Header()

        config = Create_Config_From_JSON(json_path)

        simulation, fluid   = Set_Simulation_Params(config)

        grid                = Create_Grid(simulation)
        mutable_grid        = Initialize_Distributions(simulation, fluid, grid)

        object              = Create_Object(config, simulation, fluid, grid)

        boundary_conditions = Create_Boundary_Conditions(config, fluid, object)

        Log_Discretization_Settings(simulation.delta_x, simulation.delta_t, object.lattice_Reynolds, grid.gridX)

        plot_data = Set_Up_Figures(config, simulation, mutable_grid)

        Log_Simulation_Start()

        Run_Simulation!(config, simulation, fluid, grid, mutable_grid, object, boundary_conditions, plot_data)

        Log_Simulation_Tail()
    end#run
end#JuLattice