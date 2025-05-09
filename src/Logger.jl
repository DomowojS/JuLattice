module Logger
export Log_Simulation_Runtime, Log_Discretization_Settings, Log_Simulation_Header, Log_Simulation_Tail

    function Log_Simulation_Runtime(i::Int64, simulationTime::Int64)
        println("Time Step: $i / $simulationTime")
    end#Log_Simulation_Runtime

    function Log_Discretization_Settings(delta_x::Float64, delta_t::Float64, lattice_Reynolds::Float64, gridX::Matrix{Int64})
        println("Δx: $delta_x m")
        println("Δt: $delta_t s")
        println("Re: $(floor(Int,lattice_Reynolds))")
        println("Number of lattice nodes: $(length(gridX))")
    end#Log_Discretization_Settings

    function Log_Simulation_Header()
        println("#############################")
        println("      Running JuLattice      ")
        println("#############################")
    end#Log_Simulation_Header

    function Log_Simulation_Tail()
        println("Simulation finished.")
        println("#############################")
    end#Log_Simulation_Tail

end#Logger