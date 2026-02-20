module Logger
export Log_Simulation_Runtime, Log_Discretization_Settings,
       Log_Simulation_Header, Log_Simulation_Start, Log_Simulation_Tail

function Log_Simulation_Header()
    println("#############################")
    println("      Running JuLattice      ")
    println("#############################")
end

function Log_Discretization_Settings(deltaX::Float64, deltaT::Float64, omegaBGK::Float64, Re::Int64)
    println("dx: $deltaX m")
    println("dt: $deltaT s")
    println("omegaBGK: $omegaBGK")
    println("Re: $Re")
end

function Log_Simulation_Start()
    println("#############################")
    println("      Starting Simulation    ")
    println("#############################")
end

function Log_Simulation_Runtime(i::Int64, nSteps::Int64, nups::Float64)
    mnups = nups / 1e6
    println("Step: $i / $nSteps  |  MNUPS: $(round(mnups, digits=2))")
end

function Log_Simulation_Tail()
    println("Simulation finished.")
    println("#############################")
end

end#Logger
