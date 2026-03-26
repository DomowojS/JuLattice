module Logger
export Log_Grid_Dimensions, Log_Fluid_Parameters, Log_Simulation_Header, 
        Log_Simulation_Runtime, Log_Simulation_Tail, Log_Simulation_Start


function Log_Simulation_Header()
    println("#############################")
    println("      Running JuLattice      ")
    println("#############################")
end#Log_Simulation_Header

function Log_Grid_Dimensions(nx::Int, ny::Int, nz::Int, delta_x, delta_t)
    println("Grid dimensions (in nodes): x=$nx, y=$ny, z=$nz") 
    println("Δx: $delta_x m")
    println("Δt: $delta_t s")
end#Log_Grid_Dimensions

function Log_Fluid_Parameters(τ, ω, u_phys, u_lat, Re_phys, Re_lat)
    println("##### Computed Fluid Values #####")
    println("   τ  = ", round(τ, digits=10))
    println("   ω = ", round(ω, digits=10))
    println("   u (phys.) = ", round(u_phys, digits=10))
    println("   u (lat.) = ", round(u_lat, digits= 10))
    println("Reynolds number check:")
    #Re_phys = Inflow_Velocity * 2 * Radius / Kinematic_Viscosity
    println("   Re (physical) = ", round(Re_phys, digits=2))
    println("   Re (lattice) = ", round(Re_lat, digits=2))
    println("#################################")
end#Log_Fluid_Parameters


function Log_Simulation_Start()
    println("#################################")
    println("Starting Simulation:")
end

function Log_Simulation_Runtime(i::Int64, simulationTime::Int64)
    println("Time Step: $i / $simulationTime")
end#Log_Simulation_Runtime

function Log_Simulation_Tail()
    println("Simulation finished.")
    println("#############################")
end#Log_Simulation_Tail

end#Logger