module IO
using Printf
export Log_Simulation_Header, Log_Discretization_Settings,
       Log_Simulation_Start, Log_Simulation_Runtime, Log_Simulation_Tail,
       Save_Forces!

function Log_Simulation_Header()
    println("#############################")
    println("      Running JuLattice      ")
    println("#############################")
end

function Log_Discretization_Settings(deltaX::Float64, deltaT::Float64, omegaBGK::Float64, Re::Int64, u::Float64)
    println("dx: $deltaX m")
    println("dt: $deltaT s")
    println("omegaBGK: $omegaBGK")
    println("Re: $Re")
    println("Inflow velocity: $u m/s")
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

function Save_Forces!(times::Vector{Float64}, cds::Vector{Float64}, cls::Vector{Float64};
                      dir::String = "./output")
    isdir(dir) || mkpath(dir)
    path = joinpath(dir, "forces.txt")
    open(path, "w") do io
        println(io, "# time[s]    cL    cD")
        for k in eachindex(times)
            @printf(io, "%.6e  %.6e  %.6e\n", times[k], cls[k], cds[k])
        end
    end
end

end # module IO
