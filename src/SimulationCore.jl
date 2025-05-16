module SimulationCore
using ..ConfigReader, ..Simulation_SetUp, ..Object_SetUp, ..Boundary_Conditions_SetUp, ..Plotter
export Run_Simulation!

function Run_Simulation!(config::Config, simulation::Simulation_Params, fluid::Fluid, grid::Grid, mutable_grid::Mutable_Grid, object::Geometry, boundary_conditions::NamedTuple, plot_data::Dict{Symbol, Any})
    time_steps = simulation.time_steps   

   # Run Simulation Loop
    for i in 1:time_steps
        Get_Macroscopic_Values!(simulation, mutable_grid)
    
        Apply_Collision!(simulation, mutable_grid)
        
        Stream!(simulation, mutable_grid)
  
        Apply_Boundary_Conditions!(simulation, grid, mutable_grid, fluid, object, boundary_conditions)

        if ((i % 10 == 0)) || (i == time_steps)
        Update_Figures!(config, simulation, mutable_grid, object, plot_data, i, time_steps)
        end
    end 

    end#Run_Simulation

    function Get_Macroscopic_Values!(simulation::Simulation_Params, mutable_grid::Mutable_Grid)
        # Assign pointers
        velocity_vector_x = simulation.velocity_vector_x
        velocity_vector_y = simulation.velocity_vector_y

        distributions = mutable_grid.distributions

        # Compute Values
        mutable_grid.densityGrid = sum(distributions, dims=3);
        mutable_grid.velocityX .= (1 ./ mutable_grid.densityGrid) .* sum(distributions.*velocity_vector_x, dims=3); 
        mutable_grid.velocityY .= (1 ./ mutable_grid.densityGrid) .* sum(distributions.*velocity_vector_y, dims=3); 
    end#Get_Macroscopic_Values

    function Apply_Collision!(simulation::Simulation_Params, mutable_grid::Mutable_Grid)
        # Assign pointers
        densityGrid = mutable_grid.densityGrid
        distributions = mutable_grid.distributions
        velocityX = mutable_grid.velocityX
        velocityY = mutable_grid.velocityY

        velocity_vector_x = simulation.velocity_vector_x
        velocity_vector_y = simulation.velocity_vector_y
        weights = simulation.weights
        τ = simulation.τ

        # Compute equilibrium state
        mutable_grid.dotprod_velocity .= (velocity_vector_x .* velocityX) .+ (velocity_vector_y .* velocityY);
        mutable_grid.distributions_equilibrium .= weights .* densityGrid .*(1 .+ 3 .*mutable_grid.dotprod_velocity .+ 4.5 .*mutable_grid.dotprod_velocity.^2 .- 1.5 .*(velocityX.^2 .+ velocityY.^2));
        # Relax towards equilibrium
        mutable_grid.distributions .+= -(1/τ) .* (distributions .- mutable_grid.distributions_equilibrium);
    end#Apply_Collision

    function Stream!(simulation::Simulation_Params, mutable_grid::Mutable_Grid)
        #Assign pointers
        velocity_vector_x = simulation.velocity_vector_x
        velocity_vector_y = simulation.velocity_vector_y

        # Stream 
        for j in 1:simulation.Q
            mutable_grid.distributions[:,:,j] = circshift(mutable_grid.distributions[:,:,j], (velocity_vector_x[j], velocity_vector_y[j]))
        end
    
    end#Stream

    function Apply_Boundary_Conditions!(simulation::Simulation_Params, grid::Grid, mutable_grid::Mutable_Grid, fluid::Fluid, object::Geometry, boundary_conditions::NamedTuple)

        for side in (:Left, :Right, :Top, :Bottom, :Object)
            if side == :Object
                mask = object.mask
            else
                mask = getproperty(grid, Symbol(string(side), "_Boundary"))
            end

            if isdefined(boundary_conditions, side)
                Boundary_Condition!(simulation, grid, mutable_grid, fluid, mask, boundary_conditions[side], side)
            end
        end

    end#Apply_Boundary_Conditions

    function Boundary_Condition!(simulation::Simulation_Params, grid::Grid, mutable_grid::Mutable_Grid, fluid::Fluid, mask::BitMatrix, type::NoSlip, side::Symbol)
        mutable_grid.distributions[mask, 1:simulation.Q] .= mutable_grid.distributions[mask, [1,4,5,2,3,8,9,6,7]];
    end#Boundary_Condition

    function Boundary_Condition!(simulation::Simulation_Params, grid::Grid, mutable_grid::Mutable_Grid, fluid::Fluid, mask::BitMatrix, type::Velocity, side::Symbol)
        if side == :Left
        #Inlet velocity bc (unknown: f_4, f_8, f_9)
            mutable_grid.densityGrid[mask, :] .= (sum(mutable_grid.distributions[mask, [1,3,5]], dims=2).+ 2 .*sum(mutable_grid.distributions[mask, [2,6,7]], dims=2)) ./ (1-fluid.lattice_inflow_velocity)
            mutable_grid.distributions[mask, 4] .= mutable_grid.distributions[mask, 2] .+ (2/3 .* mutable_grid.densityGrid[mask,:] .* fluid.lattice_inflow_velocity)
            mutable_grid.distributions[mask, 8] .= mutable_grid.distributions[mask, 6] .+ (1/6 .* mutable_grid.densityGrid[mask,:] .* fluid.lattice_inflow_velocity) .- (1/2 .* (mutable_grid.distributions[mask, 3] .- mutable_grid.distributions[mask, 5]))
            mutable_grid.distributions[mask, 9] .= mutable_grid.distributions[mask, 7] .+ (1/6 .* mutable_grid.densityGrid[mask,:] .* fluid.lattice_inflow_velocity) .+ (1/2 .* (mutable_grid.distributions[mask, 3] .- mutable_grid.distributions[mask, 5]))
        elseif side == :Right 
            #Inlet velocity bc (unknown: f_2, f_6, f_7)
            mutable_grid.densityGrid[mask, :] .= (sum(mutable_grid.distributions[mask, [1,3,5]], dims=2).+ 2 .*sum(mutable_grid.distributions[mask, [4,8,9]], dims=2)) ./ (1-fluid.lattice_inflow_velocity)
            mutable_grid.distributions[mask, 2] .= mutable_grid.distributions[mask, 4] .+ (2/3 .* mutable_grid.densityGrid[mask,:] .* fluid.lattice_inflow_velocity)
            mutable_grid.distributions[mask, 6] .= mutable_grid.distributions[mask, 8] .+ (1/6 .* mutable_grid.densityGrid[mask,:] .* fluid.lattice_inflow_velocity) .- (1/2 .* (mutable_grid.distributions[mask, 5] .- mutable_grid.distributions[mask, 3]))
            mutable_grid.distributions[mask, 7] .= mutable_grid.distributions[mask, 9] .+ (1/6 .* mutable_grid.densityGrid[mask,:] .* fluid.lattice_inflow_velocity) .+ (1/2 .* (mutable_grid.distributions[mask, 5] .- mutable_grid.distributions[mask, 3]))
        elseif side == :Top 
            #Inlet velocity bc (unknown: f_5, f_6, f_9)
            mutable_grid.densityGrid[mask, :] .= (sum(mutable_grid.distributions[mask, [1,2,4]], dims=2).+ 2 .*sum(mutable_grid.distributions[mask, [3,7,8]], dims=2)) ./ (1-fluid.lattice_inflow_velocity)
            mutable_grid.distributions[mask, 5] .= mutable_grid.distributions[mask, 3] .+ (2/3 .* mutable_grid.densityGrid[mask,:] .* fluid.lattice_inflow_velocity)
            mutable_grid.distributions[mask, 6] .= mutable_grid.distributions[mask, 8] .+ (1/6 .* mutable_grid.densityGrid[mask,:] .* fluid.lattice_inflow_velocity) .- (1/2 .* (mutable_grid.distributions[mask, 2] .- mutable_grid.distributions[mask, 4]))
            mutable_grid.distributions[mask, 9] .= mutable_grid.distributions[mask, 7] .+ (1/6 .* mutable_grid.densityGrid[mask,:] .* fluid.lattice_inflow_velocity) .+ (1/2 .* (mutable_grid.distributions[mask, 2] .- mutable_grid.distributions[mask, 4]))
        elseif side == :Bottom
            #Inlet velocity bc (unknown: f_3, f_7, f_8)
            mutable_grid.densityGrid[mask, :] .= (sum(mutable_grid.distributions[mask, [1,2,4]], dims=2).+ 2 .*sum(mutable_grid.distributions[mask, [5,6,9]], dims=2)) ./ (1-fluid.lattice_inflow_velocity)
            mutable_grid.distributions[mask, 3] .= mutable_grid.distributions[mask, 5] .+ (2/3 .* mutable_grid.densityGrid[mask,:] .* fluid.lattice_inflow_velocity)
            mutable_grid.distributions[mask, 7] .= mutable_grid.distributions[mask, 9] .+ (1/6 .* mutable_grid.densityGrid[mask,:] .* fluid.lattice_inflow_velocity) .- (1/2 .* (mutable_grid.distributions[mask, 2] .- mutable_grid.distributions[mask, 4]))
            mutable_grid.distributions[mask, 8] .= mutable_grid.distributions[mask, 6] .+ (1/6 .* mutable_grid.densityGrid[mask,:] .* fluid.lattice_inflow_velocity) .+ (1/2 .* (mutable_grid.distributions[mask, 2] .- mutable_grid.distributions[mask, 4]))
        else
            error("Unsupported BC time for: $(side)")
        end
    end#Boundary_Condition

    function Boundary_Condition!(simulation::Simulation_Params, grid::Grid, mutable_grid::Mutable_Grid, fluid::Fluid, mask::BitMatrix, type::ZeroGradient, side::Symbol)
        if side == :Left
            mutable_grid.distributions[mask, [3, 6, 7]] .= mutable_grid.distributions[2, :, [3, 6, 7]]
        elseif side == :Right 
            mutable_grid.distributions[mask, [4, 8, 9]] .= mutable_grid.distributions[simulation.gridlengthX-1, :, [4, 8, 9]]
        elseif side == :Top 
            mutable_grid.distributions[mask, [3, 7, 8]] .= mutable_grid.distributions[:, simulation.gridlengthY-1, [3, 7, 8]]
        elseif side == :Bottom
            mutable_grid.distributions[mask, [5, 6, 9]] .= mutable_grid.distributions[:, 2, [5, 6, 9]]
        else
            error("Unsupported BC time for: $(side)")
        end
    end#Boundary_Condition

    function Boundary_Condition!(simulation::Simulation_Params, grid::Grid, mutable_grid::Mutable_Grid, fluid::Fluid, mask::BitMatrix, type::Periodic, side::Symbol)
        #Nothing happens for periodic (stream function cirshift already works periodic)
    end#Boundary_Condition
end#SimulationCore