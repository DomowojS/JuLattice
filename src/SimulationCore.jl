module SimulationCore
using ..ConfigReader, .. Simulation_SetUp, ..Object_SetUp, ..Boundary_Conditions_SetUp
export Run_Simulation

function Run_Simulation(simulation::Simulation_Params, fluid::Fluid, grid::Grid, mutable_grid::Mutable_Grid, object::Geometry, boundary_conditions)
    time_steps = simulation.time_steps
    τ = simulation.τ
    velocity_vector_x = simulation.velocity_vector_x
    velocity_vector_y = simulation.velocity_vector_y
    weights = simulation.weights

    densityGrid = mutable_grid.densityGrid
    distributions = mutable_grid.distributions
    distributions_equilibrium= mutable_grid.distributions_equilibrium
    velocityX = mutable_grid.velocityX
    velocityY = mutable_grid.velocityY
    dotprod_velocities = mutable_grid.dotprod_velocities

    

   # Run Simulation Loop
    for i in 1:time_steps
        
        Get_Macroscopic_Values!(simulation, mutable_grid)
    
        Apply_Collision!(simulation, mutable_grid)
        
        Stream!(simulation, mutable_grid)
  
        Apply_Boundary_Conditions!(simulation, mutable_grid, grid, object, boundary_conditions)
  

            # Plot of the field
            if ((i % 10 == 0)) || (i == simulationTime)
                # Set velocities inside the cylinder to zero
                velocityX[cylinder] .= NaN
                velocityY[cylinder] .= NaN

                # Compute vorticity
                fill!(vorticity, 0.0)
                dv_dx = circshift(velocityY, (-1, 0)) .- circshift(velocityY, (1, 0))
                du_dy = circshift(velocityX, (0, -1)) .- circshift(velocityX, (0, 1))
                vorticity .= dv_dx .- du_dy
                vorticity[inlet] .= 0.0
                vorticity[outlet] .= 0.0

                # Mask the cylinder region
                vorticity[cylinder] .= NaN

                if ((i % 100 == 0)) || (i == simulationTime)
                    Log_Simulation_Runtime(i, simulationTime)
                end
                # Update the observables
                if Plotvorticity==true 
                    vorticity_obs[] = copy(vorticity) 
                    step_text[] = "Time step: $i, $(floor(Int, i*delta_t))s"
                end
                if Plotvx==true 
                    velocityX_obs[] = copy(velocityX) 
                    step_text_vx[] = "Time step: $i, $(floor(Int, i*delta_t))s"
                end
                if Plotvy==true 
                    velocityY_obs[] = copy(velocityY) 
                    step_text_vy[] = "Time step: $i, $(floor(Int, i*delta_t))s"
                end

                yield()
                sleep(0.05)
            end

    end 


    end#Run_Simulation

    function Get_Macroscopic_Values!(simulation::Simulation_Params, mutable_grid::Mutable_Grid)
        # Assign pointers
        velocity_vector_x = simulation.velocity_vector_x
        velocity_vector_y = simulation.velocity_vector_y

        densityGrid = mutable_grid.densityGrid
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
        mutable_grid.dotprod_velocities .= (velocity_vector_x .* velocityX) .+ (velocity_vector_y .* velocityY);
        mutable_grid.distributions_equilibrium .= weights .* densityGrid .*(1 .+ 3 .*mutable_grid.dotprod_velocities .+ 4.5 .*mutable_grid.dotprod_velocities.^2 .- 1.5 .*(velocityX.^2 .+ velocityY.^2));
        # Relax towards equilibrium
        mutable_grid.distributions .+= -(1/τ) .* (distributions .- mutable_grid.distributions_equilibrium);
    end#Apply_Collision

    function Stream!(simulation::Simulation_Params, mutable_grid::Mutable_Grid)
        #Assign pointers
        velocity_vector_x = simulation.velocity_vector_x
        velocity_vector_y = simulation.velocity_vector_y

        # Stream 
        for j in 1:Q
            mutable_grid.distributions[:,:,j] = circshift(mutable_grid.distributions[:,:,j], (velocity_vector_x[j], velocity_vector_y[j]))
        end
    
    end#Stream

    function Apply_Boundary_Conditions!(simulation::Simulation_Params, mutable_grid::Mutable_Grid, grid::Grid, object::Geometry, boundary_conditions)

        for side in (:Left, :Right, :Top, :Bottom, :Object)
            if side == :Object
                mask = object.mask
            else
                mask = getproperty(grid, Symbol(string(side), "_Boundary"))
            end

            if isdefined(boundary_conditions, side)
                Boundary_Condition!(mutable_grid, mask, boundary_conditions[side])
            end
        end

        ## Apply Boundary conditions
        #Inlet velocity bc (unknown: f_1, f_8, f_9)
        densityGrid[inlet, :] .= (sum(distributions[inlet, [1,3,5]], dims=2).+ 2 .*sum(distributions[inlet, [2,6,7]], dims=2)) ./ (1-lattice_inflow_velocity)
        distributions[inlet, 4] .= distributions[inlet, 2] .+ (2/3 .* densityGrid[inlet,:] .* lattice_inflow_velocity)
        distributions[inlet, 8] .= distributions[inlet, 6] .+ (1/6 .* densityGrid[inlet,:] .* lattice_inflow_velocity) .- (1/2 .* (distributions[inlet, 3] .- distributions[inlet, 5]))
        distributions[inlet, 9] .= distributions[inlet, 7] .+ (1/6 .* densityGrid[inlet,:] .* lattice_inflow_velocity) .+ (1/2 .* (distributions[inlet, 3] .- distributions[inlet, 5]))

    end#Apply_Boundary_Conditions

    function Boundary_Condition!(mutable_grid::Mutable_Grid, mask::BitMatrix, type::NoSlip)
        mutable_grid.distributions[mask, 1:Q] .= mutable_grid.distributions[mask, [1,4,5,2,3,8,9,6,7]];
    end#Boundary_Condition

    function Boundary_Condition!(mutable_grid::Mutable_Grid, mask::BitMatrix, object::Geometry, type::Velocity)

    end#Boundary_Condition

    function Boundary_Condition!(mutable_grid::Mutable_Grid, mask::BitMatrix, object::Geometry, type::Velocity)

    end#Boundary_Condition
    function Boundary_Condition!(mutable_grid::Mutable_Grid, mask::BitMatrix, object::Geometry, type::Velocity)

    end#Boundary_Condition

    function Boundary_Condition!(mutable_grid::Mutable_Grid, mask::BitMatrix, object::Geometry, type::Velocity)

    end#Boundary_Condition

    function Boundary_Condition!(mutable_grid::Mutable_Grid, mask::BitMatrix, object::Geometry, type::ZeroGradient)
        distributions[mask, [4, 8, 9]] .= distributions[gridlengthX-1, :, [4, 8, 9]]
    end#Boundary_Condition

    function Boundary_Condition!(mutable_grid::Mutable_Grid, mask::BitMatrix, object::Geometry, type::ZeroGradient)

    end#Boundary_Condition

    function Boundary_Condition!(mutable_grid::Mutable_Grid, mask::BitMatrix, object::Geometry, type::ZeroGradient)

    end#Boundary_Condition

    function Boundary_Condition!(mutable_grid::Mutable_Grid, mask::BitMatrix, object::Geometry, type::ZeroGradient)
    end#Boundary_Condition

end#SimulationCore