module Plotter
using GLMakie
using ..ConfigReader, ..Logger, ..Simulation_SetUp, ..Object_SetUp
export Create_Plot, Set_Up_Figures, Update_Figures!

## Create Plot 
    function Set_Up_Figures(config::Config, simulation::Simulation_Params, mutable_grid::Mutable_Grid)
        Plotvorticity=config.Plotvorticity 
        Plotvx=config.Plotvx
        Plotvy=config.Plotvy
        
        observables = Dict{Symbol, Any}()
        
        if Plotvorticity
            vorticity_obs, step_text, fig_vorticity = Create_Plot(simulation.gridlengthX, simulation.gridlengthY, mutable_grid.vorticity)
            screen1 = GLMakie.Screen()
            GLMakie.display(screen1, fig_vorticity)

            observables[:vorticity_obs] = vorticity_obs
            observables[:step_text_vorticity] = step_text
            observables[:screen1] = screen1
        end

        if Plotvx
            velocityX_obs, step_text_vx, fig_vx = Create_Plot(simulation.gridlengthX, simulation.gridlengthY, mutable_grid.velocityX, "X")
            screen2 = GLMakie.Screen()
            GLMakie.display(screen2, fig_vx)

            observables[:velocityX_obs] = velocityX_obs
            observables[:step_text_vx] = step_text_vx
            observables[:screen2] = screen2
        end

        if Plotvy
            velocityY_obs, step_text_vy, fig_vy = Create_Plot(simulation.gridlengthX, simulation.gridlengthY, mutable_grid.velocityY, "Y")
            screen3 = GLMakie.Screen()
            GLMakie.display(screen3, fig_vy)

            observables[:velocityY_obs] = velocityY_obs
            observables[:step_text_vy] = step_text_vy
            observables[:screen3] = screen3
        end

        return observables
    end#Set_Up_Figures

    function Create_Plot(gridlengthX::Int64, gridlengthY::Int64, vorticity::Matrix{Float64})
        ## Set up the figure and axis
        # Initialize Plot arrays
        vorticity_obs = Observable(vorticity);       
        # Set up the figure and axis with explicit sizing
        fig = Figure(size = (1000, 400))
        ax = Axis(fig[1, 1], aspect = DataAspect(), title = "Vorticity")

        # Display the vorticity field using a heatmap with dynamic color range
        hm = heatmap!(ax, 1:gridlengthX, 1:gridlengthY, vorticity_obs, 
                    colormap = :curl, 
                    nan_color = :black,
                    colorrange = (-0.2, 0.2))
        Colorbar(fig[1, 2], hm, label = "Lattice_Vorticity")
        rowsize!(fig.layout, 1, ax.scene.viewport[].widths[2])
        # Set axis limits explicitly
        xlims!(ax, 1, gridlengthX)
        ylims!(ax, 1, gridlengthY)

        # Create a text element for time step display
        step_text = Observable("Time step: 0, 0s")
        text!(ax, ceil(Int, (gridlengthX*2/100)), ceil(Int, (gridlengthY*2/100)), text = step_text, 
                color = :black, fontsize = 14)
                return vorticity_obs, step_text, fig
    end#Plot_Vorticity

    function Create_Plot(gridlengthX::Int64, gridlengthY::Int64, velocityX::Array{Float64, 2}, Direction::String)
                ## Set up the figure and axis
        # Initialize Plot arrays
        velocity_obs = Observable(velocityX);       
        # Set up the figure and axis with explicit sizing
        fig = Figure(size = (1000, 400))
        ax = Axis(fig[1, 1], aspect = DataAspect(), title = "Velocity_$Direction")

        # Display the vorticity field using a heatmap with dynamic color range
        hm = heatmap!(ax, 1:gridlengthX, 1:gridlengthY, velocity_obs, 
                    colormap = :inferno, 
                    nan_color = :black,
                    colorrange = (-0.2, 0.2))
        Colorbar(fig[1, 2], hm, label = "Lattice_Velocity_$Direction")
        rowsize!(fig.layout, 1, ax.scene.viewport[].widths[2])
        # Set axis limits explicitly
        xlims!(ax, 1, gridlengthX)
        ylims!(ax, 1, gridlengthY)

        # Create a text element for time step display
        step_text = Observable("Time step: 0, 0s")
        text!(ax, ceil(Int, (gridlengthX*2/100)), ceil(Int, (gridlengthY*2/100)), text = step_text, 
                color = :black, fontsize = 14)

                return velocity_obs, step_text, fig
    end#Plot_vx

## Update Plot
    function Update_Figures!(config::Config, simulation::Simulation_Params, mutable_grid::Mutable_Grid, object::Union{Cylinder, Rectangle}, plot_data::Dict{Symbol, Any}, i, time_steps)
            mutable_grid.velocityX[object.mask] .= NaN
            mutable_grid.velocityY[object.mask] .= NaN
            mutable_grid.vorticity[object.mask] .= NaN

            if ((i % 100 == 0)) || (i == time_steps)
                Log_Simulation_Runtime(i, time_steps)
            end

            # Update the observables
            if config.Plotvorticity
                Compute_Vorticity!(mutable_grid)            
                plot_data[:vorticity_obs][] = mutable_grid.vorticity
                plot_data[:step_text_vorticity][] = "Time step: $i, $(floor(Int, i*simulation.delta_t))s"
            end
            if config.Plotvx 
                plot_data[:velocityX_obs][] = mutable_grid.velocityX
                plot_data[:step_text_vx][] = "Time step: $i, $(floor(Int, i*simulation.delta_t))s"
            end
            if config.Plotvy 
                plot_data[:velocityY_obs][] = mutable_grid.velocityY 
                plot_data[:step_text_vy][] = "Time step: $i, $(floor(Int, i*simulation.delta_t))s"
            end

            yield()
            sleep(0.05)
    end#Update_Figures

    function Update_Figures!(config::Config, simulation::Simulation_Params, mutable_grid::Mutable_Grid, object::none, plot_data::Dict{Symbol, Any}, i, time_steps)

        if ((i % 100 == 0)) || (i == time_steps)
            Log_Simulation_Runtime(i, time_steps)
        end

        # Update the observables
        if config.Plotvorticity
            Compute_Vorticity!(mutable_grid)            
            plot_data[:vorticity_obs][] = mutable_grid.vorticity
            plot_data[:step_text_vorticity][] = "Time step: $i, $(floor(Int, i*simulation.delta_t))s"
        end
        if config.Plotvx 
            plot_data[:velocityX_obs][] = mutable_grid.velocityX
            plot_data[:step_text_vx][] = "Time step: $i, $(floor(Int, i*simulation.delta_t))s"
        end
        if config.Plotvy 
            plot_data[:velocityY_obs][] = mutable_grid.velocityY 
            plot_data[:step_text_vy][] = "Time step: $i, $(floor(Int, i*simulation.delta_t))s"
        end

        yield()
        sleep(0.05)
end#Update_Figures

    #Computations for plots
    function Compute_Vorticity!(mutable_grid::Mutable_Grid)
        dv_dx = circshift(mutable_grid.velocityY, (-1, 0)) .- circshift(mutable_grid.velocityY, (1, 0))
        du_dy = circshift(mutable_grid.velocityX, (0, -1)) .- circshift(mutable_grid.velocityX, (0, 1))
        mutable_grid.vorticity .= dv_dx .- du_dy  
    end#Compute_Vorticity

end#Plotter
