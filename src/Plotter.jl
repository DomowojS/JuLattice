module Plotter
using GLMakie
export Create_Plot, Create_Plot3D

## Set up Plot 
    function Create_Plot(gridlengthX::Int64, gridlengthY::Int64)
        ## Set up the figure and axis
        # Initialize Plot arrays
        vorticity = zeros(gridlengthX, gridlengthY);
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
        text_obj = text!(ax, ceil(Int, (gridlengthX*2/100)), ceil(Int, (gridlengthY*2/100)), text = step_text, 
                color = :black, fontsize = 14)
                return vorticity, vorticity_obs, text_obj, step_text, fig
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
        text_obj = text!(ax, ceil(Int, (gridlengthX*2/100)), ceil(Int, (gridlengthY*2/100)), text = step_text, 
                color = :black, fontsize = 14)

                return velocity_obs, text_obj, step_text, fig
    end#Plot_vx

    function Create_Plot3D(nx::Int, ny::Int, nz::Int, field3d::Array{<:Real,3}; title="Test123")
        # Initialize Plot arrays
        vol_obs = Observable(Float32.(field3d))
        # Set up the figure and axis with explicit sizing
        fig = Figure(size = (900, 700))
        ax  = Axis3(fig[1, 1], title = title)

        plt = volume!(ax, vol_obs; colormap = :turbo, transparency = true)
        Colorbar(fig[1, 2], plt, label = title)

        #Create a text element for time step display
        step_text = Observable("Time step: 0, 0s")
        text!(ax, 1, 1, nz, text = step_text, color = :white, fontsize = 22, align = (:left, :bottom))

        return vol_obs, step_text, fig
    end
## Update Plot

end#module
