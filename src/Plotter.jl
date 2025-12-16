module Plotter
using GLMakie
using Colors
export Create_Plot_XY, Create_Plot_XZ, Create_Vorticity_XY, Create_Vorticity_XZ

## Set up Plot 
######################################################

    function transparency_map(n::Int=256)
        t = range(-1, 1; length = n)
        alpha = 0.05 .+ 0.95 .* abs.(t)
        Makie.cgrad(:balance, n, alpha=alpha)
    end

    function transparency_map_vorticity(n::Int=256)
        t=range(0, 1; length=n)
        alpha = 0.15 .+ 0.85 .* (t .^ 0.6)
        Makie.cgrad(:turbo, n, alpha=alpha)
    end
######################################################

    function Create_Plot_XY(nx::Int, ny::Int, field2d::Array{<:Real, 2}; title="vx slice XY")
        vx_obs = Observable(field2d)
        fig = Figure(size = (900, 400))
        ax = Axis(fig[1,1], aspect = DataAspect(), title = title)
        colorrange =(-0.2, 0.2)
        hm = heatmap!(ax, 1:nx, 1:ny, vx_obs; 
                        colormap = transparency_map(), 
                        nan_color= :black, colorrange = colorrange, 
                        transparency=true)

        Colorbar(fig[1, 2], hm, label = "Lattice_Velocity")
        xlims!(ax, 1, nx); ylims!(ax, 1, ny)
        step_text = Observable("Time step: 0, 0s")
        text!(ax, ceil(Int, nx*0.05), ceil(Int, ny*0.05), text = step_text, color = :black, fontsize = 14)
        return vx_obs, step_text, fig
    end

    function Create_Plot_XZ(nx::Int, nz::Int, field2d::Array{<:Real, 2}; title="vx slice XZ")
        vx_obs = Observable(field2d)
        fig = Figure(size = (900, 400))
        ax = Axis(fig[1,1], aspect = DataAspect(), title = title)
        hm = heatmap!(ax, 1:nx, 1:nz, vx_obs; colormap = :inferno, nan_color= :white, colorrange = (-0.2, 0.2))
        Colorbar(fig[1, 2], hm, label = "Lattice_Velocity")
        xlims!(ax, 1, nx); ylims!(ax, 1, nz)
        step_text = Observable("Time step: 0, 0s")
        text!(ax, ceil(Int, nx*0.05), ceil(Int, nz*0.05), text = step_text, color = :black, fontsize = 14)
        return vx_obs, step_text, fig
    end

    function Create_Vorticity_XY(nx::Int, ny::Int, vorticity_XY::Array{<:Real, 2}; title="|ω| slice XY")
        omega_obs = Observable(vorticity_XY)
        fig = Figure(size = (900, 400))
        ax = Axis(fig[1,1], aspect = DataAspect(), title = title)
        hm = heatmap!(ax, 1:nx, 1:ny, omega_obs;
                      colormap = transparency_map_vorticity(),
                      nan_color = :black,
                      colorrange = (0.0, 0.05),
                      transparency = true)
        Colorbar(fig[1,2], hm, label = "|ω| (lattice)")
        xlims!(ax, 1, nx); ylims!(ax, 1, ny)
        step_text = Observable("Time step: 0, 0s")
        text!(ax, ceil(Int, nx*0.05), ceil(Int, ny*0.05), text = step_text, color = :black, fontsize = 14)
        
        return omega_obs, step_text, fig
    end

    function Create_Vorticity_XZ(nx::Int, nz::Int, vorticity_XZ::Array{<:Real, 2}; title="|ω| slice XZ")
        omega_obs = Observable(vorticity_XZ)
        fig = Figure(size = (900, 400))
        ax = Axis(fig[1,1], aspect = DataAspect(), title = title)
        hm = heatmap!(ax, 1:nx, 1:nz, omega_obs;
                      colormap = transparency_map_vorticity(),
                      nan_color = :black,
                      colorrange = (0.0, 0.05),
                      transparency = true)
        Colorbar(fig[1,2], hm, label = "|ω| (lattice)")
        xlims!(ax, 1, nx); ylims!(ax, 1, nz)
        step_text = Observable("Time step: 0, 0s")
        text!(ax, ceil(Int, nx*0.05), ceil(Int, nz*0.05), text = step_text, color = :black, fontsize = 14)
        
        return omega_obs, step_text, fig
    end









    # function Create_Plot3D(nx::Int, ny::Int, nz::Int, field3d::Array{<:Real,3}; title="Test123")
    #     # Initialize Plot arrays
    #     vol_obs = Observable(Float32.(field3d))
    #     # Set up the figure and axis with explicit sizing
    #     fig = Figure(size = (900, 700))
    #     ax  = Axis3(fig[1, 1], title = title)

    #     plt = volume!(ax, vol_obs; colormap = :turbo, transparency = true, colorrange = (0f0, 50f0))
    #     Colorbar(fig[1, 2], plt, label = title)

    #     #Create a text element for time step display
    #     step_text = Observable("Time step: 0, 0s")
    #     text!(ax, 1, 1, nz, text = step_text, color = :white, fontsize = 22, align = (:left, :bottom))

    #     return vol_obs, step_text, fig
    # end

    #############################################################################
    # function Create_Plot(gridlengthX::Int64, gridlengthY::Int64)
    #     ## Set up the figure and axis
    #     # Initialize Plot arrays
    #     vorticity = zeros(gridlengthX, gridlengthY);
    #     vorticity_obs = Observable(vorticity);       
    #     # Set up the figure and axis with explicit sizing
    #     fig = Figure(size = (1000, 400))
    #     ax = Axis(fig[1, 1], aspect = DataAspect(), title = "Vorticity")

    #     # Display the vorticity field using a heatmap with dynamic color range
    #     hm = heatmap!(ax, 1:gridlengthX, 1:gridlengthY, vorticity_obs, 
    #                 colormap = :curl, 
    #                 nan_color = :black,
    #                 colorrange = (-0.2, 0.2))
    #     Colorbar(fig[1, 2], hm, label = "Lattice_Vorticity")
    #     rowsize!(fig.layout, 1, ax.scene.viewport[].widths[2])
    #     # Set axis limits explicitly
    #     xlims!(ax, 1, gridlengthX)
    #     ylims!(ax, 1, gridlengthY)

    #     # Create a text element for time step display
    #     step_text = Observable("Time step: 0, 0s")
    #     text_obj = text!(ax, ceil(Int, (gridlengthX*2/100)), ceil(Int, (gridlengthY*2/100)), text = step_text, 
    #             color = :black, fontsize = 14)
    #             return vorticity, vorticity_obs, text_obj, step_text, fig
    # end#Plot_Vorticity

    # function Create_Plot(gridlengthX::Int64, gridlengthY::Int64, velocityX::Array{Float64, 2}, Direction::String)
    #             ## Set up the figure and axis
    #     # Initialize Plot arrays
    #     velocity_obs = Observable(velocityX);       
    #     # Set up the figure and axis with explicit sizing
    #     fig = Figure(size = (1000, 400))
    #     ax = Axis(fig[1, 1], aspect = DataAspect(), title = "Velocity_$Direction")

    #     # Display the vorticity field using a heatmap with dynamic color range
    #     hm = heatmap!(ax, 1:gridlengthX, 1:gridlengthY, velocity_obs, 
    #                 colormap = :inferno, 
    #                 nan_color = :black,
    #                 colorrange = (-0.2, 0.2))
    #     Colorbar(fig[1, 2], hm, label = "Lattice_Velocity_$Direction")
    #     rowsize!(fig.layout, 1, ax.scene.viewport[].widths[2])
    #     # Set axis limits explicitly
    #     xlims!(ax, 1, gridlengthX)
    #     ylims!(ax, 1, gridlengthY)

    #     # Create a text element for time step display
    #     step_text = Observable("Time step: 0, 0s")
    #     text_obj = text!(ax, ceil(Int, (gridlengthX*2/100)), ceil(Int, (gridlengthY*2/100)), text = step_text, 
    #             color = :black, fontsize = 14)

    #             return velocity_obs, text_obj, step_text, fig
    # end#Plot_vx

## Update Plot

end#module
