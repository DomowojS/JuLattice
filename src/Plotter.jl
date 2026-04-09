module Plotter
using GLMakie
using Colors
using ColorSchemes

export Create_Plot_XY, Create_Plot_XZ, Create_Vorticity_XY, Create_Vorticity_XZ, Create_Plot_Mag_XY, Create_Plot_Mag_XZ,
       setup_vx_plot, setup_mag_plot, setup_debug_plots, update_plots!

## Custom Heatmap settings 
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

    function custom_rdbu_with_zero(n::Int=256)
        base = get(colorschemes[:RdBu], LinRange(0, 1, n))

        base[div(n,2)] = RGB(0,1,0)
        return cgrad(base)
    end
######################################################

## Plot setup functions
######################################################
function setup_vx_plot(gridlengthX, gridlengthY, gridlengthZ, velocityX, midY, midZ)
    fig = Figure(size = (1400, 1200))
    ax1 = Axis(fig[1,1])
    ax2 = Axis(fig[2,1])
    vx_xy_obs, step_text_vx_xy, hm1 = Create_Plot_XY(gridlengthX-2, gridlengthY-2,
                                                           velocityX[2:gridlengthX-1, 2:gridlengthY-1, midZ];
                                                           title="v_x at z=$(midZ)", ax=ax1)
    vx_xz_obs, step_text_vx_xz, hm2 = Create_Plot_XZ(gridlengthX-2, gridlengthZ-2,
                                                        velocityX[2:gridlengthX-1, midY, 2:gridlengthZ-1];
                                                        title="v_x at y=$(midY)", ax=ax2)
    Colorbar(fig[1, 2], hm1, label = "Lattice Velocity")
    Colorbar(fig[2, 2], hm2, label = "Lattice Velocity")
    display(fig)
    return vx_xy_obs, step_text_vx_xy, vx_xz_obs, step_text_vx_xz
end

function setup_mag_plot(gridlengthX, gridlengthY, gridlengthZ, velocityMag, midY, midZ)
    ny = gridlengthY - 2
    nz = gridlengthZ - 2
    base_height = 280
    row2_height = round(Int, base_height * nz / ny)

    fig_mag = Figure(size = (1000, base_height + row2_height + 30 + 150))
        ax1_mag = Axis(fig_mag[1,1])
        ax2_mag = Axis(fig_mag[2,1])
        mag_xy_obs, step_text_mag_xy, hm1_mag = Create_Plot_Mag_XY(gridlengthX-2, gridlengthY-2,
                                                                     velocityMag[2:gridlengthX-1, 2:gridlengthY-1, midZ];
                                                                     title="|v| at z=$(midZ)", ax=ax1_mag)
        mag_xz_obs, step_text_mag_xz, hm2_mag = Create_Plot_Mag_XZ(gridlengthX-2, gridlengthZ-2,
                                                                     velocityMag[2:gridlengthX-1, midY, 2:gridlengthZ-1];
                                                                     title="|v| at y=$(midY)", ax=ax2_mag)
        Colorbar(fig_mag[1, 2], hm1_mag, label = "Lattice Velocity Magnitude")
        Colorbar(fig_mag[2, 2], hm2_mag, label = "Lattice Velocity Magnitude")
        Label(fig_mag[3,1:2], text=step_text_mag_xy)
        rowsize!(fig_mag.layout, 1, Fixed(base_height))
        rowsize!(fig_mag.layout, 2, Fixed(row2_height))
        rowsize!(fig_mag.layout, 3, Fixed(30))
        colsize!(fig_mag.layout, 1, Auto(0.9))
        colsize!(fig_mag.layout, 2, Auto(0.1))
        rowgap!(fig_mag.layout, 0)
        display(fig_mag)
        return mag_xy_obs, step_text_mag_xy, mag_xz_obs, step_text_mag_xz
end

function setup_debug_plots(gridlengthX, gridlengthZ, velocityX, frontY, backY)
    vx_xz_front_obs, step_text_vx_xz_front, fig_front = Create_Plot_XZ(gridlengthX, gridlengthZ,
                                                                              velocityX[:, frontY, :];
                                                                              title="v_x at y=$(frontY) [FRONT]")
    screen_front = GLMakie.Screen()
    display(screen_front, fig_front)

    vx_xz_back_obs, step_text_vx_xz_back, fig_back = Create_Plot_XZ(gridlengthX, gridlengthZ,
                                                                        velocityX[:, backY, :];
                                                                        title="v_x at y=$(backY) [BACK]")
    screen_back = GLMakie.Screen()
    display(screen_back, fig_back)
    return vx_xz_front_obs, step_text_vx_xz_front, vx_xz_back_obs, step_text_vx_xz_back
end
######################################################

## update_plots
######################################################
function update_plots!(Plotmag, Plotvx, Plotdebug,
                        velocityMag, velocityX,
                        gridlengthX, gridlengthY, gridlengthZ, midY, midZ, frontY, backY,
                        i, simulationTime, delta_t,
                        mag_xy_obs=nothing, step_text_mag_xy=nothing,
                        mag_xz_obs=nothing, step_text_mag_xz=nothing,
                        vx_xy_obs=nothing, step_text_vx_xy=nothing,
                        vx_xz_obs=nothing, step_text_vx_xz=nothing,
                        vx_xz_front_obs=nothing, step_text_vx_xz_front=nothing,
                        vx_xz_back_obs=nothing, step_text_vx_xz_back=nothing)

    if Plotmag
        mag_xy_obs[] = velocityMag[2:gridlengthX-1, 2:gridlengthY-1, midZ]
        step_text_mag_xy[] = "Time step: $i / $simulationTime \n $(floor(Int, i*delta_t))s"
        mag_xz_obs[] = velocityMag[2:gridlengthX-1, midY, 2:gridlengthZ-1]
        step_text_mag_xz[] = "Time step: $i / $simulationTime \n $(floor(Int, i*delta_t))s"
    end
    if Plotvx
        vx_xy_obs[] = velocityX[2:gridlengthX-1, 2:gridlengthY-1, midZ]
        step_text_vx_xy[] = "Time step: $i / $simulationTime \n $(floor(Int, i*delta_t))s"
        vx_xz_obs[] = velocityX[2:gridlengthX-1, midY, 2:gridlengthZ-1]
        step_text_vx_xz[] = "Time step: $i / $simulationTime \n $(floor(Int, i*delta_t))s"
    end
    if Plotdebug
        vx_xz_back_obs[] = velocityX[:,backY,:]
        step_text_vx_xz_back[] = "Time step: $i / $simulationTime \n $(floor(Int, i*delta_t))s"
        vx_xz_front_obs[] = velocityX[:,frontY,:]
        step_text_vx_xz_front[] = "Time step: $i / $simulationTime \n $(floor(Int, i*delta_t))s"
    end
    
end
######################################################


## Plotting functions
######################################################
    function Create_Plot_XY(nx::Int, ny::Int, field2d::Array{<:Real, 2}; title="vx slice XY", ax=nothing)
        vx_obs = Observable(field2d)          
        colorrange = (-0.02, 0.02)
        hm = heatmap!(ax, 1:nx, 1:ny, vx_obs;
                        colormap = :RdBu,
                        nan_color = :black, colorrange = colorrange,
                        interpolate = false)
        xlims!(ax, 1, nx); ylims!(ax, 1, ny)
        ax.title = title
        ax.aspect = DataAspect()
        step_text = Observable("Time step: 0, 0s")
        text!(ax, 10, 10, text=step_text, color=:black, fontsize=14, align=(:left, :top))
        return vx_obs, step_text, hm
    end

    function Create_Plot_XZ(nx::Int, nz::Int, field2d::Array{<:Real, 2}; title="vx slice XZ", ax=nothing)
        vx_obs = Observable(field2d)
        colorrange = (-0.02, 0.02)
        hm = heatmap!(ax, 1:nx, 1:nz, vx_obs; 
                    colormap = :RdBu, 
                    nan_color= :black, colorrange = colorrange,
                    interpolate=false)
        xlims!(ax, 1, nx); ylims!(ax, 1, nz)
        ax.title = title
        ax.aspect = DataAspect()
        step_text = Observable("Time step: 0, 0s")
        text!(ax, 10, 10, text=step_text, color=:black, fontsize=14, align=(:left, :top))
        return vx_obs, step_text, hm    
    end

    function Create_Plot_Mag_XY(nx::Int, ny::Int, field2d::Array{<:Real, 2}; title="Vel. Magnitude slice XY", ax=nothing)
        vx_obs = Observable(field2d)          
        colorrange = (-0.02, 0.02)
        hm = heatmap!(ax, 1:nx, 1:ny, vx_obs;
                        colormap = :RdBu,
                        nan_color = :black, colorrange = colorrange,
                        interpolate = false)
        xlims!(ax, 1, nx); ylims!(ax, 1, ny)
        ax.title = title
        ax.aspect = DataAspect()
        step_text = Observable("Time step: 0 / 0, \n0s")
        #text!(ax, 10, 10, text=step_text, color=:black, fontsize=14, align=(:left, :top))
        return vx_obs, step_text, hm
    end

    function Create_Plot_Mag_XZ(nx::Int, nz::Int, field2d::Array{<:Real, 2}; title="Vel. Magnitude slice XZ", ax=nothing)
        vx_obs = Observable(field2d)
        colorrange = (-0.03, 0.03)
        hm = heatmap!(ax, 1:nx, 1:nz, vx_obs; 
                    colormap = :RdBu, 
                    nan_color= :black, colorrange = colorrange,
                    interpolate=false)
        xlims!(ax, 1, nx); ylims!(ax, 1, nz)
        ax.title = title
        ax.aspect = DataAspect()
        step_text = Observable("Time step: 0, 0s")
        #text!(ax, 10, 10, text=step_text, color=:black, fontsize=14, align=(:left, :top))
        return vx_obs, step_text, hm 
    end
######################################################
    

end#module
