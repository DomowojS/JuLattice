module Plotter
using GLMakie
export Create_Plot, Update_Plot!, Create_Force_Plot, Update_Force_Plot!

function Create_Plot(Nx::Int, Ny::Int, plotU::Bool, plotV::Bool, plotVorticity::Bool,
                     rangeU::Tuple{Float64,Float64},
                     rangeV::Tuple{Float64,Float64},
                     rangeVort::Tuple{Float64,Float64},
                     deltaX::Float64)
    nplots = Int(plotU) + Int(plotV) + Int(plotVorticity)
    fig = Figure(size = (900, 280 * max(nplots, 1)))

    # Physical coordinates: ghost nodes at index 1 sit at -deltaX, fluid domain starts at 0
    xs = range(-deltaX, step = deltaX, length = Nx)
    ys = range(-deltaX, step = deltaX, length = Ny)
    text_x = 0.02 * (Nx - 2) * deltaX
    text_y = 0.93 * (Ny - 2) * deltaX

    step_text = Observable("Step: 0  |  t = 0.00 s")
    obs_u    = nothing
    obs_v    = nothing
    obs_vort = nothing
    row = 1

    if plotU
        ax = Axis(fig[row, 1], title = "U  [m/s]", aspect = DataAspect(),
                  xlabel = "x [m]", ylabel = "y [m]")
        obs_u = Observable(zeros(Nx, Ny))
        hm = heatmap!(ax, xs, ys, obs_u, colormap = :inferno, colorrange = rangeU, nan_color = :dimgray)
        Colorbar(fig[row+1, 1], hm, vertical = false)
        text!(ax, text_x, text_y, text = step_text, color = :dimgray, fontsize = 11)
        row += 2
    end

    if plotV
        ax = Axis(fig[row, 1], title = "V  [m/s]", aspect = DataAspect(),
                  xlabel = "x [m]", ylabel = "y [m]")
        obs_v = Observable(zeros(Nx, Ny))
        hm = heatmap!(ax, xs, ys, obs_v, colormap = :inferno, colorrange = rangeV, nan_color = :dimgray)
        Colorbar(fig[row+1, 1], hm, vertical = false)
        row += 2
    end

    if plotVorticity
        ax = Axis(fig[row, 1], title = "Vorticity  [1/s]", aspect = DataAspect(),
                  xlabel = "x [m]", ylabel = "y [m]")
        obs_vort = Observable(zeros(Nx, Ny))
        hm = heatmap!(ax, xs, ys, obs_vort, colormap = :curl, colorrange = rangeVort, nan_color = :dimgray)
        Colorbar(fig[row+1, 1], hm, vertical = false)
        row += 2
    end

    return fig, obs_u, obs_v, obs_vort, step_text
end

function Update_Plot!(obs_u, obs_v, obs_vort, step_text,
                      velocityX, velocityY,
                      i, deltaT, deltaX,
                      plotU::Bool, plotV::Bool, plotVorticity::Bool,
                      isFluid::BitMatrix, isObject::BitMatrix)

    velScale  = deltaX / deltaT          # [lu/ts] -> [m/s]
    vortScale = 1.0 / (2.0 * deltaT)    # central-diff (no /2) -> [1/s]

    step_text[] = "Step: $i  |  t = $(round(i*deltaT, digits=2)) s"

    if plotU && obs_u !== nothing
        data = velocityX .* velScale
        data[isObject] .= NaN
        obs_u[] = data
    end

    if plotV && obs_v !== nothing
        data = velocityY .* velScale
        data[isObject] .= NaN
        obs_v[] = data
    end

    if plotVorticity && obs_vort !== nothing
        vort = fill(NaN, size(velocityX))
        for I in CartesianIndices(isFluid)
            if isFluid[I]
                ix, iy = I.I
                vort[I] = ((velocityY[ix+1, iy] - velocityY[ix-1, iy]) -
                            (velocityX[ix, iy+1] - velocityX[ix, iy-1])) * vortScale
            end
        end
        obs_vort[] = vort
    end
end

function Create_Force_Plot()
    fig = Figure(size = (700, 350))
    ax  = Axis(fig[1, 1],
               title  = "Aerodynamic Coefficients",
               xlabel = "Time [s]",
               ylabel = "Coefficient [-]")

    obs_time = Observable(Float64[])
    obs_cd   = Observable(Float64[])
    obs_cl   = Observable(Float64[])

    lines!(ax, obs_time, obs_cd, color = :steelblue,  label = "CD")
    lines!(ax, obs_time, obs_cl, color = :orangered,  label = "CL")
    axislegend(ax, position = :rt)

    return fig, ax, obs_time, obs_cd, obs_cl
end

function Update_Force_Plot!(ax, obs_time, obs_cd, obs_cl, t, cd, cl)
    push!(obs_time[], t)
    push!(obs_cd[],   cd)
    push!(obs_cl[],   cl)
    notify(obs_time)
    notify(obs_cd)
    notify(obs_cl)
    autolimits!(ax)
end

end#Plotter
