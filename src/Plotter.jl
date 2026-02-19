module Plotter
using GLMakie
export Create_Plot, Update_Plot!

function Create_Plot(Nx::Int, Ny::Int, plotU::Bool, plotV::Bool, plotVorticity::Bool,
                     rangeU::Tuple{Float64,Float64},
                     rangeV::Tuple{Float64,Float64},
                     rangeVort::Tuple{Float64,Float64})
    nplots = Int(plotU) + Int(plotV) + Int(plotVorticity)
    fig = Figure(size = (900 * max(nplots, 1), 420))

    step_text = Observable("Step: 0  |  t = 0.00 s")
    obs_u    = nothing
    obs_v    = nothing
    obs_vort = nothing
    col = 1

    if plotU
        ax = Axis(fig[1, col], title = "U  [m/s]", aspect = DataAspect())
        obs_u = Observable(zeros(Nx, Ny))
        hm = heatmap!(ax, obs_u, colormap = :inferno, colorrange = rangeU, nan_color = :black)
        Colorbar(fig[2, col], hm, vertical = false)
        text!(ax, 2, 2, text = step_text, color = :white, fontsize = 11)
        col += 1
    end

    if plotV
        ax = Axis(fig[1, col], title = "V  [m/s]", aspect = DataAspect())
        obs_v = Observable(zeros(Nx, Ny))
        hm = heatmap!(ax, obs_v, colormap = :inferno, colorrange = rangeV, nan_color = :black)
        Colorbar(fig[2, col], hm, vertical = false)
        col += 1
    end

    if plotVorticity
        ax = Axis(fig[1, col], title = "Vorticity  [1/s]", aspect = DataAspect())
        obs_vort = Observable(zeros(Nx, Ny))
        hm = heatmap!(ax, obs_vort, colormap = :curl, colorrange = rangeVort, nan_color = :black)
        Colorbar(fig[2, col], hm, vertical = false)
        col += 1
    end

    return fig, obs_u, obs_v, obs_vort, step_text
end

function Update_Plot!(obs_u, obs_v, obs_vort, step_text,
                      velocityX, velocityY,
                      i, deltaT, deltaX,
                      plotU::Bool, plotV::Bool, plotVorticity::Bool,
                      isFluid::BitMatrix)

    velScale  = deltaX / deltaT          # [lu/ts] -> [m/s]
    vortScale = 1.0 / (2.0 * deltaT)    # central-diff (no /2) -> [1/s]

    step_text[] = "Step: $i  |  t = $(round(i*deltaT, digits=2)) s"

    if plotU && obs_u !== nothing
        obs_u[] = velocityX .* velScale
    end

    if plotV && obs_v !== nothing
        obs_v[] = velocityY .* velScale
    end

    if plotVorticity && obs_vort !== nothing
        vort = zeros(size(velocityX))
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

end#Plotter
