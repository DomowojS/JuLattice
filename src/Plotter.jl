module Plotter
using GLMakie
export Create_Plot, Update_Plot!

function Create_Plot(Nx::Int, Ny::Int, plotU::Bool, plotV::Bool, plotVorticity::Bool,
                     rangeU::Tuple{Float64,Float64},
                     rangeV::Tuple{Float64,Float64},
                     rangeVort::Tuple{Float64,Float64})
    nplots = Int(plotU) + Int(plotV) + Int(plotVorticity)
    fig = Figure(size = (900, 280 * max(nplots, 1)))

    step_text = Observable("Step: 0  |  t = 0.00 s")
    obs_u    = nothing
    obs_v    = nothing
    obs_vort = nothing
    row = 1

    if plotU
        ax = Axis(fig[row, 1], title = "U  [m/s]", aspect = DataAspect())
        obs_u = Observable(zeros(Nx, Ny))
        hm = heatmap!(ax, obs_u, colormap = :inferno, colorrange = rangeU, nan_color = :dimgray)
        Colorbar(fig[row+1, 1], hm, vertical = false)
        text!(ax, 2, 2, text = step_text, color = :dimgray, fontsize = 11)
        row += 2
    end

    if plotV
        ax = Axis(fig[row, 1], title = "V  [m/s]", aspect = DataAspect())
        obs_v = Observable(zeros(Nx, Ny))
        hm = heatmap!(ax, obs_v, colormap = :inferno, colorrange = rangeV, nan_color = :dimgray)
        Colorbar(fig[row+1, 1], hm, vertical = false)
        row += 2
    end

    if plotVorticity
        ax = Axis(fig[row, 1], title = "Vorticity  [1/s]", aspect = DataAspect())
        obs_vort = Observable(zeros(Nx, Ny))
        hm = heatmap!(ax, obs_vort, colormap = :curl, colorrange = rangeVort, nan_color = :dimgray)
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

end#Plotter
