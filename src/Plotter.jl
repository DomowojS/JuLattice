module Plotter
using GLMakie
export Create_Plot, Update_Plot!

function Create_Plot(Nx::Int, Ny::Int, plotU::Bool, plotV::Bool, plotVorticity::Bool)
    nplots = Int(plotU) + Int(plotV) + Int(plotVorticity)
    fig = Figure(size = (900 * max(nplots, 1), 420))

    step_text = Observable("Step: 0  |  t = 0.00 s")
    obs_u    = nothing
    obs_v    = nothing
    obs_vort = nothing
    col = 1

    if plotU
        ax = Axis(fig[1, col], title = "U  (x-velocity)", aspect = DataAspect())
        obs_u = Observable(zeros(Nx, Ny))
        hm = heatmap!(ax, obs_u, colormap = :inferno, colorrange = (-0.08, 0.08), nan_color = :black)
        Colorbar(fig[2, col], hm, vertical = false)
        text!(ax, 2, 2, text = step_text, color = :white, fontsize = 11)
        col += 1
    end

    if plotV
        ax = Axis(fig[1, col], title = "V  (y-velocity)", aspect = DataAspect())
        obs_v = Observable(zeros(Nx, Ny))
        hm = heatmap!(ax, obs_v, colormap = :inferno, colorrange = (-0.04, 0.04), nan_color = :black)
        Colorbar(fig[2, col], hm, vertical = false)
        col += 1
    end

    if plotVorticity
        ax = Axis(fig[1, col], title = "Vorticity", aspect = DataAspect())
        obs_vort = Observable(zeros(Nx, Ny))
        hm = heatmap!(ax, obs_vort, colormap = :curl, colorrange = (-0.05, 0.05), nan_color = :black)
        Colorbar(fig[2, col], hm, vertical = false)
        col += 1
    end

    return fig, obs_u, obs_v, obs_vort, step_text
end

function Update_Plot!(obs_u, obs_v, obs_vort, step_text,
                      velocityX, velocityY,
                      i, deltaT,
                      plotU::Bool, plotV::Bool, plotVorticity::Bool,
                      isFluid::BitMatrix)

    step_text[] = "Step: $i  |  t = $(round(i*deltaT, digits=2)) s"

    if plotU && obs_u !== nothing
        obs_u[] = copy(velocityX)
    end

    if plotV && obs_v !== nothing
        obs_v[] = copy(velocityY)
    end

    if plotVorticity && obs_vort !== nothing
        vort = zeros(size(velocityX))
        # Central difference over fluid nodes (works for arbitrary fluid masks)
        for I in CartesianIndices(isFluid)
            if isFluid[I]
                ix, iy = I.I
                vort[I] = (velocityY[ix+1, iy] - velocityY[ix-1, iy]) -
                           (velocityX[ix, iy+1] - velocityX[ix, iy-1])
            end
        end
        obs_vort[] = vort
    end
end

end#Plotter
