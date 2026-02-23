module Plotter
using GLMakie
export Create_Plot, Update_Plot!, Create_Force_Plot, Update_Force_Plot!

function Create_Plot(Nx::Int, Ny::Int,
                     NxFine::Int, NyFine::Int,
                     deltaX::Float64, deltaXFine::Float64,
                     originXFine::Float64, originYFine::Float64,
                     plotU::Bool, plotV::Bool, plotVorticity::Bool,
                     rangeU::Tuple{Float64,Float64},
                     rangeV::Tuple{Float64,Float64},
                     rangeVort::Tuple{Float64,Float64})
    nplots = Int(plotU) + Int(plotV) + Int(plotVorticity)
    fig = Figure(size = (900, 280 * max(nplots, 1)))

    # Coarse grid coordinates (ghost nodes at index 1 sit at -deltaX)
    xs = range(-deltaX, step = deltaX, length = Nx)
    ys = range(-deltaX, step = deltaX, length = Ny)

    # Fine grid coordinates — interior only (ixF=2:NxFine-1, skip ghosts)
    xs_fine = range(originXFine, step = deltaXFine, length = NxFine - 2)
    ys_fine = range(originYFine, step = deltaXFine, length = NyFine - 2)
    nxF = NxFine - 2
    nyF = NyFine - 2

    text_x = 0.02 * (Nx - 2) * deltaX
    text_y = 0.93 * (Ny - 2) * deltaX

    step_text     = Observable("Step: 0  |  t = 0.00 s")
    obs_u         = nothing
    obs_v         = nothing
    obs_vort      = nothing
    obs_u_fine    = nothing
    obs_v_fine    = nothing
    obs_vort_fine = nothing
    row = 1

    if plotU
        ax = Axis(fig[row, 1], title = "U  [m/s]", aspect = DataAspect(),
                  xlabel = "x [m]", ylabel = "y [m]")
        obs_u = Observable(zeros(Nx, Ny))
        heatmap!(ax, xs, ys, obs_u, colormap = :inferno, colorrange = rangeU, nan_color = :dimgray)
        obs_u_fine = Observable(fill(NaN, nxF, nyF))
        hm = heatmap!(ax, xs_fine, ys_fine, obs_u_fine, colormap = :inferno, colorrange = rangeU, nan_color = :dimgray)
        Colorbar(fig[row+1, 1], hm, vertical = false)
        text!(ax, text_x, text_y, text = step_text, color = :dimgray, fontsize = 11)
        row += 2
    end

    if plotV
        ax = Axis(fig[row, 1], title = "V  [m/s]", aspect = DataAspect(),
                  xlabel = "x [m]", ylabel = "y [m]")
        obs_v = Observable(zeros(Nx, Ny))
        heatmap!(ax, xs, ys, obs_v, colormap = :inferno, colorrange = rangeV, nan_color = :dimgray)
        obs_v_fine = Observable(fill(NaN, nxF, nyF))
        hm = heatmap!(ax, xs_fine, ys_fine, obs_v_fine, colormap = :inferno, colorrange = rangeV, nan_color = :dimgray)
        Colorbar(fig[row+1, 1], hm, vertical = false)
        row += 2
    end

    if plotVorticity
        ax = Axis(fig[row, 1], title = "Vorticity  [1/s]", aspect = DataAspect(),
                  xlabel = "x [m]", ylabel = "y [m]")
        obs_vort = Observable(zeros(Nx, Ny))
        heatmap!(ax, xs, ys, obs_vort, colormap = :curl, colorrange = rangeVort, nan_color = :dimgray)
        obs_vort_fine = Observable(fill(NaN, nxF, nyF))
        hm = heatmap!(ax, xs_fine, ys_fine, obs_vort_fine, colormap = :curl, colorrange = rangeVort, nan_color = :dimgray)
        Colorbar(fig[row+1, 1], hm, vertical = false)
        row += 2
    end

    return fig, obs_u, obs_v, obs_vort, obs_u_fine, obs_v_fine, obs_vort_fine, step_text
end

function Update_Plot!(obs_u, obs_v, obs_vort,
                      obs_u_fine, obs_v_fine, obs_vort_fine,
                      step_text,
                      velocityX, velocityY,
                      velocityXFine, velocityYFine,
                      i, deltaT, deltaX,
                      plotU::Bool, plotV::Bool, plotVorticity::Bool,
                      isFluid::BitMatrix,
                      isFluidFine::BitMatrix, isObjectFine::BitMatrix)

    velScale      = deltaX / deltaT
    vortScale     = 1.0 / (2.0 * deltaT)   # coarse: 1/(2*deltaT)
    vortScaleFine = 1.0 / deltaT            # fine:   1/(2*deltaTFine) = 1/(2*(deltaT/2)) = 1/deltaT

    step_text[] = "Step: $i  |  t = $(round(i*deltaT, digits=2)) s"

    NxFine, NyFine = size(velocityXFine)
    obj_mask_fine = @view isObjectFine[2:NxFine-1, 2:NyFine-1]

    if plotU
        if obs_u !== nothing
            obs_u[] = velocityX .* velScale
        end
        if obs_u_fine !== nothing
            data = velocityXFine[2:NxFine-1, 2:NyFine-1] .* velScale
            data[obj_mask_fine] .= NaN
            obs_u_fine[] = data
        end
    end

    if plotV
        if obs_v !== nothing
            obs_v[] = velocityY .* velScale
        end
        if obs_v_fine !== nothing
            data = velocityYFine[2:NxFine-1, 2:NyFine-1] .* velScale
            data[obj_mask_fine] .= NaN
            obs_v_fine[] = data
        end
    end

    if plotVorticity
        if obs_vort !== nothing
            vort = fill(NaN, size(velocityX))
            for I in CartesianIndices(isFluid)
                if isFluid[I]
                    ix, iy = I.I
                    vort[I] = ((velocityY[ix+1,iy] - velocityY[ix-1,iy]) -
                                (velocityX[ix,iy+1] - velocityX[ix,iy-1])) * vortScale
                end
            end
            obs_vort[] = vort
        end
        if obs_vort_fine !== nothing
            vort_fine = fill(NaN, NxFine-2, NyFine-2)
            # skip outermost row/col (ixF=2 and NxFine-1): ghost neighbors needed for central diff
            for iyF in 3:NyFine-2, ixF in 3:NxFine-2
                if isFluidFine[ixF, iyF] && !isObjectFine[ixF, iyF]
                    vort_fine[ixF-1, iyF-1] = ((velocityYFine[ixF+1,iyF] - velocityYFine[ixF-1,iyF]) -
                                                  (velocityXFine[ixF,iyF+1] - velocityXFine[ixF,iyF-1])) * vortScaleFine
                end
            end
            obs_vort_fine[] = vort_fine
        end
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
