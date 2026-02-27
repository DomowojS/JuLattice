module Plotter
using GLMakie
using Printf
export Create_Plot, Update_Plot!, Create_Force_Plot, Update_Force_Plot!, Save_Contour_Plot!, Save_Contour_Images!

function Create_Plot(Nx::Int, Ny::Int,
                     NxFine::Int, NyFine::Int,
                     deltaX::Float64, deltaXFine::Float64,
                     originXFine::Float64, originYFine::Float64,
                     plotU::Bool, plotV::Bool, plotVorticity::Bool, plotVmag::Bool, plotGridBoundary::Bool,
                     rangeU::Tuple{Float64,Float64},
                     rangeV::Tuple{Float64,Float64},
                     rangeVort::Tuple{Float64,Float64},
                     rangeVmag::Tuple{Float64,Float64})
    nplots = Int(plotU) + Int(plotV) + Int(plotVorticity) + Int(plotVmag)
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

    # Fine box boundary in physical coords (coarse nodes at box edge)
    x_box_left   = originXFine - 0.5 * deltaXFine
    x_box_right  = x_box_left + (NxFine - 2) * deltaXFine
    y_box_bottom = originYFine - 0.5 * deltaXFine
    y_box_top    = y_box_bottom + (NyFine - 2) * deltaXFine
    box_xs = [x_box_left, x_box_right, x_box_right, x_box_left, x_box_left]
    box_ys = [y_box_bottom, y_box_bottom, y_box_top, y_box_top, y_box_bottom]

    step_text      = Observable("Step: 0  |  t = 0.00 s")
    obs_u          = nothing
    obs_v          = nothing
    obs_vort       = nothing
    obs_vmag       = nothing
    obs_u_fine     = nothing
    obs_v_fine     = nothing
    obs_vort_fine  = nothing
    obs_vmag_fine  = nothing
    row = 1

    if plotU
        ax = Axis(fig[row, 1], title = "U  [m/s]", aspect = DataAspect(),
                  xlabel = "x [m]", ylabel = "y [m]")
        obs_u = Observable(zeros(Nx, Ny))
        heatmap!(ax, xs, ys, obs_u, colormap = :inferno, colorrange = rangeU, nan_color = :dimgray)
        obs_u_fine = Observable(fill(NaN, nxF, nyF))
        hm = heatmap!(ax, xs_fine, ys_fine, obs_u_fine, colormap = :inferno, colorrange = rangeU, nan_color = :dimgray)
        plotGridBoundary && lines!(ax, box_xs, box_ys, color = :black, linestyle = :dot, linewidth = 1.5)
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
        plotGridBoundary && lines!(ax, box_xs, box_ys, color = :black, linestyle = :dot, linewidth = 1.5)
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
        plotGridBoundary && lines!(ax, box_xs, box_ys, color = :black, linestyle = :dot, linewidth = 1.5)
        Colorbar(fig[row+1, 1], hm, vertical = false)
        row += 2
    end

    if plotVmag
        ax = Axis(fig[row, 1], title = "|V|  [m/s]", aspect = DataAspect(),
                  xlabel = "x [m]", ylabel = "y [m]")
        obs_vmag = Observable(zeros(Nx, Ny))
        heatmap!(ax, xs, ys, obs_vmag, colormap = :viridis, colorrange = rangeVmag, nan_color = :dimgray)
        obs_vmag_fine = Observable(fill(NaN, nxF, nyF))
        hm = heatmap!(ax, xs_fine, ys_fine, obs_vmag_fine, colormap = :viridis, colorrange = rangeVmag, nan_color = :dimgray)
        plotGridBoundary && lines!(ax, box_xs, box_ys, color = :black, linestyle = :dot, linewidth = 1.5)
        Colorbar(fig[row+1, 1], hm, vertical = false)
        row += 2
    end

    return fig, obs_u, obs_v, obs_vort, obs_vmag, obs_u_fine, obs_v_fine, obs_vort_fine, obs_vmag_fine, step_text,
           xs, ys, xs_fine, ys_fine
end

function Update_Plot!(obs_u, obs_v, obs_vort, obs_vmag,
                      obs_u_fine, obs_v_fine, obs_vort_fine, obs_vmag_fine,
                      step_text,
                      velocityX, velocityY,
                      velocityXFine, velocityYFine,
                      i, deltaT, deltaX,
                      plotU::Bool, plotV::Bool, plotVorticity::Bool, plotVmag::Bool,
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

    if plotVmag
        if obs_vmag !== nothing
            obs_vmag[] = @. sqrt((velocityX * velScale)^2 + (velocityY * velScale)^2)
        end
        if obs_vmag_fine !== nothing
            ux = velocityXFine[2:NxFine-1, 2:NyFine-1] .* velScale
            uy = velocityYFine[2:NxFine-1, 2:NyFine-1] .* velScale
            data = @. sqrt(ux^2 + uy^2)
            data[obj_mask_fine] .= NaN
            obs_vmag_fine[] = data
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
    ylims!(ax, -5.0, 10.0)
end

function Save_Contour_Images!(obs_u, obs_vmag, obs_vort, obs_u_fine, obs_vmag_fine, obs_vort_fine,
                               xs, ys, xs_fine, ys_fine,
                               rangeU::Tuple{Float64,Float64},
                               rangeVmag::Tuple{Float64,Float64},
                               rangeVort::Tuple{Float64,Float64},
                               t::Float64; dir::String = "./output")
    isdir(dir) || mkpath(dir)

    function _bare_heatmap(xs_, ys_, data_, cmap, crange)
        px_h = 600
        phys_w = Float64(xs_[end] - xs_[1])
        phys_h = Float64(ys_[end] - ys_[1])
        px_w = round(Int, px_h * phys_w / phys_h)
        fig_ = Figure(size = (px_w, px_h), figure_padding = 0)
        ax_  = Axis(fig_[1,1];
                    xautolimitmargin = (0f0, 0f0),
                    yautolimitmargin = (0f0, 0f0))
        hidedecorations!(ax_)
        hidespines!(ax_)
        heatmap!(ax_, xs_, ys_, data_; colormap = cmap, colorrange = crange, nan_color = :dimgray)
        return fig_
    end

    if obs_u !== nothing && obs_u_fine !== nothing
        fig_u = _bare_heatmap(xs, ys, obs_u[], :inferno, rangeU)
        heatmap!(fig_u.content[1], xs_fine, ys_fine, obs_u_fine[];
                 colormap = :inferno, colorrange = rangeU, nan_color = :dimgray)
        save(joinpath(dir, @sprintf("contour_u_t%07.1f.png", t)), fig_u; px_per_unit = 4)
    end

    if obs_vmag !== nothing && obs_vmag_fine !== nothing
        fig_m = _bare_heatmap(xs, ys, obs_vmag[], :viridis, rangeVmag)
        heatmap!(fig_m.content[1], xs_fine, ys_fine, obs_vmag_fine[];
                 colormap = :viridis, colorrange = rangeVmag, nan_color = :dimgray)
        save(joinpath(dir, @sprintf("contour_vmag_t%07.1f.png", t)), fig_m; px_per_unit = 4)
    end

    if obs_vort !== nothing && obs_vort_fine !== nothing
        fig_v = _bare_heatmap(xs, ys, obs_vort[], :curl, rangeVort)
        heatmap!(fig_v.content[1], xs_fine, ys_fine, obs_vort_fine[];
                 colormap = :curl, colorrange = rangeVort, nan_color = :dimgray)
        save(joinpath(dir, @sprintf("contour_vort_t%07.1f.png", t)), fig_v; px_per_unit = 4)
    end
end

function Save_Contour_Plot!(fig, t::Float64; dir::String = "./output")
    isdir(dir) || mkpath(dir)
    fname = @sprintf("contour_t%07.1f.png", t)
    save(joinpath(dir, fname), fig; px_per_unit = 4)
end

end # module Plotter
