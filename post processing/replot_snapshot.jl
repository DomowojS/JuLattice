include(joinpath(@__DIR__, "../src/Plotter.jl"))
using .Plotter
using Serialization, GLMakie
using NativeFileDialog

function replot()
    ## -------- Settings  -------- ##
    # snapshot_file = "../simulation_data/"
    snapshot_file = pick_file(; filterlist="jls")
    isempty(snapshot_file) && error("no file selected!!!")
    # :mag, :vx, :vorticity
    plot_mode = :mag       
    # plot_mode = :vorticity
    # :component, :magnitude
    vorticity_mode = :magnitude

    colormap_mag    = :Blues
    colorrange_mag  = (0.0, 0.06) 

    ## -------- Plotting  -------- ##
    s = deserialize(snapshot_file)
    println("Loaded snapshot $(s.run_tag)")

    if plot_mode == :mag
        setup_mag_plot(s.gridlengthX, s.gridlengthY, s.gridlengthZ,
                        s.velocityMag, s.midY, s.midZ;
                        colormap=colormap_mag, colorrange=colorrange_mag)

    elseif plot_mode == :vx
        setup_vx_plot(s.gridlengthX, s.gridlengthY, s.gridlengthZ,
        s.velocityX, s.midY, s.midZ)

    elseif plot_mode == :vorticity
        setup_vorticity_plot(s.gridlengthX, s.gridlengthY, s.gridlengthZ,
        s.vortZ, s.vortY, s.midY, s.midZ; mode=vorticity_mode)
    end

end #function

