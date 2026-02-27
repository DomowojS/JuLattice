module JuLattice
    ############################
    ## Main file for JuLattice #
    ############################
    include("src/Kernel.jl")
    include("src/GridSetup.jl")
    include("src/GridRefinement.jl")
    include("src/Plotter.jl")
    include("src/IO.jl")

    using GLMakie
    using .Kernel, .GridSetup, .GridRefinement, .Plotter, .IO

    function run()
        ## User Settings
        # Domain Settings
        lengthX = 8.0             # m
        lengthY = 4.0             # m

        # Fine Grid Settings
        lengthXFine       = 4.0     # m  (width of fine region)
        lengthYFine       = 2.4     # m  (height of fine region)
        positionFineGridX = 1.5     # m  (lower-left anchor; snapped to nearest coarse node below)
        positionFineGridY = 0.8    # m

        # Object reference length (for reynoldsNumber; object itself added later)
        d = 0.5                   # m
        angleDeg = 30.0           # degrees
        positionX = 3.0
        positionY = lengthY/2
        # Fluid Properties
        reynoldsNumber  = 300
        machNumber      = 0.1        # Ma = U / c_s
        viscosity       = 0.0001      # m^2/s

        # Simulation Settings
        simulationTime = 10000   # s
        deltaX = 0.02             # m per lattice unit
        positionFineGridX, positionFineGridY =
            _snap_fine_grid_position(positionFineGridX, positionFineGridY, deltaX)

        # Plot Requests
        plotU             = true
        plotV             = true
        plotVorticity     = true
        plotVmag          = true
        plotGridBoundary  = true

        plotUMin     =  -0.1      # m/s
        plotUMax     =  0.16      # m/s
        plotVMin     = -0.1       # m/s
        plotVMax     =  0.08      # m/s
        plotVortMin  = -1.1       # 1/s
        plotVortMax  =  1.1       # 1/s
        plotVmagMin  =  0.0       # m/s
        plotVmagMax  =  0.16      # m/s

        #### Run Simulation #####
        Log_Simulation_Header()

        ##-------- LBM Parameters --------##
        latticeSpeedOfSound   = 1.0 / sqrt(3)
        speedOfSound          = viscosity * reynoldsNumber / (machNumber * d)
        deltaT                = deltaX * latticeSpeedOfSound / speedOfSound
        latticeViscosity      = viscosity * deltaT / deltaX^2
        latticeInflowVelocity = machNumber * latticeSpeedOfSound
        latticeDensity        = 1.0

        nSteps = ceil(Int, simulationTime / deltaT)

        # Fine grid spacing (acoustic scaling: same sound speed)
        deltaXFine = deltaX / 2
        deltaTFine = deltaT / 2

        ##-------- Grid Setup --------##
        Nx = ceil(Int, lengthX / deltaX) + 2   # +2 for ghost ring
        Ny = ceil(Int, lengthY / deltaX) + 2

        # Fine grid dimensions.
        # Snap fine box to the nearest integer number of coarse cells so that
        # ghost nodes land exactly at xi ± deltaX/4 (required for C↔F coupling).
        nCoarseX    = round(Int, lengthXFine / deltaX)
        nCoarseY    = round(Int, lengthYFine / deltaX)
        lengthXFine = nCoarseX * deltaX
        lengthYFine = nCoarseY * deltaX
        # 2 fine interior nodes per coarse cell + 2 ghost nodes
        NxFine = 2 * nCoarseX + 2
        NyFine = 2 * nCoarseY + 2

        ##-------- Node Classification --------##
        (; isInlet, isOutlet, isWall, isFluid, isObject, isSolid,
           isOuterInterfaceNode, isInnerInterfaceNode, isFineInterior,
           ix_left, ix_right, iy_bottom, iy_top) =
            _classify_coarse_nodes(Nx, Ny, deltaX,
                                   positionFineGridX, positionFineGridY,
                                   lengthXFine, lengthYFine)

        (; isFluidFine, isObjectFine,
           isOuterInterfaceNodeFine, isInnerInterfaceNodeFine,
           originXFine, originYFine) =
            _classify_fine_nodes(NxFine, NyFine, deltaXFine,
                                 positionFineGridX, positionFineGridY,
                                 positionX, positionY, d, angleDeg)

        ##-------- Precomputed Node Index Lists --------##
        (; fluidNodes, solidNodes, objectNodes,
           outerInterfaceNodes, innerInterfaceNodes,
           fluidNodesFine, objectNodesFine,
           outerInterfaceNodesFine, innerInterfaceNodesFine,
           boundaryNodesAndDistances, boundaryNodesAndDistancesFine) =
            _compute_node_lists(
                isFluid, isObject, isInlet, isOutlet, isWall,
                isOuterInterfaceNode, isInnerInterfaceNode,
                isFluidFine, isObjectFine, isOuterInterfaceNodeFine, isInnerInterfaceNodeFine,
                deltaX, deltaXFine, positionX, positionY, d, angleDeg,
                originXFine, originYFine)

        ##-------- MRT Setup --------##
        omegaBGK      = 1.0 / (3.0 * latticeViscosity + 0.5)
        omegaAcoustic = 1.0
        ##-------- Fine Grid MRT Setup --------##
        omegaBGKFine = 1.0 / (2.0 / omegaBGK - 0.5)   # acoustic scaling: (τF-0.5)=2(τC-0.5)

        Log_Discretization_Settings(deltaX, deltaT, omegaBGK, reynoldsNumber, latticeInflowVelocity*deltaX/deltaT)

        ##-------- Array Allocation --------##
        # Current distributions
        f00 = zeros(Nx,Ny); fp0 = zeros(Nx,Ny); fm0 = zeros(Nx,Ny)
        f0p = zeros(Nx,Ny); f0m = zeros(Nx,Ny)
        fpp = zeros(Nx,Ny); fpm = zeros(Nx,Ny); fmp = zeros(Nx,Ny); fmm = zeros(Nx,Ny)

        # Post collision distributions
        f00S = zeros(Nx,Ny); fp0S = zeros(Nx,Ny); fm0S = zeros(Nx,Ny)
        f0pS = zeros(Nx,Ny); f0mS = zeros(Nx,Ny)
        fppS = zeros(Nx,Ny); fpmS = zeros(Nx,Ny); fmpS = zeros(Nx,Ny); fmmS = zeros(Nx,Ny)

        # Macroscopic fields
        densityGrid = zeros(Nx,Ny)
        velocityX   = zeros(Nx,Ny)
        velocityY   = zeros(Nx,Ny)

        ##-------- Fine Grid Array Allocation --------##
        # Current distributions (Fine)
        f00Fine = zeros(NxFine,NyFine); fp0Fine = zeros(NxFine,NyFine); fm0Fine = zeros(NxFine,NyFine)
        f0pFine = zeros(NxFine,NyFine); f0mFine = zeros(NxFine,NyFine)
        fppFine = zeros(NxFine,NyFine); fpmFine = zeros(NxFine,NyFine)
        fmpFine = zeros(NxFine,NyFine); fmmFine = zeros(NxFine,NyFine)

        # Post-collision distributions (Fine)
        f00SFine = zeros(NxFine,NyFine); fp0SFine = zeros(NxFine,NyFine); fm0SFine = zeros(NxFine,NyFine)
        f0pSFine = zeros(NxFine,NyFine); f0mSFine = zeros(NxFine,NyFine)
        fppSFine = zeros(NxFine,NyFine); fpmSFine = zeros(NxFine,NyFine)
        fmpSFine = zeros(NxFine,NyFine); fmmSFine = zeros(NxFine,NyFine)

        # Macroscopic fields (Fine)
        densityGridFine = zeros(NxFine,NyFine)
        velocityXFine   = zeros(NxFine,NyFine)
        velocityYFine   = zeros(NxFine,NyFine)

        ##-------- Initialise via Equilibrium --------##
        u0   = latticeInflowVelocity
        rho0 = latticeDensity

        f00 .= getEquilibrium(rho0, u0, 0.0,  0,  0)
        fp0 .= getEquilibrium(rho0, u0, 0.0,  1,  0)
        fm0 .= getEquilibrium(rho0, u0, 0.0, -1,  0)
        f0p .= getEquilibrium(rho0, u0, 0.0,  0,  1)
        f0m .= getEquilibrium(rho0, u0, 0.0,  0, -1)
        fpp .= getEquilibrium(rho0, u0, 0.0,  1,  1)
        fpm .= getEquilibrium(rho0, u0, 0.0,  1, -1)
        fmp .= getEquilibrium(rho0, u0, 0.0, -1,  1)
        fmm .= getEquilibrium(rho0, u0, 0.0, -1, -1)

        ##-------- Initialise Fine Grid via Equilibrium --------##
        f00Fine .= getEquilibrium(rho0, u0, 0.0,  0,  0)
        fp0Fine .= getEquilibrium(rho0, u0, 0.0,  1,  0)
        fm0Fine .= getEquilibrium(rho0, u0, 0.0, -1,  0)
        f0pFine .= getEquilibrium(rho0, u0, 0.0,  0,  1)
        f0mFine .= getEquilibrium(rho0, u0, 0.0,  0, -1)
        fppFine .= getEquilibrium(rho0, u0, 0.0,  1,  1)
        fpmFine .= getEquilibrium(rho0, u0, 0.0,  1, -1)
        fmpFine .= getEquilibrium(rho0, u0, 0.0, -1,  1)
        fmmFine .= getEquilibrium(rho0, u0, 0.0, -1, -1)

        # Wall equilibrium: at (latticeDensity, 0, 0) — no-slip
        f00_eq_wall = getEquilibrium(rho0, 0.0, 0.0,  0,  0)
        fp0_eq_wall = getEquilibrium(rho0, 0.0, 0.0,  1,  0)
        fm0_eq_wall = getEquilibrium(rho0, 0.0, 0.0, -1,  0)
        f0p_eq_wall = getEquilibrium(rho0, 0.0, 0.0,  0,  1)
        f0m_eq_wall = getEquilibrium(rho0, 0.0, 0.0,  0, -1)
        fpp_eq_wall = getEquilibrium(rho0, 0.0, 0.0,  1,  1)
        fpm_eq_wall = getEquilibrium(rho0, 0.0, 0.0,  1, -1)
        fmp_eq_wall = getEquilibrium(rho0, 0.0, 0.0, -1,  1)
        fmm_eq_wall = getEquilibrium(rho0, 0.0, 0.0, -1, -1)

        # Inflow momentum for fp0, fpp and fpm
        inlet_add_fp0       = (2.0 / (9.0  * latticeSpeedOfSound^2)) * latticeInflowVelocity
        inlet_add_fdiagonal = (2.0 / (36.0 * latticeSpeedOfSound^2)) * latticeInflowVelocity

        ##-------- Plotting Setup --------##
        rangeU    = (plotUMin,    plotUMax)
        rangeV    = (plotVMin,    plotVMax)
        rangeVort = (plotVortMin, plotVortMax)
        rangeVmag = (plotVmagMin, plotVmagMax)
        fig, obs_u, obs_v, obs_vort, obs_vmag, obs_u_fine, obs_v_fine, obs_vort_fine, obs_vmag_fine, step_text,
            xs_plot, ys_plot, xs_fine_plot, ys_fine_plot =
            Create_Plot(Nx, Ny, NxFine, NyFine, deltaX, deltaXFine, originXFine, originYFine,
                        plotU, plotV, plotVorticity, plotVmag, plotGridBoundary, rangeU, rangeV, rangeVort, rangeVmag)
        screen = GLMakie.Screen()
        GLMakie.display(screen, fig)

        force_fig, force_ax, obs_time, obs_cd, obs_cl = Create_Force_Plot()
        force_screen = GLMakie.Screen()
        GLMakie.display(force_screen, force_fig)
        # Projected frontal length: shadow cast by the object onto the y-axis (⊥ to flow)
        proj_frontal = 2.0 * (1.5*d * sind(angleDeg) + 0.5*d * cosd(angleDeg)) / deltaX
        coeff_denom  = 1.0 / (0.5 * latticeDensity * latticeInflowVelocity^2 * proj_frontal)

        ##-------- Main Loop --------##
        Log_Simulation_Start()
        t_start = time()
        next_plot_save_time = 200.0   # physical seconds

        for i in 1:nSteps

            ##-- Fine Sub-steps (2 per coarse step) --##
            forceX = 0.0
            forceY = 0.0
            for _ in 1:2
                f00SFine .= 0.0; fp0SFine .= 0.0; fm0SFine .= 0.0
                f0pSFine .= 0.0; f0mSFine .= 0.0
                fppSFine .= 0.0; fpmSFine .= 0.0; fmpSFine .= 0.0; fmmSFine .= 0.0

                _collide_and_stream!(
                    fluidNodesFine,
                    f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine,
                    f00SFine, fp0SFine, fm0SFine, f0pSFine, f0mSFine, fppSFine, fpmSFine, fmpSFine, fmmSFine,
                    densityGridFine, velocityXFine, velocityYFine,
                    omegaBGKFine, omegaAcoustic)

                f00Fine, f00SFine = f00SFine, f00Fine
                fp0Fine, fp0SFine = fp0SFine, fp0Fine
                fm0Fine, fm0SFine = fm0SFine, fm0Fine
                f0pFine, f0pSFine = f0pSFine, f0pFine
                f0mFine, f0mSFine = f0mSFine, f0mFine
                fppFine, fppSFine = fppSFine, fppFine
                fpmFine, fpmSFine = fpmSFine, fpmFine
                fmpFine, fmpSFine = fmpSFine, fmpFine
                fmmFine, fmmSFine = fmmSFine, fmmFine

                dfx, dfy = _apply_bouzidi_bc!(boundaryNodesAndDistancesFine,
                    f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)
                forceX += dfx
                forceY += dfy
            end
            forceX *= 0.5
            forceY *= 0.5

            ##-- Coarse Step --##
            f00S .= 0.0; fp0S .= 0.0; fm0S .= 0.0
            f0pS .= 0.0; f0mS .= 0.0
            fppS .= 0.0; fpmS .= 0.0; fmpS .= 0.0; fmmS .= 0.0

            _collide_and_stream!(
                fluidNodes,
                f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm,
                f00S, fp0S, fm0S, f0pS, f0mS, fppS, fpmS, fmpS, fmmS,
                densityGrid, velocityX, velocityY,
                omegaBGK, omegaAcoustic)

            f00, f00S = f00S, f00
            fp0, fp0S = fp0S, fp0
            fm0, fm0S = fm0S, fm0
            f0p, f0pS = f0pS, f0p
            f0m, f0mS = f0mS, f0m
            fpp, fppS = fppS, fpp
            fpm, fpmS = fpmS, fpm
            fmp, fmpS = fmpS, fmp
            fmm, fmmS = fmmS, fmm

            # Inlet: bounce-back at first fluid column (ix=2), reading from inlet ghost (ix=1)
            for iy in 2:Ny-1
                ix = 2
                fp0[ix, iy] = fm0[ix-1, iy]   + inlet_add_fp0
                fpp[ix, iy] = fmm[ix-1, iy-1] + inlet_add_fdiagonal
                fpm[ix, iy] = fmp[ix-1, iy+1] + inlet_add_fdiagonal
            end

            # Outlet (x=Nx): zero-gradient -- copy eastward populations from last fluid node
            @views begin
                fm0[Nx-1, :] .= fm0[Nx-2, :]
                fmp[Nx-1, :] .= fmp[Nx-2, :]
                fmm[Nx-1, :] .= fmm[Nx-2, :]
            end

            # Walls: equilibrium at (latticeDensity, 0, 0) injected into first fluid row
            @inbounds @views begin
                f0p[2:Nx-1, 2] .= f0p_eq_wall
                fpp[2:Nx-1, 2] .= fpp_eq_wall
                fmp[2:Nx-1, 2] .= fmp_eq_wall
                f0m[2:Nx-1, Ny-1] .= f0m_eq_wall
                fpm[2:Nx-1, Ny-1] .= fpm_eq_wall
                fmm[2:Nx-1, Ny-1] .= fmm_eq_wall
            end

            ##-- Synchronization --##
            _synchronize!(
                innerInterfaceNodes, outerInterfaceNodes,
                innerInterfaceNodesFine, outerInterfaceNodesFine,
                ix_left, iy_bottom,
                omegaBGK, omegaBGKFine,
                f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm,
                f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)

            ##-- Logging & Plotting --##
            if (i % 100 == 0) || (i == nSteps)
                nups = (length(fluidNodes) + 2 * length(fluidNodesFine)) * i / (time() - t_start)
                Log_Simulation_Runtime(i, nSteps, nups)
                println("  CD = $(round(forceX*coeff_denom, digits=4))  CL = $(round(forceY*coeff_denom, digits=4))")
            end

            if (i % 10 == 0) || (i == nSteps)
                Update_Plot!(obs_u, obs_v, obs_vort, obs_vmag,
                             obs_u_fine, obs_v_fine, obs_vort_fine, obs_vmag_fine,
                             step_text,
                             velocityX, velocityY, velocityXFine, velocityYFine,
                             i, deltaT, deltaX,
                             plotU, plotV, plotVorticity, plotVmag,
                             isFluid, isFluidFine, isObjectFine)
                Update_Force_Plot!(force_ax, obs_time, obs_cd, obs_cl,
                                   i * deltaT,
                                   forceX * coeff_denom,
                                   forceY * coeff_denom)
                yield()
                t_now = i * deltaT
                if t_now >= next_plot_save_time || i == nSteps
                    Save_Contour_Plot!(fig, t_now)
                    Save_Contour_Images!(obs_u, obs_vmag, obs_vort, obs_u_fine, obs_vmag_fine, obs_vort_fine,
                                        xs_plot, ys_plot, xs_fine_plot, ys_fine_plot,
                                        rangeU, rangeVmag, rangeVort, t_now)
                    next_plot_save_time = floor(t_now / 200.0) * 200.0 + 200.0
                end
            end

        end#loop

        Save_Forces!(obs_time[], obs_cd[], obs_cl[])
        Log_Simulation_Tail()
    end#run

end#JuLattice
