module JuLattice
    ############################
    ## Main file for JuLattice #
    ############################
    include("src/Kernel.jl")
    include("src/GridSetup.jl")
    include("src/Plotter.jl")
    include("src/IO.jl")

    using GLMakie
    using .Kernel, .GridSetup, .Plotter, .IO

    function run()
        ## User Settings
        # Domain Settings
        lengthX = 8.0             # m
        lengthY = 6.0             # m

        # Object reference length (for reynoldsNumber; object itself added later)
        d = 0.5                   # m
        angleDeg = 30.0           # degrees
        positionX = 3.0
        positionY = lengthY/2
        # Fluid Properties
        reynoldsNumber  = 300
        machNumber      = 0.1        # Ma = U / c_s  (keep < 0.1 for incompressible)
        viscosity       = 0.0001      # m^2/s

        # Simulation Settings
        simulationTime = 3600.0   # s
        deltaX = 0.01              # m per lattice unit

        # Plot Requests
        plotU         = true
        plotV         = true
        plotVorticity = true

        plotUMin    =  -0.04      # m/s
        plotUMax    =  0.1     # m/s
        plotVMin    = -0.1     # m/s
        plotVMax    =  0.08    # m/s
        plotVortMin = -1.1      # 1/s
        plotVortMax =  1.1      # 1/s

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

        ##-------- Grid Setup --------##
        Nx = ceil(Int, lengthX / deltaX) + 2   # +2 for ghost ring
        Ny = ceil(Int, lengthY / deltaX) + 2

        ##-------- Node Classification --------##
        (; isInlet, isOutlet, isWall, isFluid, isObject, isSolid,
           fluidNodes, solidNodes, objectNodes, boundaryNodesAndDistances) =
            _classify_nodes(Nx, Ny, deltaX, positionX, positionY, d, angleDeg)

        ##-------- MRT Setup --------##
        omegaBGK      = 1.0 / (3.0 * latticeViscosity + 0.5)
        omegaAcoustic = 1.0

        Log_Discretization_Settings(deltaX, deltaT, omegaBGK, reynoldsNumber)

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
        fig, obs_u, obs_v, obs_vort, step_text = Create_Plot(Nx, Ny, plotU, plotV, plotVorticity,
                                                              rangeU, rangeV, rangeVort, deltaX)
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

        for i in 1:nSteps

            ##-- Collision + Streaming --##
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

            ##-- Boundary Conditions --##
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

            # Object bounce-back (Bouzidi) + momentum exchange
            forceX, forceY = _apply_bouzidi_bc!(boundaryNodesAndDistances,
                f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm)

            ##-- Logging & Plotting --##
            if (i % 100 == 0) || (i == nSteps)
                nups = length(fluidNodes) * i / (time() - t_start)
                Log_Simulation_Runtime(i, nSteps, nups)
                println("  CD = $(round(forceX*coeff_denom, digits=4))  CL = $(round(forceY*coeff_denom, digits=4))")
            end

            if (i % 10 == 0) || (i == nSteps)
                Update_Plot!(obs_u, obs_v, obs_vort, step_text,
                             velocityX, velocityY,
                             i, deltaT, deltaX,
                             plotU, plotV, plotVorticity,
                             isFluid, isObject)
                Update_Force_Plot!(force_ax, obs_time, obs_cd, obs_cl,
                                   i * deltaT,
                                   forceX * coeff_denom,
                                   forceY * coeff_denom)
                yield()
            end

        end#loop

        Save_Forces!(obs_time[], obs_cd[], obs_cl[])
        Log_Simulation_Tail()
    end#run

end#JuLattice
