module JuLattice
    ############################
    ## Main file for JuLattice #
    ############################
    include("src/Plotter.jl")
    include("src/Logger.jl")
    include("src/Equilibrium.jl")

    using GLMakie
    using .Plotter, .Logger, .Equilibrium

    function run()
        ## User Settings
        # Domain Settings
        lengthX = 8.0             # m
        lengthY = 2.0             # m

        # Object reference length (for reynoldsNumber; object itself added later)
        d = 0.5                   # m

        # Fluid Properties
        reynoldsNumber          = 300
        machNumber  = 0.025        # Ma = U / c_s  (keep < 0.1 for incompressible)
        viscosity   = 0.001       # m^2/s

        # Simulation Settings
        simulationTime = 8000.0   # s
        deltaX = 0.05              # m per lattice unit

        # Plot Requests
        plotU         = true
        plotV         = true
        plotVorticity = true

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

        # Node identifiers (Boolean masks, Nx x Ny)
        isInlet  = falses(Nx, Ny);  isInlet[1,    :]         .= true
        isOutlet = falses(Nx, Ny);  isOutlet[Nx,  :]         .= true
        isWall   = falses(Nx, Ny);  isWall[2:Nx-1, [1, Ny]]  .= true
        isFluid  = falses(Nx, Ny);  isFluid[2:Nx-1, 2:Ny-1]  .= true

        ##-------- MRT Setup --------##
        omegaBGK      = 1.0 / (3.0 * latticeViscosity + 0.5)
        omegaAcoustic = 1.0

        Log_Discretization_Settings(deltaX, deltaT, omegaBGK, reynoldsNumber)
        # Distribution ordering: [f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm]
        # Moment ordering:       [m00, m10, m01, m11, mP,  mxx, m21, m12, m22]
        #   mP  = (m20+m02) - 2cs²·m00  (weighted trace,   relaxes with omegaAcoustic)
        #   mxx = m20 - m02             (stress diff,      relaxes with omegaBGK)
        #   m21 = m21_raw - cs²·m01    (weighted 3rd-order ghost)
        #   m12 = m12_raw - cs²·m10    (weighted 3rd-order ghost)
        #   m22 = m22_raw - cs²·(m20+m02) + cs⁴·m00  (weighted 4th-order ghost)
        M = Float64[
            1      1      1      1      1      1      1      1      1   ;  # m00
            0      1     -1      0      0      1      1     -1     -1   ;  # m10
            0      0      0      1     -1      1     -1      1     -1   ;  # m01
            0      0      0      0      0      1     -1     -1      1   ;  # m11
           -2/3    1/3    1/3    1/3    1/3    4/3    4/3    4/3    4/3 ;  # mP  = (m20+m02) - 2cs²·m00
            0      1      1     -1     -1      0      0      0      0   ;  # mxx = m20-m02
            0      0      0     -1/3   1/3    2/3   -2/3    2/3   -2/3  ;  # m21 = m21 - cs²·m01
            0     -1/3   1/3     0      0     2/3    2/3   -2/3   -2/3  ;  # m12 = m12 - cs²·m10
           1/9   -2/9   -2/9   -2/9   -2/9   4/9    4/9    4/9    4/9  ]  # m22 = m22 - cs²·(m20+m02) + cs⁴·m00

        M_inv = inv(M)

        ##-------- Array Allocation --------##
        # Current distributions
        f00 = zeros(Nx,Ny); fp0 = zeros(Nx,Ny); fm0 = zeros(Nx,Ny)
        f0p = zeros(Nx,Ny); f0m = zeros(Nx,Ny)
        fpp = zeros(Nx,Ny); fpm = zeros(Nx,Ny); fmp = zeros(Nx,Ny); fmm = zeros(Nx,Ny)

        # Post-collision distributions (ping-pong buffer for streaming)
        f00s = zeros(Nx,Ny); fp0s = zeros(Nx,Ny); fm0s = zeros(Nx,Ny)
        f0ps = zeros(Nx,Ny); f0ms = zeros(Nx,Ny)
        fpps = zeros(Nx,Ny); fpms = zeros(Nx,Ny); fmps = zeros(Nx,Ny); fmms = zeros(Nx,Ny)

        # Moments (overwritten each step, fluid nodes only)
        m00 = zeros(Nx,Ny); m10 = zeros(Nx,Ny); m01 = zeros(Nx,Ny)
        m11 = zeros(Nx,Ny); mP  = zeros(Nx,Ny); mxx = zeros(Nx,Ny)
        m21 = zeros(Nx,Ny); m12 = zeros(Nx,Ny); m22 = zeros(Nx,Ny)

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

        f00s .= f00; fp0s .= fp0; fm0s .= fm0
        f0ps .= f0p; f0ms .= f0m
        fpps .= fpp; fpms .= fpm; fmps .= fmp; fmms .= fmm

        densityGrid .= rho0
        velocityX   .= u0

        # Inlet equilibrium: at (latticeDensity, latticeInflowVelocity, 0)
        fp0_eq = getEquilibrium(rho0, latticeInflowVelocity, 0.0,  1,  0)
        fpp_eq = getEquilibrium(rho0, latticeInflowVelocity, 0.0,  1,  1)
        fpm_eq = getEquilibrium(rho0, latticeInflowVelocity, 0.0,  1, -1)

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

        ##-------- Plotting Setup --------##
        fig, obs_u, obs_v, obs_vort, step_text = Create_Plot(Nx, Ny, plotU, plotV, plotVorticity)
        screen = GLMakie.Screen()
        GLMakie.display(screen, fig)

        ##-------- Main Loop --------##
        Log_Simulation_Start()

        for i in 1:nSteps

            ##-- 1. Collision (fluid nodes only) --##
            @views begin
                # f -> m  (M applied row by row)
                m00[isFluid] .= f00[isFluid] .+ fp0[isFluid] .+ fm0[isFluid] .+ f0p[isFluid] .+ f0m[isFluid] .+ fpp[isFluid] .+ fpm[isFluid] .+ fmp[isFluid] .+ fmm[isFluid]
                m10[isFluid] .= fp0[isFluid] .- fm0[isFluid] .+ fpp[isFluid] .+ fpm[isFluid] .- fmp[isFluid] .- fmm[isFluid]
                m01[isFluid] .= f0p[isFluid] .- f0m[isFluid] .+ fpp[isFluid] .- fpm[isFluid] .+ fmp[isFluid] .- fmm[isFluid]
                m11[isFluid] .= fpp[isFluid] .- fpm[isFluid] .- fmp[isFluid] .+ fmm[isFluid]
                mP[isFluid]  .= fp0[isFluid] .+ fm0[isFluid] .+ f0p[isFluid] .+ f0m[isFluid] .+ 2.0.*(fpp[isFluid] .+ fpm[isFluid] .+ fmp[isFluid] .+ fmm[isFluid])
                mxx[isFluid] .= fp0[isFluid] .+ fm0[isFluid] .- f0p[isFluid] .- f0m[isFluid]
                m21[isFluid] .= fpp[isFluid] .- fpm[isFluid] .+ fmp[isFluid] .- fmm[isFluid]
                m12[isFluid] .= fpp[isFluid] .+ fpm[isFluid] .- fmp[isFluid] .- fmm[isFluid]
                m22[isFluid] .= fpp[isFluid] .+ fpm[isFluid] .+ fmp[isFluid] .+ fmm[isFluid]

                # Macroscopic variables (m00=ρ, m10=ρu, m01=ρv are conserved throughout relaxation)
                densityGrid[isFluid] .= m00[isFluid]
                velocityX[isFluid]   .= m10[isFluid] ./ m00[isFluid]
                velocityY[isFluid]   .= m01[isFluid] ./ m00[isFluid]

                # Relaxation  m_s = m - omega*(m - m_eq)
                # Equilibria expressed via conserved moments only: m00, m10, m01
                m11[isFluid] .-= omegaBGK      .* (m11[isFluid] .- m10[isFluid] .* m01[isFluid] ./ m00[isFluid])
                mP[isFluid]  .-= omegaAcoustic .* (mP[isFluid]  .- (m10[isFluid].^2 .+ m01[isFluid].^2) ./ m00[isFluid])
                mxx[isFluid] .-= omegaBGK      .* (mxx[isFluid] .- (m10[isFluid].^2 .- m01[isFluid].^2) ./ m00[isFluid])
                m21[isFluid] .-= 1.0           .* (m21[isFluid] .- m10[isFluid].^2 .* m01[isFluid] ./ m00[isFluid].^2)
                m12[isFluid] .-= 1.0           .* (m12[isFluid] .- m10[isFluid] .* m01[isFluid].^2 ./ m00[isFluid].^2)
                m22[isFluid] .-= 1.0           .* (m22[isFluid] .- m10[isFluid].^2 .* m01[isFluid].^2 ./ m00[isFluid].^3)

                # m -> fs  (M_inv applied row by row)
                f00s[isFluid] .= M_inv[1,1].*m00[isFluid] .+ M_inv[1,2].*m10[isFluid] .+ M_inv[1,3].*m01[isFluid] .+ M_inv[1,4].*m11[isFluid] .+ M_inv[1,5].*mP[isFluid] .+ M_inv[1,6].*mxx[isFluid] .+ M_inv[1,7].*m21[isFluid] .+ M_inv[1,8].*m12[isFluid] .+ M_inv[1,9].*m22[isFluid]
                fp0s[isFluid] .= M_inv[2,1].*m00[isFluid] .+ M_inv[2,2].*m10[isFluid] .+ M_inv[2,3].*m01[isFluid] .+ M_inv[2,4].*m11[isFluid] .+ M_inv[2,5].*mP[isFluid] .+ M_inv[2,6].*mxx[isFluid] .+ M_inv[2,7].*m21[isFluid] .+ M_inv[2,8].*m12[isFluid] .+ M_inv[2,9].*m22[isFluid]
                fm0s[isFluid] .= M_inv[3,1].*m00[isFluid] .+ M_inv[3,2].*m10[isFluid] .+ M_inv[3,3].*m01[isFluid] .+ M_inv[3,4].*m11[isFluid] .+ M_inv[3,5].*mP[isFluid] .+ M_inv[3,6].*mxx[isFluid] .+ M_inv[3,7].*m21[isFluid] .+ M_inv[3,8].*m12[isFluid] .+ M_inv[3,9].*m22[isFluid]
                f0ps[isFluid] .= M_inv[4,1].*m00[isFluid] .+ M_inv[4,2].*m10[isFluid] .+ M_inv[4,3].*m01[isFluid] .+ M_inv[4,4].*m11[isFluid] .+ M_inv[4,5].*mP[isFluid] .+ M_inv[4,6].*mxx[isFluid] .+ M_inv[4,7].*m21[isFluid] .+ M_inv[4,8].*m12[isFluid] .+ M_inv[4,9].*m22[isFluid]
                f0ms[isFluid] .= M_inv[5,1].*m00[isFluid] .+ M_inv[5,2].*m10[isFluid] .+ M_inv[5,3].*m01[isFluid] .+ M_inv[5,4].*m11[isFluid] .+ M_inv[5,5].*mP[isFluid] .+ M_inv[5,6].*mxx[isFluid] .+ M_inv[5,7].*m21[isFluid] .+ M_inv[5,8].*m12[isFluid] .+ M_inv[5,9].*m22[isFluid]
                fpps[isFluid] .= M_inv[6,1].*m00[isFluid] .+ M_inv[6,2].*m10[isFluid] .+ M_inv[6,3].*m01[isFluid] .+ M_inv[6,4].*m11[isFluid] .+ M_inv[6,5].*mP[isFluid] .+ M_inv[6,6].*mxx[isFluid] .+ M_inv[6,7].*m21[isFluid] .+ M_inv[6,8].*m12[isFluid] .+ M_inv[6,9].*m22[isFluid]
                fpms[isFluid] .= M_inv[7,1].*m00[isFluid] .+ M_inv[7,2].*m10[isFluid] .+ M_inv[7,3].*m01[isFluid] .+ M_inv[7,4].*m11[isFluid] .+ M_inv[7,5].*mP[isFluid] .+ M_inv[7,6].*mxx[isFluid] .+ M_inv[7,7].*m21[isFluid] .+ M_inv[7,8].*m12[isFluid] .+ M_inv[7,9].*m22[isFluid]
                fmps[isFluid] .= M_inv[8,1].*m00[isFluid] .+ M_inv[8,2].*m10[isFluid] .+ M_inv[8,3].*m01[isFluid] .+ M_inv[8,4].*m11[isFluid] .+ M_inv[8,5].*mP[isFluid] .+ M_inv[8,6].*mxx[isFluid] .+ M_inv[8,7].*m21[isFluid] .+ M_inv[8,8].*m12[isFluid] .+ M_inv[8,9].*m22[isFluid]
                fmms[isFluid] .= M_inv[9,1].*m00[isFluid] .+ M_inv[9,2].*m10[isFluid] .+ M_inv[9,3].*m01[isFluid] .+ M_inv[9,4].*m11[isFluid] .+ M_inv[9,5].*mP[isFluid] .+ M_inv[9,6].*mxx[isFluid] .+ M_inv[9,7].*m21[isFluid] .+ M_inv[9,8].*m12[isFluid] .+ M_inv[9,9].*m22[isFluid]
            end

            ##-- 2. Streaming (full ranges; ghost garbage overwritten by BCs) --##
            f00 .= f00s
            fp0[2:Nx,    :]    .= fp0s[1:Nx-1, :]        # cx=+1
            fm0[1:Nx-1,  :]    .= fm0s[2:Nx,   :]        # cx=-1
            f0p[:,   2:Ny]     .= f0ps[:,   1:Ny-1]      # cy=+1
            f0m[:,   1:Ny-1]   .= f0ms[:,   2:Ny  ]      # cy=-1
            fpp[2:Nx,   2:Ny]  .= fpps[1:Nx-1, 1:Ny-1]  # cx=+1,cy=+1
            fpm[2:Nx,   1:Ny-1].= fpms[1:Nx-1, 2:Ny  ]  # cx=+1,cy=-1
            fmp[1:Nx-1, 2:Ny]  .= fmps[2:Nx,   1:Ny-1]  # cx=-1,cy=+1
            fmm[1:Nx-1, 1:Ny-1].= fmms[2:Nx,   2:Ny  ]  # cx=-1,cy=-1

            ##-- 3. Boundary Conditions --##

            # Inlet (x=1): velocity anti-bounce-back -> inject into first fluid node x=2
            # Post-streaming: fm0[1,y]=fm0s[2,y], fmm[1,y]=fmms[2,y+1], fmp[1,y]=fmps[2,y-1]
            @views begin
                fp0[2, :]      .= .-fm0[1, :]      .+ 2*fp0_eq
                fpp[2, 2:Ny]   .= .-fmm[1, 1:Ny-1] .+ 2*fpp_eq
                fpm[2, 1:Ny-1] .= .-fmp[1, 2:Ny]   .+ 2*fpm_eq
            end

            # Outlet (x=Nx): zero-gradient -- copy eastward populations from last fluid node
            @views begin
                fp0[Nx, :] .= fp0[Nx-1, :]
                fpp[Nx, :] .= fpp[Nx-1, :]
                fpm[Nx, :] .= fpm[Nx-1, :]
            end

            # Walls (y=1 bottom, y=Ny top): equilibrium at (latticeDensity, 0, 0) -- no-slip
            for iy in (1, Ny)
                @views begin
                    f00[2:Nx-1, iy] .= f00_eq_wall
                    fp0[2:Nx-1, iy] .= fp0_eq_wall;  fm0[2:Nx-1, iy] .= fm0_eq_wall
                    f0p[2:Nx-1, iy] .= f0p_eq_wall;  f0m[2:Nx-1, iy] .= f0m_eq_wall
                    fpp[2:Nx-1, iy] .= fpp_eq_wall;  fpm[2:Nx-1, iy] .= fpm_eq_wall
                    fmp[2:Nx-1, iy] .= fmp_eq_wall;  fmm[2:Nx-1, iy] .= fmm_eq_wall
                end
            end

            ##-- 4. Logging & Plotting --##
            if (i % 100 == 0) || (i == nSteps)
                Log_Simulation_Runtime(i, nSteps)
            end

            if (i % 10 == 0) || (i == nSteps)
                Update_Plot!(obs_u, obs_v, obs_vort, step_text,
                             velocityX, velocityY,
                             i, deltaT,
                             plotU, plotV, plotVorticity,
                             isFluid)
                yield()
                sleep(0.01)
            end

        end#loop

        Log_Simulation_Tail()
    end#run

end#JuLattice