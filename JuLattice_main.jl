module JuLattice
    ############################
    ## Main file for JuLattice #
    ############################
    include("src/Plotter.jl")
    include("src/Logger.jl")
    include("src/Equilibrium.jl")

    using GLMakie
    using .Plotter, .Logger, .Equilibrium

    function _collide_and_stream!(
            fluidNodes,
            f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm,
            f00S, fp0S, fm0S, f0pS, f0mS, fppS, fpmS, fmpS, fmmS,
            densityGrid, velocityX, velocityY,
            omegaBGK, omegaAcoustic)
        @inbounds Threads.@threads for k in eachindex(fluidNodes)
            idx = fluidNodes[k]
            ix, iy = Tuple(idx)
            # Macroscopics
            rho = f00[ix,iy] + fp0[ix,iy] + fm0[ix,iy] + f0p[ix,iy] + f0m[ix,iy] + fpp[ix,iy] + fpm[ix,iy] + fmp[ix,iy] + fmm[ix,iy]
            u = (-fmm[ix,iy] + fpp[ix,iy] - fmp[ix,iy] + fpm[ix,iy] - fm0[ix,iy] + fp0[ix,iy]) / rho
            v = (-fmm[ix,iy] + fpp[ix,iy] + fmp[ix,iy] - fpm[ix,iy] + f0p[ix,iy] - f0m[ix,iy]) / rho
            # Velocity Products
            uu = u * u
            vv = v * v

            densityGrid[ix,iy] = rho
            velocityX[ix,iy]   = u
            velocityY[ix,iy]   = v

            # f -> m
            m20 = (fmm[ix,iy] + fpp[ix,iy] + fmp[ix,iy] + fpm[ix,iy] + fm0[ix,iy] + fp0[ix,iy])
            m02 = (fmm[ix,iy] + fpp[ix,iy] + fmp[ix,iy] + fpm[ix,iy] + f0p[ix,iy] + f0m[ix,iy])
            m11 = (fmm[ix,iy] + fpp[ix,iy] - fmp[ix,iy] - fpm[ix,iy])
            m21 = (-fmm[ix,iy] + fpp[ix,iy] + fmp[ix,iy] - fpm[ix,iy])
            m12 = (-fmm[ix,iy] + fpp[ix,iy] - fmp[ix,iy] + fpm[ix,iy])
            m22 = (fmm[ix,iy] + fpp[ix,iy] + fmp[ix,iy] + fpm[ix,iy])

            mP  = m20 + m02
            mxx = m20 - m02

            # Relaxation
            m11 += omegaBGK * (rho*u*v - m11)
            mxx += omegaBGK * (rho*(uu - vv) - mxx)
            mP  += omegaAcoustic * (rho*(2.0/3.0 + uu + vv) - mP)

            m21 += 1.0 * (rho*(1.0/3.0 + uu)*v - m21)
            m12 += 1.0 * (rho*(1.0/3.0 + vv)*u - m12)
            m22 += 1.0 * (rho*(1.0/3.0 + uu)*(1.0/3.0 + vv) - m22)

            # Back-compute m20/m02 from relaxed mP and mxx
            m20 = 0.5 * (mP + mxx)
            m02 = 0.5 * (mP - mxx)

            # m -> f* and push-stream to destination
            fmmS[ix-1,iy-1] = 0.25*( m11 - m12 - m21 + m22 )
            f0mS[ix,  iy-1] = 0.5 *( -v*rho + m02 + m21 - m22 )
            fpmS[ix+1,iy-1] = 0.25*( -m11 + m12 - m21 + m22 )

            fm0S[ix-1,iy]   = 0.5 *( -u*rho + m20 + m12 - m22 )
            f00S[ix,  iy]   = rho - m02 - m20 + m22
            fp0S[ix+1,iy]   = 0.5 *(  u*rho + m20 - m12 - m22 )

            fmpS[ix-1,iy+1] = 0.25*( -m11 - m12 + m21 + m22 )
            f0pS[ix,  iy+1] = 0.5 *(  v*rho + m02 - m21 - m22 )
            fppS[ix+1,iy+1] = 0.25*(  m11 + m12 + m21 + m22 )
        end
    end

    function add_rectangle!(isObject, Nx, Ny, deltaX;
                            centerX, centerY, d, angleDeg)
        half_len = 1.5 * d   # half of 3d  (long axis)
        half_hgt = 0.5 * d   # half of d   (short axis)

        alpha = deg2rad(angleDeg)
        c = cos(alpha)
        s = sin(alpha)

        for ix in 1:Nx, iy in 1:Ny
            # Physical coordinates: ix=2,iy=2 → x=0,y=0
            x = (ix - 2) * deltaX
            y = (iy - 2) * deltaX

            # Translate to rectangle-local origin
            dx = x - centerX
            dy = y - centerY

            # Transform point into rectangle's local frame (inverse of clockwise rotation)
            lx = c * dx - s * dy
            ly = s * dx + c * dy

            if abs(lx) <= half_len && abs(ly) <= half_hgt
                isObject[ix, iy] = true
            end
        end
    end

    function _ray_rect_q(lx0, ly0, dlx, dly, half_len, half_hgt)
        q = Inf
        # Check lx = ±half_len sides
        if abs(dlx) > 1e-14
            for wall in (-half_len, half_len)
                t = (wall - lx0) / dlx
                if 0.0 < t <= 1.0 + 1e-10
                    ly_t = ly0 + t * dly
                    if abs(ly_t) <= half_hgt + 1e-10
                        q = min(q, t)
                    end
                end
            end
        end
        # Check ly = ±half_hgt sides
        if abs(dly) > 1e-14
            for wall in (-half_hgt, half_hgt)
                t = (wall - ly0) / dly
                if 0.0 < t <= 1.0 + 1e-10
                    lx_t = lx0 + t * dlx
                    if abs(lx_t) <= half_len + 1e-10
                        q = min(q, t)
                    end
                end
            end
        end
        return q
    end

    function find_object_boundary_nodes(isObject, fluidNodes,
                                        deltaX, center_x, center_y,
                                        half_len, half_hgt, cos_a, sin_a)
        dirs = ((1,0), (-1,0), (0,1), (0,-1), (1,1), (1,-1), (-1,1), (-1,-1))

        boundaryNodesAndDistances = Tuple{Int,Int,Int,Int,Float64}[]

        @inbounds for idx in fluidNodes
            ix, iy = Tuple(idx)
            # Fluid node in local frame
            dx0 = (ix - 2) * deltaX - center_x
            dy0 = (iy - 2) * deltaX - center_y
            lx0 =  cos_a * dx0 - sin_a * dy0
            ly0 =  sin_a * dx0 + cos_a * dy0

            for (cx, cy) in dirs
                if isObject[ix + cx, iy + cy]
                    # Object node in local frame
                    dx1 = (ix + cx - 2) * deltaX - center_x
                    dy1 = (iy + cy - 2) * deltaX - center_y
                    lx1 =  cos_a * dx1 - sin_a * dy1
                    ly1 =  sin_a * dx1 + cos_a * dy1
                    
                    # Find exact distance q to rectangle edge along ray from (lx0,ly0) to (lx1,ly1)
                    q = _ray_rect_q(lx0, ly0, lx1 - lx0, ly1 - ly0, half_len, half_hgt)
                    push!(boundaryNodesAndDistances, (ix, iy, cx, cy, q))
                end
            end
        end

        return boundaryNodesAndDistances
    end

    function run()
        ## User Settings
        # Domain Settings
        lengthX = 8.0             # m
        lengthY = 2.0             # m

        # Object reference length (for reynoldsNumber; object itself added later)
        d = 0.5                   # m
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
        plotVMin    = -0.05     # m/s
        plotVMax    =  0.03    # m/s
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

        # Node identifiers (Boolean masks, Nx x Ny)
        isInlet  = falses(Nx, Ny);  isInlet[1,    :]         .= true
        isOutlet = falses(Nx, Ny);  isOutlet[Nx,  :]         .= true
        isWall   = falses(Nx, Ny);  isWall[2:Nx-1, [1, Ny]]  .= true
        isFluid  = falses(Nx, Ny);  isFluid[2:Nx-1, 2:Ny-1]  .= true
        isObject = falses(Nx, Ny);  add_rectangle!(isObject, Nx, Ny, deltaX;
                                                   centerX=positionX, centerY=positionY, d=d, angleDeg=30.0)
        isFluid .&= .!isObject  # cut object nodes out of fluid
        isSolid  = isInlet .| isOutlet .| isWall .| isObject

        # Precomputed node index lists — rebuild after any mask change (e.g. adding a solid object)
        fluidNodes = findall(isFluid)
        solidNodes = findall(.!isFluid .& .!isInlet .& .!isOutlet .& .!isWall)
        objectNodes = findall(isObject)
        
        boundaryNodesAndDistances = find_object_boundary_nodes(isObject, fluidNodes,
                                                        deltaX, positionX, positionY,
                                                        1.5*d, 0.5*d, cosd(30), sind(30))

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
        inlet_add_fp0 = (2.0 / (9.0  * latticeSpeedOfSound^2)) * latticeInflowVelocity
        inlet_add_fdiagonal = (2.0 / (36.0 * latticeSpeedOfSound^2)) * latticeInflowVelocity

        ##-------- Plotting Setup --------##
        rangeU    = (plotUMin,    plotUMax)
        rangeV    = (plotVMin,    plotVMax)
        rangeVort = (plotVortMin, plotVortMax)
        fig, obs_u, obs_v, obs_vort, step_text = Create_Plot(Nx, Ny, plotU, plotV, plotVorticity,
                                                              rangeU, rangeV, rangeVort)
        screen = GLMakie.Screen()
        GLMakie.display(screen, fig)

        ##-------- Main Loop --------##
        Log_Simulation_Start()
        t_start = time()

        for i in 1:nSteps

            ##-- 1. Collision + Streaming --##
            # Clear post-collision buffers
            f00S .= 0.0; fp0S .= 0.0; fm0S .= 0.0
            f0pS .= 0.0; f0mS .= 0.0
            fppS .= 0.0; fpmS .= 0.0; fmpS .= 0.0; fmmS .= 0.0

            _collide_and_stream!(
                fluidNodes,
                f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm,
                f00S, fp0S, fm0S, f0pS, f0mS, fppS, fpmS, fmpS, fmmS,
                densityGrid, velocityX, velocityY,
                omegaBGK, omegaAcoustic)

            # Swap f <-> fS
            f00, f00S = f00S, f00
            fp0, fp0S = fp0S, fp0
            fm0, fm0S = fm0S, fm0
            f0p, f0pS = f0pS, f0p
            f0m, f0mS = f0mS, f0m
            fpp, fppS = fppS, fpp
            fpm, fpmS = fpmS, fpm
            fmp, fmpS = fmpS, fmp
            fmm, fmmS = fmmS, fmm

            ##-- 3. Boundary Conditions --##
            # Inlet:bounce-back at first fluid column (ix=2), reading from inlet ghost (ix=1)
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
            # Bottom wall (ghost at iy=1): set cy=+1 populations at first fluid row iy=2
            @inbounds @views begin
                f0p[2:Nx-1, 2] .= f0p_eq_wall
                fpp[2:Nx-1, 2] .= fpp_eq_wall
                fmp[2:Nx-1, 2] .= fmp_eq_wall
                # Top wall (ghost at iy=Ny): set cy=-1 populations at last fluid row iy=Ny-1
                f0m[2:Nx-1, Ny-1] .= f0m_eq_wall
                fpm[2:Nx-1, Ny-1] .= fpm_eq_wall
                fmm[2:Nx-1, Ny-1] .= fmm_eq_wall
            end

            # Object bounce-back (Bouzidi)
            @inbounds for (ix, iy, cx, cy, q) in boundaryNodesAndDistances
                q2 = 2.0 * q
                if q < 0.5
                    # f_ᾱ[ix,iy] = 2q*f_α[ix+cx,iy+cy] + (1-2q)*f_α[ix,iy]
                    if     cx ==  1 && cy ==  0;  fm0[ix,iy] = q2*fp0[ix+1,iy  ] + (1.0-q2)*fp0[ix,  iy  ]
                    elseif cx == -1 && cy ==  0;  fp0[ix,iy] = q2*fm0[ix-1,iy  ] + (1.0-q2)*fm0[ix,  iy  ]
                    elseif cx ==  0 && cy ==  1;  f0m[ix,iy] = q2*f0p[ix,  iy+1] + (1.0-q2)*f0p[ix,  iy  ]
                    elseif cx ==  0 && cy == -1;  f0p[ix,iy] = q2*f0m[ix,  iy-1] + (1.0-q2)*f0m[ix,  iy  ]
                    elseif cx ==  1 && cy ==  1;  fmm[ix,iy] = q2*fpp[ix+1,iy+1] + (1.0-q2)*fpp[ix,  iy  ]
                    elseif cx ==  1 && cy == -1;  fmp[ix,iy] = q2*fpm[ix+1,iy-1] + (1.0-q2)*fpm[ix,  iy  ]
                    elseif cx == -1 && cy ==  1;  fpm[ix,iy] = q2*fmp[ix-1,iy+1] + (1.0-q2)*fmp[ix,  iy  ]
                    elseif cx == -1 && cy == -1;  fpp[ix,iy] = q2*fmm[ix-1,iy-1] + (1.0-q2)*fmm[ix,  iy  ]
                    end
                else
                    iq2 = 1.0 / q2
                    r   = (q2 - 1.0) * iq2   # (2q-1)/(2q)
                    # f_ᾱ[ix,iy] = (1/2q)*f_α[ix+cx,iy+cy] + ((2q-1)/2q)*f_ᾱ[ix-cx,iy-cy]
                    if     cx ==  1 && cy ==  0;  fm0[ix,iy] = iq2*fp0[ix+1,iy  ] + r*fm0[ix-1,iy  ]
                    elseif cx == -1 && cy ==  0;  fp0[ix,iy] = iq2*fm0[ix-1,iy  ] + r*fp0[ix+1,iy  ]
                    elseif cx ==  0 && cy ==  1;  f0m[ix,iy] = iq2*f0p[ix,  iy+1] + r*f0m[ix,  iy-1]
                    elseif cx ==  0 && cy == -1;  f0p[ix,iy] = iq2*f0m[ix,  iy-1] + r*f0p[ix,  iy+1]
                    elseif cx ==  1 && cy ==  1;  fmm[ix,iy] = iq2*fpp[ix+1,iy+1] + r*fmm[ix-1,iy-1]
                    elseif cx ==  1 && cy == -1;  fmp[ix,iy] = iq2*fpm[ix+1,iy-1] + r*fmp[ix-1,iy+1]
                    elseif cx == -1 && cy ==  1;  fpm[ix,iy] = iq2*fmp[ix-1,iy+1] + r*fpm[ix+1,iy-1]
                    elseif cx == -1 && cy == -1;  fpp[ix,iy] = iq2*fmm[ix-1,iy-1] + r*fpp[ix+1,iy+1]
                    end
                end
            end

            ##-- 4. Logging & Plotting --##
            if (i % 100 == 0) || (i == nSteps)
                nups = length(fluidNodes) * i / (time() - t_start)
                Log_Simulation_Runtime(i, nSteps, nups)
            end

            if (i % 10 == 0) || (i == nSteps)
                Update_Plot!(obs_u, obs_v, obs_vort, step_text,
                             velocityX, velocityY,
                             i, deltaT, deltaX,
                             plotU, plotV, plotVorticity,
                             isFluid, isObject)
                yield()
            end

        end#loop

        Log_Simulation_Tail()
    end#run

end#JuLattice