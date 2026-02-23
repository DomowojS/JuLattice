module JuLattice
    ############################
    ## Main file for JuLattice #
    ############################
    include("src/Plotter.jl")
    include("src/Logger.jl")
    include("src/Equilibrium.jl")

    using GLMakie
    using .Plotter, .Logger, .Equilibrium

    @inline function _momentum_exchange(cx, cy, f_in, f_out)
        s = f_in + f_out
        return cx * s, cy * s
    end

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
        lengthY = 6.0             # m

        # Fine Grid Settings
        lengthXFine       = 4.0     # m  (width of fine region)
        lengthYFine       = 3.0     # m  (height of fine region)
        positionFineGridX = 1.5     # m  (lower-left anchor; should lie on a coarse node)
        positionFineGridY = 1.5     # m

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
        deltaX = 0.01             # m per lattice unit

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

        # Fine grid spacing (acoustic scaling: same sound speed)
        deltaXFine = deltaX / 2
        deltaTFine = deltaT / 2

        ##-------- Grid Setup --------##
        Nx = ceil(Int, lengthX / deltaX) + 2   # +2 for ghost ring
        Ny = ceil(Int, lengthY / deltaX) + 2

        # Fine grid dimensions
        # Fine node (ixF, iyF) sits at physical coordinates:
        #   x = positionFineGridX + (ixF - 1.5) * deltaXFine
        #   y = positionFineGridY + (iyF - 1.5) * deltaXFine
        # The offset of deltaXFine/2 = deltaX/4 from the anchor means fine nodes
        # are staggered: between coarse nodes at 0 and deltaX, fine nodes sit at
        # deltaX/4 and 3*deltaX/4 — never coinciding, symmetrically embedded.
        NxFine = ceil(Int, lengthXFine / deltaXFine) + 2
        NyFine = ceil(Int, lengthYFine / deltaXFine) + 2

        # Node identifiers (Boolean masks, Nx x Ny)
        isInlet  = falses(Nx, Ny);  isInlet[1,    :]         .= true
        isOutlet = falses(Nx, Ny);  isOutlet[Nx,  :]         .= true
        isWall   = falses(Nx, Ny);  isWall[2:Nx-1, [1, Ny]]  .= true
        isFluid  = falses(Nx, Ny);  isFluid[2:Nx-1, 2:Ny-1]  .= true
        isObject = falses(Nx, Ny);  add_rectangle!(isObject, Nx, Ny, deltaX;
                                                   centerX=positionX, centerY=positionY, d=d, angleDeg=angleDeg)
        isFluid .&= .!isObject  # cut object nodes out of fluid
        isSolid  = isInlet .| isOutlet .| isWall .| isObject

        ##-------- Fine Grid Interface Node Classification (Coarse Grid) --------##
        # Map fine grid physical boundaries onto coarse grid indices
        # (valid only when anchor lies on a coarse node, i.e. positionFineGridX/deltaX is integer)
        ix_fine_left   = 2 + round(Int, positionFineGridX / deltaX)
        ix_fine_right  = 2 + round(Int, (positionFineGridX + lengthXFine) / deltaX)
        iy_fine_bottom = 2 + round(Int, positionFineGridY / deltaX)
        iy_fine_top    = 2 + round(Int, (positionFineGridY + lengthYFine) / deltaX)

        isOuterInterfaceNode = falses(Nx, Ny)   # row just outside + row on fine grid boundary
        isInnerInterfaceNode = falses(Nx, Ny)   # second row inside fine grid
        isFineInterior       = falses(Nx, Ny)   # coarse nodes covered by fine grid → solid

        @inbounds for ix in 1:Nx, iy in 1:Ny
            !isFluid[ix, iy] && continue

            in_x   = ix_fine_left <= ix <= ix_fine_right
            in_y   = iy_fine_bottom <= iy <= iy_fine_top
            inside = in_x && in_y

            if inside
                minDist = min(ix - ix_fine_left, ix_fine_right - ix,
                              iy - iy_fine_bottom, iy_fine_top - iy)
                if minDist == 0
                    isOuterInterfaceNode[ix, iy] = true   # ON boundary = first embedded row
                elseif minDist == 1
                    isInnerInterfaceNode[ix, iy] = true   # second embedded row
                else
                    isFineInterior[ix, iy] = true          # deep inside → replaced by fine grid
                end
            else
                # Just outside: within Chebyshev distance 1 of the fine grid rectangle
                near_x = (ix_fine_left - 1) <= ix <= (ix_fine_right + 1)
                near_y = (iy_fine_bottom - 1) <= iy <= (iy_fine_top + 1)
                if near_x && near_y
                    isOuterInterfaceNode[ix, iy] = true   # just-outside row
                end
            end
        end

        # Remove fine-interior coarse nodes from fluid (replaced by fine grid)
        isFluid .&= .!isFineInterior

        # Precomputed node index lists — rebuild after any mask change (e.g. adding a solid object)
        fluidNodes = findall(isFluid)
        solidNodes = findall(.!isFluid .& .!isInlet .& .!isOutlet .& .!isWall)
        objectNodes = findall(isObject)
        
        boundaryNodesAndDistances = find_object_boundary_nodes(isObject, fluidNodes,
                                                        deltaX, positionX, positionY,
                                                        1.5*d, 0.5*d, cosd(angleDeg), sind(angleDeg))

        ##-------- MRT Setup --------##
        omegaBGK      = 1.0 / (3.0 * latticeViscosity + 0.5)
        omegaAcoustic = 1.0
        ##-------- Fine Grid MRT Setup --------##
        omegaBGKFine = 1.0 / (2.0 / omegaBGK - 0.5)   # acoustic scaling: (τF-0.5)=2(τC-0.5)

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
        inlet_add_fp0 = (2.0 / (9.0  * latticeSpeedOfSound^2)) * latticeInflowVelocity
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

            # Object bounce-back (Bouzidi) + momentum exchange
            forceX = 0.0
            forceY = 0.0
            @inbounds for (ix, iy, cx, cy, q) in boundaryNodesAndDistances
                q2   = 2.0 * q
                f_in = 0.0
                f_out = 0.0
                if q < 0.5
                    # f_ᾱ[ix,iy] = 2q*f_α[ix+cx,iy+cy] + (1-2q)*f_α[ix,iy]
                    if     cx ==  1 && cy ==  0;  f_in = fp0[ix+1,iy  ]; f_out = q2*f_in + (1.0-q2)*fp0[ix,  iy  ]; fm0[ix,iy] = f_out
                    elseif cx == -1 && cy ==  0;  f_in = fm0[ix-1,iy  ]; f_out = q2*f_in + (1.0-q2)*fm0[ix,  iy  ]; fp0[ix,iy] = f_out
                    elseif cx ==  0 && cy ==  1;  f_in = f0p[ix,  iy+1]; f_out = q2*f_in + (1.0-q2)*f0p[ix,  iy  ]; f0m[ix,iy] = f_out
                    elseif cx ==  0 && cy == -1;  f_in = f0m[ix,  iy-1]; f_out = q2*f_in + (1.0-q2)*f0m[ix,  iy  ]; f0p[ix,iy] = f_out
                    elseif cx ==  1 && cy ==  1;  f_in = fpp[ix+1,iy+1]; f_out = q2*f_in + (1.0-q2)*fpp[ix,  iy  ]; fmm[ix,iy] = f_out
                    elseif cx ==  1 && cy == -1;  f_in = fpm[ix+1,iy-1]; f_out = q2*f_in + (1.0-q2)*fpm[ix,  iy  ]; fmp[ix,iy] = f_out
                    elseif cx == -1 && cy ==  1;  f_in = fmp[ix-1,iy+1]; f_out = q2*f_in + (1.0-q2)*fmp[ix,  iy  ]; fpm[ix,iy] = f_out
                    elseif cx == -1 && cy == -1;  f_in = fmm[ix-1,iy-1]; f_out = q2*f_in + (1.0-q2)*fmm[ix,  iy  ]; fpp[ix,iy] = f_out
                    end
                else
                    iq2 = 1.0 / q2
                    r   = (q2 - 1.0) * iq2   # (2q-1)/(2q)
                    # f_ᾱ[ix,iy] = (1/2q)*f_α[ix+cx,iy+cy] + ((2q-1)/2q)*f_ᾱ[ix-cx,iy-cy]
                    if     cx ==  1 && cy ==  0;  f_in = fp0[ix+1,iy  ]; f_out = iq2*f_in + r*fm0[ix-1,iy  ]; fm0[ix,iy] = f_out
                    elseif cx == -1 && cy ==  0;  f_in = fm0[ix-1,iy  ]; f_out = iq2*f_in + r*fp0[ix+1,iy  ]; fp0[ix,iy] = f_out
                    elseif cx ==  0 && cy ==  1;  f_in = f0p[ix,  iy+1]; f_out = iq2*f_in + r*f0m[ix,  iy-1]; f0m[ix,iy] = f_out
                    elseif cx ==  0 && cy == -1;  f_in = f0m[ix,  iy-1]; f_out = iq2*f_in + r*f0p[ix,  iy+1]; f0p[ix,iy] = f_out
                    elseif cx ==  1 && cy ==  1;  f_in = fpp[ix+1,iy+1]; f_out = iq2*f_in + r*fmm[ix-1,iy-1]; fmm[ix,iy] = f_out
                    elseif cx ==  1 && cy == -1;  f_in = fpm[ix+1,iy-1]; f_out = iq2*f_in + r*fmp[ix-1,iy+1]; fmp[ix,iy] = f_out
                    elseif cx == -1 && cy ==  1;  f_in = fmp[ix-1,iy+1]; f_out = iq2*f_in + r*fpm[ix+1,iy-1]; fpm[ix,iy] = f_out
                    elseif cx == -1 && cy == -1;  f_in = fmm[ix-1,iy-1]; f_out = iq2*f_in + r*fpp[ix+1,iy+1]; fpp[ix,iy] = f_out
                    end
                end
                dfx, dfy = _momentum_exchange(cx, cy, f_in, f_out)
                forceX += dfx
                forceY += dfy
            end

            ##-- 4. Logging & Plotting --##
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

        Log_Simulation_Tail()
    end#run

end#JuLattice