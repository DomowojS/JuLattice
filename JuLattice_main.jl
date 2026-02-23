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

    @inline function getNonEquilibrium(rho0, ux, uy, cxx, cyy, cxy)
        m11 = (cxy + ux*uy) * rho0
        m20 = (cxx + 1.0/3.0 + ux*ux) * rho0
        m02 = (cyy + 1.0/3.0 + uy*uy) * rho0
        m12 = ux * m02
        m21 = uy * m20
        m22 = m20 * m02 / rho0
        f00 =        rho0 - m02 - m20 + m22
        fp0 = 0.5  * ( ux*rho0 + m20 - m12 - m22)
        fm0 = 0.5  * (-ux*rho0 + m20 + m12 - m22)
        f0p = 0.5  * ( uy*rho0 + m02 - m21 - m22)
        f0m = 0.5  * (-uy*rho0 + m02 + m21 - m22)
        fpp = 0.25 * ( m11 + m12 + m21 + m22)
        fpm = 0.25 * (-m11 + m12 - m21 + m22)
        fmp = 0.25 * (-m11 - m12 + m21 + m22)
        fmm = 0.25 * ( m11 - m12 - m21 + m22)
        return (; f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm)
    end

    function add_rectangle!(isObject, Nx, Ny, deltaX;
                            centerX, centerY, d, angleDeg,
                            originX=0.0, originY=0.0)
        half_len = 1.5 * d   # half of 3d  (long axis)
        half_hgt = 0.5 * d   # half of d   (short axis)

        alpha = deg2rad(angleDeg)
        c = cos(alpha)
        s = sin(alpha)

        for ix in 1:Nx, iy in 1:Ny
            # Physical coordinates: ix=2,iy=2 → x=originX, y=originY
            x = originX + (ix - 2) * deltaX
            y = originY + (iy - 2) * deltaX

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
                                        half_len, half_hgt, cos_a, sin_a;
                                        originX=0.0, originY=0.0)
        dirs = ((1,0), (-1,0), (0,1), (0,-1), (1,1), (1,-1), (-1,1), (-1,-1))

        boundaryNodesAndDistances = Tuple{Int,Int,Int,Int,Float64}[]

        @inbounds for idx in fluidNodes
            ix, iy = Tuple(idx)
            # Fluid node in local frame
            dx0 = originX + (ix - 2) * deltaX - center_x
            dy0 = originY + (iy - 2) * deltaX - center_y
            lx0 =  cos_a * dx0 - sin_a * dy0
            ly0 =  sin_a * dx0 + cos_a * dy0

            for (cx, cy) in dirs
                if isObject[ix + cx, iy + cy]
                    # Object node in local frame
                    dx1 = originX + (ix + cx - 2) * deltaX - center_x
                    dy1 = originY + (iy + cy - 2) * deltaX - center_y
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

    function _apply_bouzidi_bc!(boundaryNodesAndDistances,
                                f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm)
        forceX = 0.0
        forceY = 0.0
        @inbounds for (ix, iy, cx, cy, q) in boundaryNodesAndDistances
            q2   = 2.0 * q
            f_in = 0.0
            f_out = 0.0
            if q < 0.5
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
                r   = (q2 - 1.0) * iq2
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
        return forceX, forceY
    end

    function _classify_coarse_nodes(Nx, Ny, deltaX,
                                    positionFineGridX, positionFineGridY,
                                    lengthXFine, lengthYFine)
        isInlet  = falses(Nx, Ny);  isInlet[1,      :]         .= true
        isOutlet = falses(Nx, Ny);  isOutlet[Nx,    :]         .= true
        isWall   = falses(Nx, Ny);  isWall[2:Nx-1, [1, Ny]]    .= true
        isFluid  = falses(Nx, Ny);  isFluid[2:Nx-1,  2:Ny-1]   .= true
        isObject = falses(Nx, Ny)   # object lives on fine grid; coarse has none
        isFluid .&= .!isObject
        isSolid  = isInlet .| isOutlet .| isWall .| isObject

        # Map fine grid physical bounds onto coarse indices
        # (requires positionFineGridX/deltaX to be an integer)
        ix_left   = 2 + round(Int, positionFineGridX / deltaX)
        ix_right  = 2 + round(Int, (positionFineGridX + lengthXFine) / deltaX)
        iy_bottom = 2 + round(Int, positionFineGridY / deltaX)
        iy_top    = 2 + round(Int, (positionFineGridY + lengthYFine) / deltaX)

        isOuterInterfaceNode = falses(Nx, Ny)
        isInnerInterfaceNode = falses(Nx, Ny)
        isFineInterior       = falses(Nx, Ny)

        @inbounds for ix in 1:Nx, iy in 1:Ny
            !isFluid[ix, iy] && continue

            in_x   = ix_left <= ix <= ix_right
            in_y   = iy_bottom <= iy <= iy_top
            inside = in_x && in_y

            if inside
                minDist = min(ix - ix_left, ix_right - ix, iy - iy_bottom, iy_top - iy)
                if minDist <= 1
                    isOuterInterfaceNode[ix, iy] = true
                elseif minDist == 2
                    isInnerInterfaceNode[ix, iy] = true
                else
                    isFineInterior[ix, iy] = true
                end
            end
        end

        isFluid .&= .!isFineInterior   # coarse nodes inside fine region are no longer fluid

        return (; isInlet, isOutlet, isWall, isFluid, isObject, isSolid,
                  isOuterInterfaceNode, isInnerInterfaceNode, isFineInterior,
                  ix_left, ix_right, iy_bottom, iy_top)
    end

    function _classify_fine_nodes(NxFine, NyFine, deltaXFine,
                                   positionFineGridX, positionFineGridY,
                                   positionX, positionY, d, angleDeg)
        # Fine node (ixF,iyF) physical position: x = originX + (ixF-2)*deltaXFine
        # Offset by deltaXFine/2 from anchor → staggered w.r.t. coarse nodes
        originXFine = positionFineGridX + 0.5 * deltaXFine
        originYFine = positionFineGridY + 0.5 * deltaXFine

        isFluidFine  = falses(NxFine, NyFine);  isFluidFine[2:NxFine-1, 2:NyFine-1] .= true
        isObjectFine = falses(NxFine, NyFine)
        add_rectangle!(isObjectFine, NxFine, NyFine, deltaXFine;
                       centerX=positionX, centerY=positionY, d=d, angleDeg=angleDeg,
                       originX=originXFine, originY=originYFine)
        isFluidFine .&= .!isObjectFine

        # Outer 2 rows/cols of fine interior → C→F target (isOuterInterfaceNodeFine)
        isOuterInterfaceNodeFine = falses(NxFine, NyFine)
        isOuterInterfaceNodeFine[2:NxFine-1, 2:NyFine-1] .= true   # all interior
        isOuterInterfaceNodeFine[4:NxFine-3, 4:NyFine-3] .= false  # clear from 3rd row inward
        isOuterInterfaceNodeFine .&= isFluidFine

        # 4th and 5th rows/cols of fine interior → F→C source (isInnerInterfaceNodeFine)
        # Fill from 4th row/col inward, then clear from 6th row/col inward → ring at rows/cols 4–5
        isInnerInterfaceNodeFine = falses(NxFine, NyFine)
        isInnerInterfaceNodeFine[5:NxFine-4, 5:NyFine-4] .= true
        isInnerInterfaceNodeFine[7:NxFine-6, 7:NyFine-6] .= false
        isInnerInterfaceNodeFine .&= isFluidFine

        return (; isFluidFine, isObjectFine,
                  isOuterInterfaceNodeFine, isInnerInterfaceNodeFine,
                  originXFine, originYFine)
    end

    function _compute_node_lists(
            isFluid, isObject, isInlet, isOutlet, isWall,
            isOuterInterfaceNode, isInnerInterfaceNode,
            isFluidFine, isObjectFine, isOuterInterfaceNodeFine, isInnerInterfaceNodeFine,
            deltaX, deltaXFine, positionX, positionY, d, angleDeg,
            originXFine, originYFine)

        # Coarse grid
        fluidNodes          = findall(isFluid)
        solidNodes          = findall(.!isFluid .& .!isInlet .& .!isOutlet .& .!isWall)
        objectNodes         = findall(isObject)
        outerInterfaceNodes = findall(isOuterInterfaceNode)
        innerInterfaceNodes = findall(isInnerInterfaceNode)

        # Fine grid
        fluidNodesFine          = findall(isFluidFine)
        objectNodesFine         = findall(isObjectFine)
        outerInterfaceNodesFine = findall(isOuterInterfaceNodeFine)
        innerInterfaceNodesFine = findall(isInnerInterfaceNodeFine)

        # Bouzidi boundary — coarse (empty; object on fine grid)
        boundaryNodesAndDistances = find_object_boundary_nodes(
            isObject, fluidNodes,
            deltaX, positionX, positionY,
            1.5*d, 0.5*d, cosd(angleDeg), sind(angleDeg))

        # Bouzidi boundary — fine grid
        boundaryNodesAndDistancesFine = find_object_boundary_nodes(
            isObjectFine, fluidNodesFine,
            deltaXFine, positionX, positionY,
            1.5*d, 0.5*d, cosd(angleDeg), sind(angleDeg),
            originX=originXFine, originY=originYFine)

        return (; fluidNodes, solidNodes, objectNodes,
                  outerInterfaceNodes, innerInterfaceNodes,
                  fluidNodesFine, objectNodesFine,
                  outerInterfaceNodesFine, innerInterfaceNodesFine,
                  boundaryNodesAndDistances, boundaryNodesAndDistancesFine)
    end

    struct InterpolationCoef
        # u(x,y) = a0 + ax*x + ay*y + axy*x*y + axx*x^2 + ayy*y^2
        a0::Float64;  ax::Float64;  ay::Float64;  axy::Float64;  axx::Float64;  ayy::Float64
        # v(x,y) = b0 + bx*x + by*y + bxy*x*y + bxx*x^2 + byy*y^2
        b0::Float64;  bx::Float64;  by::Float64;  bxy::Float64;  bxx::Float64;  byy::Float64
        # rho(x,y) = c0 + cx*x + cy*y + cxy*x*y
        c0::Float64;  cx::Float64;  cy::Float64;  cxy::Float64
    end

    @inline function _node_macros_and_stress(ix, iy,
                                              f00, fp0, fm0, f0p, f0m,
                                              fpp, fpm, fmp, fmm)
        rho  = f00[ix,iy] + fp0[ix,iy] + fm0[ix,iy] + f0p[ix,iy] + f0m[ix,iy] +
               fpp[ix,iy] + fpm[ix,iy] + fmp[ix,iy] + fmm[ix,iy]
        irho = 1.0 / rho
        u    = (fp0[ix,iy] - fm0[ix,iy] + fpp[ix,iy] + fpm[ix,iy] - fmp[ix,iy] - fmm[ix,iy]) * irho
        v    = (f0p[ix,iy] - f0m[ix,iy] + fpp[ix,iy] - fpm[ix,iy] + fmp[ix,iy] - fmm[ix,iy]) * irho
        m20  = fm0[ix,iy] + fp0[ix,iy] + fmm[ix,iy] + fpp[ix,iy] + fmp[ix,iy] + fpm[ix,iy]
        m02  = f0m[ix,iy] + f0p[ix,iy] + fmm[ix,iy] + fpp[ix,iy] + fmp[ix,iy] + fpm[ix,iy]
        m11  = fmm[ix,iy] + fpp[ix,iy] - fmp[ix,iy] - fpm[ix,iy]
        c20  = m20 * irho - u*u   # cxx + 1/3
        c02  = m02 * irho - v*v   # cyy + 1/3
        c11  = m11 * irho - u*v   # cxy
        return rho, u, v, c20, c02, c11
    end

    function _get_interpolation_coef(ix, iy, omegaBGK,
                                     f00, fp0, fm0, f0p, f0m,
                                     fpp, fpm, fmp, fmm)
        rho00, u00, v00, c20_00, c02_00, c11_00 = _node_macros_and_stress(ix,   iy,   f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm)
        rho10, u10, v10, c20_10, c02_10, c11_10 = _node_macros_and_stress(ix+1, iy,   f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm)
        rho11, u11, v11, c20_11, c02_11, c11_11 = _node_macros_and_stress(ix+1, iy+1, f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm)
        rho01, u01, v01, c20_01, c02_01, c11_01 = _node_macros_and_stress(ix,   iy+1, f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm)

        # Row-averaged finite differences of stress (shear modes only; 1/3 terms cancel in c20-c02)
        DxC11 = -3.0*omegaBGK * (c11_10 + c11_11 - c11_00 - c11_01) * 0.5
        DyC11 = -3.0*omegaBGK * (c11_01 + c11_11 - c11_00 - c11_10) * 0.5
        DxC2  = -3.0*omegaBGK * 0.5 * ((c20_10 - c02_10 + c20_11 - c02_11) -
                                         (c20_00 - c02_00 + c20_01 - c02_01)) * 0.5
        DyC2  = -3.0*omegaBGK * 0.5 * ((c20_01 - c02_01 + c20_11 - c02_11) -
                                         (c20_00 - c02_00 + c20_10 - c02_10)) * 0.5

        a0  = (-DxC2  - DyC11 + 2.0*(u00 + u01 + u10 + u11)) / 8.0
        ax  = (-u00 - u01 + u10 + u11) / 2.0
        ay  = (-u00 + u01 - u10 + u11) / 2.0
        axy = u00 - u01 - u10 + u11
        axx = ( DxC2  + v00 - v01 - v10 + v11) / 2.0
        ayy = ( DyC11 - v00 + v01 + v10 - v11) / 2.0

        b0  = (-DxC11 + DyC2  + 2.0*(v00 + v01 + v10 + v11)) / 8.0
        bx  = (-v00 - v01 + v10 + v11) / 2.0
        by  = (-v00 + v01 - v10 + v11) / 2.0
        bxy = v00 - v01 - v10 + v11
        bxx = ( DxC11 - u00 + u01 + u10 - u11) / 2.0
        byy = (-DyC2  + u00 - u01 - u10 + u11) / 2.0

        c0  = (rho00 + rho01 + rho10 + rho11) / 4.0
        c_x = (-rho00 - rho01 + rho10 + rho11) / 2.0
        c_y = (-rho00 + rho01 - rho10 + rho11) / 2.0
        c_xy = rho00 - rho01 - rho10 + rho11

        return InterpolationCoef(a0, ax, ay, axy, axx, ayy,
                                 b0, bx, by, bxy, bxx, byy,
                                 c0, c_x, c_y, c_xy)
    end

    function _synchronize!(
            innerInterfaceNodes, outerInterfaceNodes,
            innerInterfaceNodesFine, outerInterfaceNodesFine,
            ix_left, iy_bottom,
            omegaBGK, omegaBGKFine,
            f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm,
            f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)

        # ── F → C: fine inner interface → coarse inner interface ────────────────
        # For each coarse inner interface node (ix,iy), the 4 surrounding fine nodes
        # form a 2×2 cell with bottom-left at ixF = 2*(ix-ix_left)+1, iyF = 2*(iy-iy_bottom)+1.
        # We call _get_interpolation_coef at that fine cell, evaluate at center (a0,b0,c0),
        # apply CE stress inversion with scaleUP=2, and reconstruct coarse distributions.
        for idx in innerInterfaceNodes
            ix, iy = Tuple(idx)
            ixF = 2*(ix - ix_left) + 1
            iyF = 2*(iy - iy_bottom) + 1

            coef = _get_interpolation_coef(ixF, iyF, omegaBGKFine,
                                           f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine,
                                           fppFine, fpmFine, fmpFine, fmmFine)

            ux   = coef.a0
            uy   = coef.b0
            rho0 = coef.c0
            cxy  = -1.0/(3.0*omegaBGKFine) * (coef.ay + coef.bx) * 2.0
            cxx  = -2.0/(3.0*omegaBGKFine) * coef.ax * 2.0
            cyy  = -2.0/(3.0*omegaBGKFine) * coef.by * 2.0

            neq = getNonEquilibrium(rho0, ux, uy, cxx, cyy, cxy)
            f00[ix,iy] = neq.f00; fp0[ix,iy] = neq.fp0; fm0[ix,iy] = neq.fm0
            f0p[ix,iy] = neq.f0p; f0m[ix,iy] = neq.f0m
            fpp[ix,iy] = neq.fpp; fpm[ix,iy] = neq.fpm; fmp[ix,iy] = neq.fmp; fmm[ix,iy] = neq.fmm
        end

        # ── C → F: coarse outer interface → fine outer interface (stub) ─────────
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