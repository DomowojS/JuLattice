module Kernel
export getEquilibrium, _collide_and_stream!, _momentum_exchange, _apply_bouzidi_bc!

function getEquilibrium(rho::Float64, u::Float64, v::Float64, cx::Int, cy::Int)
    ux2 = u * u
    uy2 = v * v
    if cx == -1 && cy == -1
        return rho * (1 - 3*u + 3*ux2) * (1 - 3*v + 3*uy2) / 36
    elseif cx ==  0 && cy == -1
        return -rho * (-2 + 3*ux2) * (1 + 3*uy2 - 3*v) / 18
    elseif cx ==  1 && cy == -1
        return rho * (1 + 3*ux2 + 3*u) * (1 + 3*uy2 - 3*v) / 36
    elseif cx == -1 && cy ==  0
        return -rho * (-2 + 3*uy2) * (1 + 3*ux2 - 3*u) / 18
    elseif cx ==  0 && cy ==  0
        return rho * (-2 + 3*ux2) * (-2 + 3*uy2) / 9
    elseif cx ==  1 && cy ==  0
        return -rho * (-2 + 3*uy2) * (1 + 3*ux2 + 3*u) / 18
    elseif cx == -1 && cy ==  1
        return rho * (1 + 3*ux2 - 3*u) * (1 + 3*uy2 + 3*v) / 36
    elseif cx ==  0 && cy ==  1
        return -rho * (-2 + 3*ux2) * (1 + 3*uy2 + 3*v) / 18
    elseif cx ==  1 && cy ==  1
        return rho * (1 + 3*ux2 + 3*u) * (1 + 3*uy2 + 3*v) / 36
    else
        error("Invalid lattice direction (cx,cy)=($cx,$cy)")
    end
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

@inline function _momentum_exchange(cx, cy, f_in, f_out)
    s = f_in + f_out
    return cx * s, cy * s
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

end # module Kernel
