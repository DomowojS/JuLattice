module GridRefinement
using ..Kernel: getNonEquilibrium, _node_macros_and_stress
export InterpolationCoef, _synchronize!

struct InterpolationCoef
    # u(x,y) = a0 + ax*x + ay*y + axy*x*y + axx*x^2 + ayy*y^2
    a0::Float64;  ax::Float64;  ay::Float64;  axy::Float64;  axx::Float64;  ayy::Float64
    # v(x,y) = b0 + bx*x + by*y + bxy*x*y + bxx*x^2 + byy*y^2
    b0::Float64;  bx::Float64;  by::Float64;  bxy::Float64;  bxx::Float64;  byy::Float64
    # rho(x,y) = c0 + cx*x + cy*y + cxy*x*y
    c0::Float64;  cx::Float64;  cy::Float64;  cxy::Float64
end

function _get_interpolation_coef(ix, iy, omegaBGK,
                                 f00, fp0, fm0, f0p, f0m,
                                 fpp, fpm, fmp, fmm)
    rho00, u00, v00, c20_00, c02_00, c11_00 = _node_macros_and_stress(ix,   iy,   f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm)
    rho10, u10, v10, c20_10, c02_10, c11_10 = _node_macros_and_stress(ix+1, iy,   f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm)
    rho11, u11, v11, c20_11, c02_11, c11_11 = _node_macros_and_stress(ix+1, iy+1, f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm)
    rho01, u01, v01, c20_01, c02_01, c11_01 = _node_macros_and_stress(ix,   iy+1, f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm)

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

@inline function _eval_and_write_fine!(ixF, iyF, coef::InterpolationCoef,
                                        xx::Float64, yy::Float64, omegaBGK::Float64,
                                        f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine,
                                        fppFine, fpmFine, fmpFine, fmmFine)
    ux   = coef.a0 + coef.ax*xx + coef.ay*yy + coef.axx*xx*xx + coef.ayy*yy*yy + coef.axy*xx*yy
    uy   = coef.b0 + coef.bx*xx + coef.by*yy + coef.bxx*xx*xx + coef.byy*yy*yy + coef.bxy*xx*yy
    rho0 = coef.c0 + coef.cx*xx + coef.cy*yy + coef.cxy*xx*yy
    cxy  = -1.0/(3.0*omegaBGK) * ((coef.ay + 2.0*coef.ayy*yy + coef.axy*xx) +
                                    (coef.bx + 2.0*coef.bxx*xx + coef.bxy*yy)) * 0.5
    cxx  = -2.0/(3.0*omegaBGK) *   (coef.ax + 2.0*coef.axx*xx + coef.axy*yy) * 0.5
    cyy  = -2.0/(3.0*omegaBGK) *   (coef.by + 2.0*coef.byy*yy + coef.bxy*xx) * 0.5
    neq  = getNonEquilibrium(rho0, ux, uy, cxx, cyy, cxy)
    f00Fine[ixF,iyF] = neq.f00; fp0Fine[ixF,iyF] = neq.fp0; fm0Fine[ixF,iyF] = neq.fm0
    f0pFine[ixF,iyF] = neq.f0p; f0mFine[ixF,iyF] = neq.f0m
    fppFine[ixF,iyF] = neq.fpp; fpmFine[ixF,iyF] = neq.fpm
    fmpFine[ixF,iyF] = neq.fmp; fmmFine[ixF,iyF] = neq.fmm
end

function _synchronize!(
        innerInterfaceNodes, outerInterfaceNodes,
        innerInterfaceNodesFine, outerInterfaceNodesFine,
        ix_left, iy_bottom,
        omegaBGK, omegaBGKFine,
        f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm,
        f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)

    # ── F → C: fine inner interface → coarse inner interface ────────────────
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

    # ── C → F: coarse outer interface → fine outer interface ────────────────
    NxFine, NyFine = size(f00Fine)
    nCoarseX = (NxFine - 2) ÷ 2
    nCoarseY = (NyFine - 2) ÷ 2
    ix_right = ix_left + nCoarseX
    iy_top   = iy_bottom + nCoarseY

    # Bottom edge: coarse iy=iy_bottom → fine iyF=2,3
    for ix in ix_left:ix_right-1
        coef  = _get_interpolation_coef(ix, iy_bottom, omegaBGK,
                                        f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm)
        ixF_L = 2*(ix - ix_left) + 2
        _eval_and_write_fine!(ixF_L,   2, coef, -0.25, -0.25, omegaBGK, f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)
        _eval_and_write_fine!(ixF_L+1, 2, coef,  0.25, -0.25, omegaBGK, f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)
        _eval_and_write_fine!(ixF_L+1, 3, coef,  0.25,  0.25, omegaBGK, f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)
        _eval_and_write_fine!(ixF_L,   3, coef, -0.25,  0.25, omegaBGK, f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)
    end

    # Top edge: coarse iy=iy_top-1 → fine iyF=NyFine-2,NyFine-1
    for ix in ix_left:ix_right-1
        coef  = _get_interpolation_coef(ix, iy_top-1, omegaBGK,
                                        f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm)
        ixF_L = 2*(ix - ix_left) + 2
        _eval_and_write_fine!(ixF_L,   NyFine-2, coef, -0.25, -0.25, omegaBGK, f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)
        _eval_and_write_fine!(ixF_L+1, NyFine-2, coef,  0.25, -0.25, omegaBGK, f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)
        _eval_and_write_fine!(ixF_L+1, NyFine-1, coef,  0.25,  0.25, omegaBGK, f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)
        _eval_and_write_fine!(ixF_L,   NyFine-1, coef, -0.25,  0.25, omegaBGK, f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)
    end

    # Left edge: coarse ix=ix_left → fine ixF=2,3
    for iy in iy_bottom:iy_top-1
        coef  = _get_interpolation_coef(ix_left, iy, omegaBGK,
                                        f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm)
        iyF_B = 2*(iy - iy_bottom) + 2
        _eval_and_write_fine!(2, iyF_B,   coef, -0.25, -0.25, omegaBGK, f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)
        _eval_and_write_fine!(3, iyF_B,   coef,  0.25, -0.25, omegaBGK, f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)
        _eval_and_write_fine!(3, iyF_B+1, coef,  0.25,  0.25, omegaBGK, f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)
        _eval_and_write_fine!(2, iyF_B+1, coef, -0.25,  0.25, omegaBGK, f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)
    end

    # Right edge: coarse ix=ix_right-1 → fine ixF=NxFine-2,NxFine-1
    for iy in iy_bottom:iy_top-1
        coef  = _get_interpolation_coef(ix_right-1, iy, omegaBGK,
                                        f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm)
        iyF_B = 2*(iy - iy_bottom) + 2
        _eval_and_write_fine!(NxFine-2, iyF_B,   coef, -0.25, -0.25, omegaBGK, f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)
        _eval_and_write_fine!(NxFine-1, iyF_B,   coef,  0.25, -0.25, omegaBGK, f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)
        _eval_and_write_fine!(NxFine-1, iyF_B+1, coef,  0.25,  0.25, omegaBGK, f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)
        _eval_and_write_fine!(NxFine-2, iyF_B+1, coef, -0.25,  0.25, omegaBGK, f00Fine, fp0Fine, fm0Fine, f0pFine, f0mFine, fppFine, fpmFine, fmpFine, fmmFine)
    end
end

end # module GridRefinement
