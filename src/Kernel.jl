module Kernel

using ..TurbulenceModel: compute_pi_norm, smagorinsky_omega

export collision_stream!
export Q000, QM00, QP00, Q0M0, Q0P0, Q00M, Q00P
export QMM0, QMP0, QPM0, QPP0
export QM0M, QM0P, QP0M, QP0P
export Q0MM, Q0MP, Q0PM, Q0PP
export NQ

# D3Q19 direction indices — first index of f[q, x, y, z]
const Q000 = 1   # rest         ( 0, 0, 0)
const QM00 = 2   # (-1, 0, 0)
const QP00 = 3   # (+1, 0, 0)
const Q0M0 = 4   # ( 0,-1, 0)
const Q0P0 = 5   # ( 0,+1, 0)
const Q00M = 6   # ( 0, 0,-1)
const Q00P = 7   # ( 0, 0,+1)
const QMM0 = 8   # (-1,-1, 0)
const QMP0 = 9   # (-1,+1, 0)
const QPM0 = 10  # (+1,-1, 0)
const QPP0 = 11  # (+1,+1, 0)
const QM0M = 12  # (-1, 0,-1)
const QM0P = 13  # (-1, 0,+1)
const QP0M = 14  # (+1, 0,-1)
const QP0P = 15  # (+1, 0,+1)
const Q0MM = 16  # ( 0,-1,-1)
const Q0MP = 17  # ( 0,-1,+1)
const Q0PM = 18  # ( 0,+1,-1)
const Q0PP = 19  # ( 0,+1,+1)
const NQ   = 19

#with rho
# function collision_stream!(
#     gridlengthX::Int, gridlengthY::Int, gridlengthZ::Int,
#     τ::Float64, CS::Float64, is_fluid,
#     rho, u, v, w,
#     f, fS, disc_nodes, F_x_lat
# )

# without rho
function collision_stream!(
    gridlengthX::Int, gridlengthY::Int, gridlengthZ::Int,
    τ::Float64, CS::Float64, is_fluid,
    u, v, w,
    f, fS, disc_nodes, F_x_lat
)

    # Iterate over all cells except boundary cells
    @inbounds Threads.@threads :static for z in 2:gridlengthZ-1
        for y in 2:gridlengthY-1
            for x in 2:gridlengthX-1

                if !is_fluid[x,y,z]
                    continue
                end

                # Load all 19 populations from local node
                f000 = f[Q000, x, y, z]
                fm00 = f[QM00, x, y, z]
                fp00 = f[QP00, x, y, z]
                f0m0 = f[Q0M0, x, y, z]
                f0p0 = f[Q0P0, x, y, z]
                f00m = f[Q00M, x, y, z]
                f00p = f[Q00P, x, y, z]
                fmm0 = f[QMM0, x, y, z]
                fmp0 = f[QMP0, x, y, z]
                fpm0 = f[QPM0, x, y, z]
                fpp0 = f[QPP0, x, y, z]
                fm0m = f[QM0M, x, y, z]
                fm0p = f[QM0P, x, y, z]
                fp0m = f[QP0M, x, y, z]
                fp0p = f[QP0P, x, y, z]
                f0mm = f[Q0MM, x, y, z]
                f0mp = f[Q0MP, x, y, z]
                f0pm = f[Q0PM, x, y, z]
                f0pp = f[Q0PP, x, y, z]

                # Compute macroscopic quantities
                rho_loc = f000 +
                            (fm00 + fp00 + f0m0 + f0p0 + f00m + f00p) +
                            (fmm0 + fmp0 + fpm0 + fpp0 +
                            fm0m + fm0p + fp0m + fp0p +
                            f0mm + f0mp + f0pm + f0pp)

                #rho[x,y,z] = rho_loc
                inv_rho = 1.0 / rho_loc

                u_loc = ((-fm00 + fp00) +
                            (-fmm0 - fmp0 + fpm0 + fpp0) +
                            (-fm0m - fm0p + fp0m + fp0p)) * inv_rho

                v_loc = ((-f0m0 + f0p0) +
                            (-fmm0 + fmp0 - fpm0 + fpp0) +
                            (-f0mm - f0mp + f0pm + f0pp)) * inv_rho

                w_loc = ((-f00m + f00p) +
                            (-fm0m + fm0p - fp0m + fp0p) +
                            (-f0mm + f0mp - f0pm + f0pp)) * inv_rho

                u[x,y,z] = u_loc
                v[x,y,z] = v_loc
                w[x,y,z] = w_loc

                # Compute equilibrium
                # Pre-compute polynomial factors
                u2 = u_loc * u_loc
                v2 = v_loc * v_loc
                w2 = w_loc * w_loc

                Pm_u = 1 - 3 * u_loc + 3*u2
                P0_u = 1 - 1.5*u2
                Pp_u = 1 + 3 * u_loc + 3*u2

                Pm_v = 1 - 3 * v_loc + 3*v2
                P0_v = 1 - 1.5*v2
                Pp_v = 1 + 3*v_loc + 3*v2

                Pm_w = 1 - 3*w_loc + 3*w2
                P0_w = 1 - 1.5*w2
                Pp_w = 1 + 3*w_loc + 3*w2

                rho_inv3 = rho_loc * (1.0/3.0)
                rho_inv18 = rho_loc * (1.0/18.0)
                rho_inv36 = rho_loc * (1.0/36.0)

                # Compute equilibrium distribution feq
                feq000 = P0_u * P0_v * P0_w * rho_inv3

                # Face directions
                feqm00 = Pm_u * P0_v * P0_w * rho_inv18
                feqp00 = Pp_u * P0_v * P0_w * rho_inv18
                feq0m0 = P0_u * Pm_v * P0_w * rho_inv18
                feq0p0 = P0_u * Pp_v * P0_w * rho_inv18
                feq00m = P0_u * P0_v * Pm_w * rho_inv18
                feq00p = P0_u * P0_v * Pp_w * rho_inv18

                # XY edges
                feqmm0 = Pm_u * Pm_v * P0_w * rho_inv36
                feqmp0 = Pm_u * Pp_v * P0_w * rho_inv36
                feqpm0 = Pp_u * Pm_v * P0_w * rho_inv36
                feqpp0 = Pp_u * Pp_v * P0_w * rho_inv36
                # XZ edges
                feqm0m = Pm_u * P0_v * Pm_w * rho_inv36
                feqm0p = Pm_u * P0_v * Pp_w * rho_inv36
                feqp0m = Pp_u * P0_v * Pm_w * rho_inv36
                feqp0p = Pp_u * P0_v * Pp_w * rho_inv36
                # YZ edges
                feq0mm = P0_u * Pm_v * Pm_w * rho_inv36
                feq0mp = P0_u * Pm_v * Pp_w * rho_inv36
                feq0pm = P0_u * Pp_v * Pm_w * rho_inv36
                feq0pp = P0_u * Pp_v * Pp_w * rho_inv36

                # Compute Non-equilibrium distribution parts

                # not needed
                # fneq000 = f000 - feq000

                fneqm00 = fm00 - feqm00
                fneqp00 = fp00 - feqp00
                fneq0m0 = f0m0 - feq0m0
                fneq0p0 = f0p0 - feq0p0
                fneq00m = f00m - feq00m
                fneq00p = f00p - feq00p
                fneqmm0 = fmm0 - feqmm0
                fneqmp0 = fmp0 - feqmp0
                fneqpm0 = fpm0 - feqpm0
                fneqpp0 = fpp0 - feqpp0
                fneqm0m = fm0m - feqm0m
                fneqm0p = fm0p - feqm0p
                fneqp0m = fp0m - feqp0m
                fneqp0p = fp0p - feqp0p
                fneq0mm = f0mm - feq0mm
                fneq0mp = f0mp - feq0mp
                fneq0pm = f0pm - feq0pm
                fneq0pp = f0pp - feq0pp

                # compute Smagorinsky variables
                pi_neq_norm = compute_pi_norm(
                    fneqm00, fneqp00, fneq0m0, fneq0p0, fneq00m, fneq00p,
                    fneqmm0, fneqmp0, fneqpm0, fneqpp0,
                    fneqm0m, fneqm0p, fneqp0m, fneqp0p,
                    fneq0mm, fneq0mp, fneq0pm, fneq0pp
                )

                omega_local = smagorinsky_omega(τ, CS, pi_neq_norm, rho_loc)

                # Push scheme: Stream+Collision
                # Rest particle
                fS[Q000, x,   y,   z  ] = f000 + omega_local * (feq000 - f000)
                # Face neighbors
                fS[QM00, x-1, y,   z  ] = fm00 + omega_local * (feqm00 - fm00)
                fS[QP00, x+1, y,   z  ] = fp00 + omega_local * (feqp00 - fp00)
                fS[Q0M0, x,   y-1, z  ] = f0m0 + omega_local * (feq0m0 - f0m0)
                fS[Q0P0, x,   y+1, z  ] = f0p0 + omega_local * (feq0p0 - f0p0)
                fS[Q00M, x,   y,   z-1] = f00m + omega_local * (feq00m - f00m)
                fS[Q00P, x,   y,   z+1] = f00p + omega_local * (feq00p - f00p)
                # XY-plane edges
                fS[QMM0, x-1, y-1, z  ] = fmm0 + omega_local * (feqmm0 - fmm0)
                fS[QMP0, x-1, y+1, z  ] = fmp0 + omega_local * (feqmp0 - fmp0)
                fS[QPM0, x+1, y-1, z  ] = fpm0 + omega_local * (feqpm0 - fpm0)
                fS[QPP0, x+1, y+1, z  ] = fpp0 + omega_local * (feqpp0 - fpp0)
                # XZ-plane edges
                fS[QM0M, x-1, y,   z-1] = fm0m + omega_local * (feqm0m - fm0m)
                fS[QM0P, x-1, y,   z+1] = fm0p + omega_local * (feqm0p - fm0p)
                fS[QP0M, x+1, y,   z-1] = fp0m + omega_local * (feqp0m - fp0m)
                fS[QP0P, x+1, y,   z+1] = fp0p + omega_local * (feqp0p - fp0p)
                # YZ-plane edges
                fS[Q0MM, x,   y-1, z-1] = f0mm + omega_local * (feq0mm - f0mm)
                fS[Q0MP, x,   y-1, z+1] = f0mp + omega_local * (feq0mp - f0mp)
                fS[Q0PM, x,   y+1, z-1] = f0pm + omega_local * (feq0pm - f0pm)
                fS[Q0PP, x,   y+1, z+1] = f0pp + omega_local * (feq0pp - f0pp)





            end #end x
        end #end y
    end #end z

    @inbounds for idx in disc_nodes
        x, y, z = idx[1], idx[2], idx[3]

        # Reload f populations
        f000 = f[Q000, x, y, z]
        fm00 = f[QM00, x, y, z]
        fp00 = f[QP00, x, y, z]
        f0m0 = f[Q0M0, x, y, z]
        f0p0 = f[Q0P0, x, y, z]
        f00m = f[Q00M, x, y, z]
        f00p = f[Q00P, x, y, z]
        fmm0 = f[QMM0, x, y, z]
        fmp0 = f[QMP0, x, y, z]
        fpm0 = f[QPM0, x, y, z]
        fpp0 = f[QPP0, x, y, z]
        fm0m = f[QM0M, x, y, z]
        fm0p = f[QM0P, x, y, z]
        fp0m = f[QP0M, x, y, z]
        fp0p = f[QP0P, x, y, z]
        f0mm = f[Q0MM, x, y, z]
        f0mp = f[Q0MP, x, y, z]
        f0pm = f[Q0PM, x, y, z]
        f0pp = f[Q0PP, x, y, z]

        # Recompute macroscopic
        rho_loc = f000 +
            (fm00 + fp00 + f0m0 + f0p0 + f00m + f00p) +
            (fmm0 + fmp0 + fpm0 + fpp0 +
            fm0m + fm0p + fp0m + fp0p +
            f0mm + f0mp + f0pm + f0pp)

        inv_rho = 1.0 / rho_loc

        u_loc = ((-fm00 + fp00) +
                    (-fmm0 - fmp0 + fpm0 + fpp0) +
                    (-fm0m - fm0p + fp0m + fp0p)) * inv_rho

        u_loc = u_loc + 0.5 * F_x_lat * inv_rho

        v_loc = ((-f0m0 + f0p0) +
                    (-fmm0 + fmp0 - fpm0 + fpp0) +
                    (-f0mm - f0mp + f0pm + f0pp)) * inv_rho

        w_loc = ((-f00m + f00p) +
                    (-fm0m + fm0p - fp0m + fp0p) +
                    (-f0mm + f0mp - f0pm + f0pp)) * inv_rho

        # Recompute equilibrium
        # Pre-compute polynomial factors
        u2 = u_loc * u_loc
        v2 = v_loc * v_loc
        w2 = w_loc * w_loc

        Pm_u = 1 - 3 * u_loc + 3*u2
        P0_u = 1 - 1.5*u2
        Pp_u = 1 + 3 * u_loc + 3*u2

        Pm_v = 1 - 3 * v_loc + 3*v2
        P0_v = 1 - 1.5*v2
        Pp_v = 1 + 3*v_loc + 3*v2

        Pm_w = 1 - 3*w_loc + 3*w2
        P0_w = 1 - 1.5*w2
        Pp_w = 1 + 3*w_loc + 3*w2

        rho_inv3 = rho_loc * (1.0/3.0)
        rho_inv18 = rho_loc * (1.0/18.0)
        rho_inv36 = rho_loc * (1.0/36.0)

        # Recompute equilibrium distribution feq
        feq000 = P0_u * P0_v * P0_w * rho_inv3

        # Face directions
        feqm00 = Pm_u * P0_v * P0_w * rho_inv18
        feqp00 = Pp_u * P0_v * P0_w * rho_inv18
        feq0m0 = P0_u * Pm_v * P0_w * rho_inv18
        feq0p0 = P0_u * Pp_v * P0_w * rho_inv18
        feq00m = P0_u * P0_v * Pm_w * rho_inv18
        feq00p = P0_u * P0_v * Pp_w * rho_inv18

        # XY edges
        feqmm0 = Pm_u * Pm_v * P0_w * rho_inv36
        feqmp0 = Pm_u * Pp_v * P0_w * rho_inv36
        feqpm0 = Pp_u * Pm_v * P0_w * rho_inv36
        feqpp0 = Pp_u * Pp_v * P0_w * rho_inv36
        # XZ edges
        feqm0m = Pm_u * P0_v * Pm_w * rho_inv36
        feqm0p = Pm_u * P0_v * Pp_w * rho_inv36
        feqp0m = Pp_u * P0_v * Pm_w * rho_inv36
        feqp0p = Pp_u * P0_v * Pp_w * rho_inv36
        # YZ edges
        feq0mm = P0_u * Pm_v * Pm_w * rho_inv36
        feq0mp = P0_u * Pm_v * Pp_w * rho_inv36
        feq0pm = P0_u * Pp_v * Pm_w * rho_inv36
        feq0pp = P0_u * Pp_v * Pp_w * rho_inv36

        # Compute Non-equilibrium distribution parts

        # not needed
        # fneq000 = f000 - feq000

        fneqm00 = fm00 - feqm00
        fneqp00 = fp00 - feqp00
        fneq0m0 = f0m0 - feq0m0
        fneq0p0 = f0p0 - feq0p0
        fneq00m = f00m - feq00m
        fneq00p = f00p - feq00p
        fneqmm0 = fmm0 - feqmm0
        fneqmp0 = fmp0 - feqmp0
        fneqpm0 = fpm0 - feqpm0
        fneqpp0 = fpp0 - feqpp0
        fneqm0m = fm0m - feqm0m
        fneqm0p = fm0p - feqm0p
        fneqp0m = fp0m - feqp0m
        fneqp0p = fp0p - feqp0p
        fneq0mm = f0mm - feq0mm
        fneq0mp = f0mp - feq0mp
        fneq0pm = f0pm - feq0pm
        fneq0pp = f0pp - feq0pp

        # compute Smagorinsky variables
        pi_neq_norm = compute_pi_norm(
            fneqm00, fneqp00, fneq0m0, fneq0p0, fneq00m, fneq00p,
            fneqmm0, fneqmp0, fneqpm0, fneqpp0,
            fneqm0m, fneqm0p, fneqp0m, fneqp0p,
            fneq0mm, fneq0mp, fneq0pm, fneq0pp
        )

        omega_local = smagorinsky_omega(τ, CS, pi_neq_norm, rho_loc)
        # tau_local = 1.0 / omega_local          # τ_eff (already = τ̄ = τ_raw + 0.5)
        # tau_dash = tau_local + 0.5             # WRONG: adds 0.5 again — τ̄ already includes it
        # omega_dash = 1.0 / tau_dash            # WRONG

        # Overwrite collision for disc nodes
        # force_factor = (1.0 - 1.0/ (2.0 * tau_dash))   # WRONG: used tau_dash instead of tau_local
        force_factor = 1.0 - 0.5 * omega_local            # correct: (1 - Δt/(2τ̄)) = (1 - ω/2)
        w_edge = 1.0 / 36.0
        w_face = 1.0 / 18.0
        # Fi calculation Krüger (6.14) S.236
        # +x
        F_p00 = w_face * (3.0 * (1.0 - u_loc) + 9.0 * 1.0 * u_loc) * F_x_lat

        F_pp0 = w_edge * (3.0 * (1.0 - u_loc) + 9.0 * 1.0 * (u_loc + v_loc)) * F_x_lat
        F_pm0 = w_edge * (3.0 * (1.0 - u_loc) + 9.0 * 1.0 * (u_loc - v_loc)) * F_x_lat
        F_p0p = w_edge * (3.0 * (1.0 - u_loc) + 9.0 * 1.0 * (u_loc + w_loc)) * F_x_lat
        F_p0m = w_edge * (3.0 * (1.0 - u_loc) + 9.0 * 1.0 * (u_loc - w_loc)) * F_x_lat
        # -x
        F_m00 = w_face * (3.0 * (-1.0 - u_loc) + 9.0 * (-1.0) * (-u_loc)) * F_x_lat

        F_mp0 = w_edge * (3.0 * (-1.0 - u_loc) + 9.0 * (-1.0) * (-u_loc + v_loc)) * F_x_lat
        F_mm0 = w_edge * (3.0 * (-1.0 - u_loc) + 9.0 * (-1.0) * (-u_loc - v_loc)) * F_x_lat
        F_m0p = w_edge * (3.0 * (-1.0 - u_loc) + 9.0 * (-1.0) * (-u_loc + w_loc)) * F_x_lat
        F_m0m = w_edge * (3.0 * (-1.0 - u_loc) + 9.0 * (-1.0) * (-u_loc - w_loc)) * F_x_lat

        fS[Q000, x,   y,   z  ] = f000 - omega_local * (-feq000 + f000)
        fS[QP00, x+1, y,   z  ] = fp00 - omega_local * (-feqp00 + fp00) + force_factor * F_p00
        fS[QM00, x-1, y,   z  ] = fm00 - omega_local * (-feqm00 + fm00) + force_factor * F_m00
        fS[Q0M0, x,   y-1, z  ] = f0m0 - omega_local * (-feq0m0 + f0m0)
        fS[Q0P0, x,   y+1, z  ] = f0p0 - omega_local * (-feq0p0 + f0p0)
        fS[Q00M, x,   y,   z-1] = f00m - omega_local * (-feq00m + f00m)
        fS[Q00P, x,   y,   z+1] = f00p - omega_local * (-feq00p + f00p)
        fS[QPP0, x+1, y+1, z  ] = fpp0 - omega_local * (-feqpp0 + fpp0) + force_factor * F_pp0
        fS[QMP0, x-1, y+1, z  ] = fmp0 - omega_local * (-feqmp0 + fmp0) + force_factor * F_mp0
        fS[QPM0, x+1, y-1, z  ] = fpm0 - omega_local * (-feqpm0 + fpm0) + force_factor * F_pm0
        fS[QMM0, x-1, y-1, z  ] = fmm0 - omega_local * (-feqmm0 + fmm0) + force_factor * F_mm0
        fS[QP0P, x+1, y,   z+1] = fp0p - omega_local * (-feqp0p + fp0p) + force_factor * F_p0p
        fS[QM0P, x-1, y,   z+1] = fm0p - omega_local * (-feqm0p + fm0p) + force_factor * F_m0p
        fS[QP0M, x+1, y,   z-1] = fp0m - omega_local * (-feqp0m + fp0m) + force_factor * F_p0m
        fS[QM0M, x-1, y,   z-1] = fm0m - omega_local * (-feqm0m + fm0m) + force_factor * F_m0m
        fS[Q0MM, x,   y-1, z-1] = f0mm - omega_local * (-feq0mm + f0mm)
        fS[Q0MP, x,   y-1, z+1] = f0mp - omega_local * (-feq0mp + f0mp)
        fS[Q0PM, x,   y+1, z-1] = f0pm - omega_local * (-feq0pm + f0pm)
        fS[Q0PP, x,   y+1, z+1] = f0pp - omega_local * (-feq0pp + f0pp)
        
    end

    return nothing
end #collision_stream


end # module Kernel