module Kernel

using ..TurbulenceModel: compute_pi_norm, smagorinsky_omega

export  collision_stream!

function collision_stream!(
    gridlengthX::Int, gridlengthY::Int, gridlengthZ::Int,
    τ::Float64, CS::Float64, is_fluid,
    rho, u, v, w,
    f000,  fm00,  fp00,  f0m0,  f0p0,  f00m,  f00p,
    fmm0,  fmp0,  fpm0,  fpp0,
    fm0m,  fm0p,  fp0m,  fp0p,
    f0mm,  f0mp,  f0pm,  f0pp,
    f000S, fm00S, fp00S, f0m0S, f0p0S, f00mS, f00pS,
    fmm0S, fmp0S, fpm0S, fpp0S,
    fm0mS, fm0pS, fp0mS, fp0pS,
    f0mmS, f0mpS, f0pmS, f0ppS    
)
    


    # Iterate over all cells except boundary cells
    @inbounds Threads.@threads for z in 2:gridlengthZ-1
        for y in 2:gridlengthY-1
            for x in 2:gridlengthX-1

                if !is_fluid[x,y,z]
                    continue
                end


                # Compute macroscopic quantities
                rho_loc = f000[x,y,z] + 
                            (fm00[x,y,z] + fp00[x,y,z] + f0m0[x,y,z] + f0p0[x,y,z] + f00m[x,y,z] + f00p[x,y,z]) +
                            (fmm0[x,y,z] + fmp0[x,y,z] + fpm0[x,y,z] + fpp0[x,y,z] + 
                            fm0m[x,y,z] + fm0p[x,y,z] + fp0m[x,y,z] + fp0p[x,y,z] +
                            f0mm[x,y,z] + f0mp[x,y,z] + f0pm[x,y,z] + f0pp[x,y,z])

                rho[x,y,z] = rho_loc
                inv_rho = 1.0 / rho_loc
                
                u_loc = ((-fm00[x,y,z] + fp00[x,y,z]) +
                            (-fmm0[x,y,z] - fmp0[x,y,z] + fpm0[x,y,z] + fpp0[x,y,z]) +
                            (-fm0m[x,y,z] - fm0p[x,y,z] + fp0m[x,y,z] + fp0p[x,y,z])) * inv_rho
                
                v_loc = ((-f0m0[x,y,z] + f0p0[x,y,z]) +
                            (-fmm0[x,y,z] + fmp0[x,y,z] - fpm0[x,y,z] + fpp0[x,y,z]) +
                            (-f0mm[x,y,z] - f0mp[x,y,z] + f0pm[x,y,z] + f0pp[x,y,z])) * inv_rho
                
                w_loc = ((-f00m[x,y,z] + f00p[x,y,z]) +
                            (-fm0m[x,y,z] + fm0p[x,y,z] - fp0m[x,y,z] + fp0p[x,y,z]) +
                            (-f0mm[x,y,z] + f0mp[x,y,z] - f0pm[x,y,z] + f0pp[x,y,z])) * inv_rho

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
                # fneq000 = f000[x,y,z] - feq000
                
                fneqm00 = fm00[x,y,z] - feqm00
                fneqp00 = fp00[x,y,z] - feqp00
                fneq0m0 = f0m0[x,y,z] - feq0m0
                fneq0p0 = f0p0[x,y,z] - feq0p0
                fneq00m = f00m[x,y,z] - feq00m
                fneq00p = f00p[x,y,z] - feq00p
                fneqmm0 = fmm0[x,y,z] - feqmm0
                fneqmp0 = fmp0[x,y,z] - feqmp0
                fneqpm0 = fpm0[x,y,z] - feqpm0
                fneqpp0 = fpp0[x,y,z] - feqpp0
                fneqm0m = fm0m[x,y,z] - feqm0m
                fneqm0p = fm0p[x,y,z] - feqm0p
                fneqp0m = fp0m[x,y,z] - feqp0m
                fneqp0p = fp0p[x,y,z] - feqp0p
                fneq0mm = f0mm[x,y,z] - feq0mm
                fneq0mp = f0mp[x,y,z] - feq0mp
                fneq0pm = f0pm[x,y,z] - feq0pm
                fneq0pp = f0pp[x,y,z] - feq0pp

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
                f000S[x,y,z] = f000[x,y,z] + omega_local * (feq000 - f000[x,y,z])
                # Face neighbors
                fm00S[x-1,y,z] = fm00[x,y,z] + omega_local * (feqm00 - fm00[x,y,z])
                fp00S[x+1,y,z] = fp00[x,y,z] + omega_local * (feqp00 - fp00[x,y,z])
                f0m0S[x,y-1,z] = f0m0[x,y,z] + omega_local * (feq0m0 - f0m0[x,y,z])
                f0p0S[x,y+1,z] = f0p0[x,y,z] + omega_local * (feq0p0 - f0p0[x,y,z])
                f00mS[x,y,z-1] = f00m[x,y,z] + omega_local * (feq00m - f00m[x,y,z])
                f00pS[x,y,z+1] = f00p[x,y,z] + omega_local * (feq00p - f00p[x,y,z])
                # XY-plane edges
                fmm0S[x-1,y-1,z] = fmm0[x,y,z] + omega_local * (feqmm0 - fmm0[x,y,z])
                fmp0S[x-1,y+1,z] = fmp0[x,y,z] + omega_local * (feqmp0 - fmp0[x,y,z])
                fpm0S[x+1,y-1,z] = fpm0[x,y,z] + omega_local * (feqpm0 - fpm0[x,y,z])
                fpp0S[x+1,y+1,z] = fpp0[x,y,z] + omega_local * (feqpp0 - fpp0[x,y,z])
                # XZ-plane edges
                fm0mS[x-1,y,z-1] = fm0m[x,y,z] + omega_local * (feqm0m - fm0m[x,y,z])
                fm0pS[x-1,y,z+1] = fm0p[x,y,z] + omega_local * (feqm0p - fm0p[x,y,z])
                fp0mS[x+1,y,z-1] = fp0m[x,y,z] + omega_local * (feqp0m - fp0m[x,y,z])
                fp0pS[x+1,y,z+1] = fp0p[x,y,z] + omega_local * (feqp0p - fp0p[x,y,z])
                # YZ-plane edges
                f0mmS[x,y-1,z-1] = f0mm[x,y,z] + omega_local * (feq0mm - f0mm[x,y,z])
                f0mpS[x,y-1,z+1] = f0mp[x,y,z] + omega_local * (feq0mp - f0mp[x,y,z])
                f0pmS[x,y+1,z-1] = f0pm[x,y,z] + omega_local * (feq0pm - f0pm[x,y,z])
                f0ppS[x,y+1,z+1] = f0pp[x,y,z] + omega_local * (feq0pp - f0pp[x,y,z])

            end #end x
        end #end y
    end #end z

    return nothing
end #collision_stream


end # module Kernel