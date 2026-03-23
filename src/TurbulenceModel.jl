module TurbulenceModel
    
export compute_pi_norm, smagorinsky_omega

function compute_pi_norm(
    fneqm00, fneqp00, fneq0m0, fneq0p0, fneq00m, fneq00p,
    fneqmm0, fneqmp0, fneqpm0, fneqpp0,
    fneqm0m, fneqm0p, fneqp0m, fneqp0p,
    fneq0mm, fneq0mp, fneq0pm, fneq0pp
    )

    # Π = √(π_ij^neq π_ij^neq)
    # π_ij^neq = ∑_a e_ai e_aj (f_a − f_a^eq)  

    # π_xx
    pi_xx = fneqm00 + fneqp00 +
            fneqmm0 + fneqmp0 + fneqpm0 + fneqpp0 +
            fneqm0m + fneqm0p + fneqp0m + fneqp0p
    
    # π_yy
    pi_yy = fneq0m0 + fneq0p0 +
            fneqmm0 + fneqmp0 + fneqpm0 + fneqpp0 +
            fneq0mm + fneq0mp + fneq0pm + fneq0pp
    
    # π_zz
    pi_zz = fneq00m + fneq00p +
            fneqm0m + fneqm0p + fneqp0m + fneqp0p +
            fneq0mm + fneq0mp + fneq0pm + fneq0pp      
    
    # π_xy, π_xz, π_yz       
    pi_xy =  fneqmm0 - fneqmp0 - fneqpm0 + fneqpp0
    pi_xz =  fneqm0m - fneqm0p - fneqp0m + fneqp0p
    pi_yz =  fneq0mm - fneq0mp - fneq0pm + fneq0pp

    # Π^2 = ... -> √Π² = Π ensures positive value
    pi_norm_sq = (pi_xx*pi_xx) + (pi_yy*pi_yy) + (pi_zz*pi_zz) + 
                    2 * ((pi_xy*pi_xy) + (pi_xz*pi_xz) + (pi_yz*pi_yz))
    
    return sqrt(max(pi_norm_sq, 0.0))
end

function smagorinsky_omega(tau, CS, pi_neq_norm, rho)
    # compute τ_eff directly from Π^neq
    tau_eff = 0.5 * (tau + sqrt((tau*tau) + 18 * (CS*CS) * pi_neq_norm / rho))

    #return ω_eff
    return 1.0 / tau_eff
end

end