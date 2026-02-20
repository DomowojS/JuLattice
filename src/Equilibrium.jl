module Equilibrium
export getEquilibrium

"""
    getEquilibrium(rho, u, v, cx, cy) -> Float64

Compute the D2Q9 equilibrium population for direction `(cx, cy)` at macroscopic
density `rho` and velocity `(u,v)`. Uses the product-form expressions consistent
with the C++ `setEq` implementation.
"""
function getEquilibrium(rho::Float64, u::Float64, v::Float64, cx::Int, cy::Int)
    # Precompute factors
    ux = u
    uy = v
    ux2 = u * u 
    uy2 = v * v

    if cx == -1 && cy == -1
        return rho * (1 - 3*ux + 3*ux2) * (1 - 3*uy + 3*uy2) / 36
    elseif cx ==  0 && cy == -1
        return -rho * (-2 + 3*ux2) * (1 + 3*uy2 - 3*uy) / 18
    elseif cx ==  1 && cy == -1
        return rho * (1 + 3*ux2 + 3*ux) * (1 + 3*uy2 - 3*uy) / 36
    elseif cx == -1 && cy ==  0
        return -rho * (-2 + 3*uy2) * (1 + 3*ux2 - 3*ux) / 18
    elseif cx ==  0 && cy ==  0
        return rho * (-2 + 3*ux2) * (-2 + 3*uy2) / 9
    elseif cx ==  1 && cy ==  0
        return -rho * (-2 + 3*uy2) * (1 + 3*ux2 + 3*ux) / 18
    elseif cx == -1 && cy ==  1
        return rho * (1 + 3*ux2 - 3*ux) * (1 + 3*uy2 + 3*uy) / 36
    elseif cx ==  0 && cy ==  1
        return -rho * (-2 + 3*ux2) * (1 + 3*uy2 + 3*uy) / 18
    elseif cx ==  1 && cy ==  1
        return rho * (1 + 3*ux2 + 3*ux) * (1 + 3*uy2 + 3*uy) / 36
    else
        error("Invalid lattice direction (cx,cy)=($cx,$cy)")
    end
end

end # module