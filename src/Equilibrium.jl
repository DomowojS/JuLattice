module Equilibrium
export getEquilibrium

# Distribution ordering: [f00, fp0, fm0, f0p, f0m, fpp, fpm, fmp, fmm]
# Moment ordering:       [m00, m10, m01, m11, mP,  mxx, m21, m12, m22]
const _M = Float64[
    1      1      1      1      1      1      1      1      1   ;  # m00
    0      1     -1      0      0      1      1     -1     -1   ;  # m10
    0      0      0      1     -1      1     -1      1     -1   ;  # m01
    0      0      0      0      0      1     -1     -1      1   ;  # m11
   -2/3    1/3    1/3    1/3    1/3    4/3    4/3    4/3    4/3 ;  # mP  = (m20+m02) - 2cs²·m00
    0      1      1     -1     -1      0      0      0      0   ;  # mxx = m20-m02
    0      0      0     -1/3   1/3    2/3   -2/3    2/3   -2/3  ;  # m21 = m21 - cs²·m01
    0     -1/3   1/3     0      0     2/3    2/3   -2/3   -2/3  ;  # m12 = m12 - cs²·m10
   1/9   -2/9   -2/9   -2/9   -2/9   4/9    4/9    4/9    4/9  ]  # m22 = m22 - cs²·(m20+m02) + cs⁴·m00

const _M_inv = inv(_M)

# Map (cx, cy) -> population index in the distribution ordering above
const _DIR_IDX = Dict(
    ( 0,  0) => 1,   # f00
    ( 1,  0) => 2,   # fp0
    (-1,  0) => 3,   # fm0
    ( 0,  1) => 4,   # f0p
    ( 0, -1) => 5,   # f0m
    ( 1,  1) => 6,   # fpp
    ( 1, -1) => 7,   # fpm
    (-1,  1) => 8,   # fmp
    (-1, -1) => 9    # fmm
)

"""
    getEquilibrium(rho, u, v, cx, cy) -> Float64

Returns the D2Q9 equilibrium population for lattice direction (cx, cy) ∈ {-1,0,+1}²
at macroscopic state (rho, u, v).

Computes via equilibrium moment vector and M⁻¹ transform:
    meq = [m00_eq, m10_eq, ..., m22_eq]
    feq = M⁻¹ * meq
"""
function getEquilibrium(rho::Float64, u::Float64, v::Float64, cx::Int, cy::Int)
    meq = [
        rho,
        rho * u,
        rho * v,
        rho * u * v,
        rho * (u^2 + v^2),     # mP  = (m20+m02) - 2cs²·ρ
        rho * (u^2 - v^2),     # mxx unchanged
        rho * u^2 * v,         # m21 = m21_raw - cs²·m01
        rho * u * v^2,         # m12 = m12_raw - cs²·m10
        rho * u^2 * v^2        # m22 = m22_raw - cs²·(m20+m02) + cs⁴·m00
    ]
    feq = _M_inv * meq
    return feq[_DIR_IDX[(cx, cy)]]
end

end#Equilibrium
