module GridSetup
export add_rectangle!, find_object_boundary_nodes, _classify_nodes

function add_rectangle!(isObject, Nx, Ny, deltaX;
                        centerX, centerY, d, angleDeg)
    half_len = 1.5 * d
    half_hgt = 0.5 * d

    alpha = deg2rad(angleDeg)
    c = cos(alpha)
    s = sin(alpha)

    for ix in 1:Nx, iy in 1:Ny
        x = (ix - 2) * deltaX
        y = (iy - 2) * deltaX

        dx = x - centerX
        dy = y - centerY

        lx = c * dx - s * dy
        ly = s * dx + c * dy

        if abs(lx) <= half_len && abs(ly) <= half_hgt
            isObject[ix, iy] = true
        end
    end
end

function _ray_rect_q(lx0, ly0, dlx, dly, half_len, half_hgt)
    q = Inf
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
        dx0 = (ix - 2) * deltaX - center_x
        dy0 = (iy - 2) * deltaX - center_y
        lx0 =  cos_a * dx0 - sin_a * dy0
        ly0 =  sin_a * dx0 + cos_a * dy0

        for (cx, cy) in dirs
            if isObject[ix + cx, iy + cy]
                dx1 = (ix + cx - 2) * deltaX - center_x
                dy1 = (iy + cy - 2) * deltaX - center_y
                lx1 =  cos_a * dx1 - sin_a * dy1
                ly1 =  sin_a * dx1 + cos_a * dy1

                q = _ray_rect_q(lx0, ly0, lx1 - lx0, ly1 - ly0, half_len, half_hgt)
                push!(boundaryNodesAndDistances, (ix, iy, cx, cy, q))
            end
        end
    end

    return boundaryNodesAndDistances
end

function _classify_nodes(Nx, Ny, deltaX, positionX, positionY, d, angleDeg)
    isInlet  = falses(Nx, Ny);  isInlet[1,      :]        .= true
    isOutlet = falses(Nx, Ny);  isOutlet[Nx,    :]        .= true
    isWall   = falses(Nx, Ny);  isWall[2:Nx-1, [1, Ny]]   .= true
    isFluid  = falses(Nx, Ny);  isFluid[2:Nx-1,  2:Ny-1]  .= true
    isObject = falses(Nx, Ny)
    add_rectangle!(isObject, Nx, Ny, deltaX;
                   centerX=positionX, centerY=positionY, d=d, angleDeg=angleDeg)
    isFluid .&= .!isObject
    isSolid  = isInlet .| isOutlet .| isWall .| isObject

    fluidNodes  = findall(isFluid)
    solidNodes  = findall(.!isFluid .& .!isInlet .& .!isOutlet .& .!isWall)
    objectNodes = findall(isObject)

    boundaryNodesAndDistances = find_object_boundary_nodes(
        isObject, fluidNodes,
        deltaX, positionX, positionY,
        1.5*d, 0.5*d, cosd(angleDeg), sind(angleDeg))

    return (; isInlet, isOutlet, isWall, isFluid, isObject, isSolid,
              fluidNodes, solidNodes, objectNodes, boundaryNodesAndDistances)
end

end # module GridSetup
