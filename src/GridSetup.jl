module GridSetup
export add_rectangle!, find_object_boundary_nodes,
       _snap_fine_grid_position, _classify_coarse_nodes,
       _classify_fine_nodes, _compute_node_lists

function add_rectangle!(isObject, Nx, Ny, deltaX;
                        centerX, centerY, d, angleDeg,
                        originX=0.0, originY=0.0)
    half_len = 1.5 * d   # half of 3d  (long axis)
    half_hgt = 0.5 * d   # half of d   (short axis)

    alpha = deg2rad(angleDeg)
    c = cos(alpha)
    s = sin(alpha)

    for ix in 1:Nx, iy in 1:Ny
        x = originX + (ix - 2) * deltaX
        y = originY + (iy - 2) * deltaX

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
                                    half_len, half_hgt, cos_a, sin_a;
                                    originX=0.0, originY=0.0)
    dirs = ((1,0), (-1,0), (0,1), (0,-1), (1,1), (1,-1), (-1,1), (-1,-1))

    boundaryNodesAndDistances = Tuple{Int,Int,Int,Int,Float64}[]

    @inbounds for idx in fluidNodes
        ix, iy = Tuple(idx)
        dx0 = originX + (ix - 2) * deltaX - center_x
        dy0 = originY + (iy - 2) * deltaX - center_y
        lx0 =  cos_a * dx0 - sin_a * dy0
        ly0 =  sin_a * dx0 + cos_a * dy0

        for (cx, cy) in dirs
            if isObject[ix + cx, iy + cy]
                dx1 = originX + (ix + cx - 2) * deltaX - center_x
                dy1 = originY + (iy + cy - 2) * deltaX - center_y
                lx1 =  cos_a * dx1 - sin_a * dy1
                ly1 =  sin_a * dx1 + cos_a * dy1

                q = _ray_rect_q(lx0, ly0, lx1 - lx0, ly1 - ly0, half_len, half_hgt)
                push!(boundaryNodesAndDistances, (ix, iy, cx, cy, q))
            end
        end
    end

    return boundaryNodesAndDistances
end

function _snap_fine_grid_position(positionFineGridX, positionFineGridY, deltaX)
    snappedX = round(positionFineGridX / deltaX) * deltaX
    snappedY = round(positionFineGridY / deltaX) * deltaX
    if snappedX != positionFineGridX || snappedY != positionFineGridY
        @info "Fine grid anchor snapped to nearest coarse node: " *
              "($(positionFineGridX), $(positionFineGridY)) → ($(snappedX), $(snappedY))"
    end
    return snappedX, snappedY
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

    isFluid .&= .!isFineInterior

    return (; isInlet, isOutlet, isWall, isFluid, isObject, isSolid,
              isOuterInterfaceNode, isInnerInterfaceNode, isFineInterior,
              ix_left, ix_right, iy_bottom, iy_top)
end

function _classify_fine_nodes(NxFine, NyFine, deltaXFine,
                               positionFineGridX, positionFineGridY,
                               positionX, positionY, d, angleDeg)
    originXFine = positionFineGridX + 0.5 * deltaXFine
    originYFine = positionFineGridY + 0.5 * deltaXFine

    isFluidFine  = falses(NxFine, NyFine);  isFluidFine[2:NxFine-1, 2:NyFine-1] .= true
    isObjectFine = falses(NxFine, NyFine)
    add_rectangle!(isObjectFine, NxFine, NyFine, deltaXFine;
                   centerX=positionX, centerY=positionY, d=d, angleDeg=angleDeg,
                   originX=originXFine, originY=originYFine)
    isFluidFine .&= .!isObjectFine

    # Outer 2 rows/cols of fine interior → C→F target
    isOuterInterfaceNodeFine = falses(NxFine, NyFine)
    isOuterInterfaceNodeFine[2:NxFine-1, 2:NyFine-1] .= true
    isOuterInterfaceNodeFine[4:NxFine-3, 4:NyFine-3] .= false
    isOuterInterfaceNodeFine .&= isFluidFine

    # 4th and 5th rows/cols of fine interior → F→C source
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

end # module GridSetup
